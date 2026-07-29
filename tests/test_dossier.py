import json
from hashlib import sha256
from io import BytesIO
from pathlib import Path

import pytest
import xlwt

from geneimpact.crispritz import CRISPRITZ_COMMIT
from geneimpact.dossier import build_research_dossier, verify_dossier_integrity
from geneimpact.housden import HOUSDEN_SERVICE_URL


CRISPRSCAN_CONTEXT = "ACCTGGATCGATGCTGATGCTAGATAAGGTTGAGC"
CRISPRITZ_HEADER = (
    "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
    "Direction\tMismatches\tBulge Size\tTotal\n"
)
CRISPRITZ_ROW = (
    "X\tAAAAAAAAAAAAAAAAAAAANNN\tAAAAAAAAAAAAAAAAAAAATGG\t"
    "chr1\t10\t10\t+\t0\t0\t0\n"
)


def _request(**context_overrides):
    context = {
        "study_id": "zf-study-001",
        "species_profile": "zebrafish",
        "strain_or_breed": "Tuebingen",
        "genome_build": "GRCz12tu",
        "assembly_accession": "GCF_049306965.2",
        "edit_class": "knockout",
        "delivery_context": "t7_in_vitro_transcription",
        "developmental_context": "zebrafish_embryo",
        "evidence_snapshot_sha256": "b" * 64,
    }
    context.update(context_overrides)
    return {
        "dossier_schema_version": "1.0",
        "study_context": context,
        "target_genes": [
            {
                "gene_id": "ENSDARG00000000001",
                "gene_symbol": "geneA",
                "intended_change": "loss_of_function",
                "evidence_signal": 0.8,
                "evidence_reference": "curated-zebrafish-record-1",
            },
            {
                "gene_id": "ENSDARG00000000002",
                "gene_symbol": "geneB",
                "intended_change": "loss_of_function",
                "evidence_signal": 0.6,
                "evidence_reference": "curated-zebrafish-record-2",
            },
        ],
        "intended_outcomes": ["Reduced activity of the declared target pathway"],
        "welfare_endpoints": ["Survival", "Developmental morphology"],
        "interaction_evidence": [
            {
                "genes": [
                    "ENSDARG00000000001",
                    "ENSDARG00000000002",
                ],
                "evidence_weight": 0.5,
                "evidence_reference": "curated-pair-record-1",
            }
        ],
        "evidence": {
            "on_target_uncertainty": 0.3,
            "off_target_evidence": 0.4,
            "network_impact_evidence": 0.5,
            "welfare_relevance": 0.6,
        },
        "predictors": {},
    }


def _housden_result():
    response = _housden_service_response()
    return {
        "request": {
            "guide_id": "drsc-guide-01",
            "species_profile": "fruit_fly",
            "genome_build": "Release 6 plus ISO1 MT",
            "assembly_accession": "GCF_000001215.4",
            "sequence_source_strain_or_isolate": "ISO-1",
            "protospacer": "ATCTGACCTCCCGGCTAATT",
            "nuclease": "SpCas9",
            "guide_expression": "u6_sgrna",
            "developmental_context": "drosophila_s2r_plus_cell_culture",
        },
        "execution": {
            "source": "flyrnai_evaluate_crispr",
            "source_url": HOUSDEN_SERVICE_URL,
            "retrieved_at": "2026-07-29T09:30:00+00:00",
            "source_response_sha256": sha256(response).hexdigest(),
        },
        "raw_output": {"housden_score": 6.93243},
    }


def _housden_service_response():
    workbook = xlwt.Workbook()
    sheet = workbook.add_sheet("EfficiencyScoreList")
    headers = (
        "Identifier",
        "Input Sequence",
        "Input Sequence Length",
        "Analyzed Sequence",
        "Analyzed Sequence Length",
        "Number of Bases Analyzed",
        "Start Index",
        "End Index",
        "Score",
        "Comments",
        "U6 Terminator",
    )
    values = (
        "1",
        "ATCTGACCTCCCGGCTAATT",
        20,
        "ATCTGACCTCCCGGCTAATT",
        20,
        20,
        1,
        20,
        "6.93243",
        "None",
        "None",
    )
    for column, value in enumerate(headers):
        sheet.write(0, column, value)
    for column, value in enumerate(values):
        sheet.write(1, column, value)
    buffer = BytesIO()
    workbook.save(buffer)
    return buffer.getvalue()


def test_builds_integrity_checkable_multigene_dossier(tmp_path):
    request = _request()
    request["predictors"]["crisprscan"] = {
        "nuclease": "SpCas9",
        "guide_expression": "t7_in_vitro_transcription",
        "developmental_context": "zebrafish_embryo",
        "guides": [
            {
                "guide_id": "zf-guide-001",
                "context_35nt": CRISPRSCAN_CONTEXT,
            }
        ],
    }

    dossier = build_research_dossier(
        request,
        attachment_base_dir=tmp_path,
        source_request_sha256="c" * 64,
    )

    assert dossier["dossier_schema_version"] == "1.0"
    assert dossier["study"]["species_profile"] == "zebrafish"
    assert dossier["study"]["assembly_accession"] == "GCF_049306965.2"
    assert dossier["declared_evidence_inputs"]["welfare_relevance"] == 0.6
    assert dossier["interaction_hypotheses"]["evidence_supported_pairs"] == 1
    assert dossier["interaction_hypotheses"]["declared_pair_evidence"][0][
        "evidence_reference"
    ] == "curated-pair-record-1"
    assert dossier["interaction_hypotheses"]["ranked_pairs"][0]["priority"] > 0
    assert dossier["model_predictions"][0]["predictor"] == "CRISPRscan"
    assert CRISPRSCAN_CONTEXT not in json.dumps(dossier)
    assert verify_dossier_integrity(dossier).matches
    assert dossier["integrity"]["signature_status"] == "unsigned"


def test_integrity_verifier_detects_report_modification(tmp_path):
    dossier = build_research_dossier(
        _request(),
        attachment_base_dir=tmp_path,
    )
    dossier["study"]["study_id"] = "tampered"

    result = verify_dossier_integrity(dossier)

    assert not result.matches
    assert result.expected_content_sha256 != result.actual_content_sha256


def test_integrates_crispritz_attachment_without_escaping_request_dir(tmp_path):
    targets = tmp_path / "rat.targets.txt"
    targets.write_text(CRISPRITZ_HEADER + CRISPRITZ_ROW, encoding="utf-8")
    request = _request(
        study_id="rat-study-001",
        species_profile="rat",
        strain_or_breed="BN/NHsdMcwi",
        genome_build="GRCr8",
        assembly_accession="GCF_036323735.1",
        delivery_context="spcas9_ribonucleoprotein",
        developmental_context="one_cell_embryo",
    )
    request["predictors"]["crispritz"] = {
        "targets_file": "rat.targets.txt",
        "pam_definition": "NNNNNNNNNNNNNNNNNNNNNGG 3",
        "max_mismatches": 4,
        "max_dna_bulge": 1,
        "max_rna_bulge": 1,
        "crispritz_commit": CRISPRITZ_COMMIT,
        "reference_fasta_sha256": "a" * 64,
        "variant_aware": False,
        "variant_snapshot_sha256": None,
    }

    dossier = build_research_dossier(request, attachment_base_dir=tmp_path)

    prediction = dossier["model_predictions"][0]
    assert prediction["predictor"] == "CRISPRitz"
    assert prediction["species_profile"] == "rat"
    assert prediction["candidate_site_count"] == 1
    assert any(
        row["predictor"] == "CRISPRitz"
        and row["execution_state"] == "included"
        for row in dossier["capability_coverage"]
    )


def test_rejects_attachment_path_escape(tmp_path):
    request = _request()
    request["predictors"]["crispritz"] = {"targets_file": "../outside.txt"}

    with pytest.raises(ValueError, match="escapes"):
        build_research_dossier(request, attachment_base_dir=tmp_path)


def test_integrates_multiple_version_locked_indelphi_results(tmp_path):
    result = Path("examples/indelphi-mouse-result.json").read_text(
        encoding="utf-8"
    )
    (tmp_path / "guide-1.json").write_text(result, encoding="utf-8")
    second = json.loads(result)
    second["request"]["target_id"] = "Tyr-guide-synthetic-02"
    (tmp_path / "guide-2.json").write_text(
        json.dumps(second),
        encoding="utf-8",
    )
    request = _request(
        study_id="mouse-indelphi-study",
        species_profile="mouse",
        strain_or_breed="C57BL/6J",
        genome_build="GRCm39",
        assembly_accession="GCF_000001635.27",
        delivery_context="spcas9_ribonucleoprotein",
        developmental_context="one_cell_embryo",
    )
    request["predictors"]["indelphi"] = {
        "result_files": ["guide-1.json", "guide-2.json"]
    }

    dossier = build_research_dossier(request, attachment_base_dir=tmp_path)

    predictions = [
        row
        for row in dossier["model_predictions"]
        if row["predictor"] == "inDelphi"
    ]
    assert len(predictions) == 2
    assert predictions[0]["source_document_sha256"]
    assert predictions[0]["external_validation"]["guide_count"] == 14
    assert any(
        row["predictor"] == "inDelphi"
        and row["execution_state"] == "included"
        for row in dossier["capability_coverage"]
    )
    assert "inDelphi" not in dossier["evidence_completeness"][
        "available_predictors_not_run"
    ]


def test_rejects_indelphi_context_mismatch(tmp_path):
    result = json.loads(
        Path("examples/indelphi-mouse-result.json").read_text(encoding="utf-8")
    )
    result["request"]["developmental_context"] = "mismatch"
    (tmp_path / "guide.json").write_text(json.dumps(result), encoding="utf-8")
    request = _request(
        species_profile="mouse",
        strain_or_breed="C57BL/6J",
        genome_build="GRCm39",
        assembly_accession="GCF_000001635.27",
        delivery_context="spcas9_ribonucleoprotein",
        developmental_context="one_cell_embryo",
    )
    request["predictors"]["indelphi"] = {"result_files": ["guide.json"]}

    with pytest.raises(ValueError, match="developmental_context"):
        build_research_dossier(request, attachment_base_dir=tmp_path)


def test_integrates_housden_only_for_matching_fruit_fly_cell_context(tmp_path):
    result_path = tmp_path / "housden.json"
    response_path = tmp_path / "housden.xls"
    result_path.write_text(json.dumps(_housden_result()), encoding="utf-8")
    response_path.write_bytes(_housden_service_response())
    request = _request(
        study_id="fruit-fly-cell-study",
        species_profile="fruit_fly",
        strain_or_breed="ISO-1",
        genome_build="Release 6 plus ISO1 MT",
        assembly_accession="GCF_000001215.4",
        delivery_context="u6_sgrna",
        developmental_context="drosophila_s2r_plus_cell_culture",
    )
    request["predictors"]["housden"] = {
        "result_files": ["housden.json"],
        "source_response_files": ["housden.xls"],
    }

    dossier = build_research_dossier(request, attachment_base_dir=tmp_path)

    prediction = next(
        row
        for row in dossier["model_predictions"]
        if row["predictor"] == "Housden"
    )
    assert prediction["housden_score"] == 6.93243
    assert "ATCTGACCTCCCGGCTAATT" not in json.dumps(dossier)
    assert any(
        row["predictor"] == "Housden"
        and row["execution_state"] == "included"
        for row in dossier["capability_coverage"]
    )


def test_rejects_housden_result_for_fruit_fly_embryo_study(tmp_path):
    result_path = tmp_path / "housden.json"
    response_path = tmp_path / "housden.xls"
    result_path.write_text(json.dumps(_housden_result()), encoding="utf-8")
    response_path.write_bytes(_housden_service_response())
    request = _request(
        study_id="fruit-fly-embryo-study",
        species_profile="fruit_fly",
        strain_or_breed="ISO-1",
        genome_build="Release 6 plus ISO1 MT",
        assembly_accession="GCF_000001215.4",
        delivery_context="u6_sgrna",
        developmental_context="fruit_fly_embryo",
    )
    request["predictors"]["housden"] = {
        "result_files": ["housden.json"],
        "source_response_files": ["housden.xls"],
    }

    with pytest.raises(ValueError, match="developmental_context"):
        build_research_dossier(request, attachment_base_dir=tmp_path)


def test_generic_external_output_cannot_impersonate_dedicated_housden_adapter(
    tmp_path,
):
    request = _request(
        study_id="fruit-fly-cell-study",
        species_profile="fruit_fly",
        strain_or_breed="ISO-1",
        genome_build="Release 6 plus ISO1 MT",
        assembly_accession="GCF_000001215.4",
        delivery_context="u6_sgrna",
        developmental_context="drosophila_s2r_plus_cell_culture",
    )
    request["predictors"]["external_concern_outputs"] = [
        {
            "predictor": "Housden",
            "predictor_version": "self-declared",
            "task": "guide_activity",
            "concern_score": 0.1,
            "confidence": 1.0,
            "supported_species": ["fruit_fly"],
            "supported_edit_classes": ["knockout"],
            "evidence_reference": "self-declared",
        }
    ]

    dossier = build_research_dossier(request, attachment_base_dir=tmp_path)

    housden = next(
        row
        for row in dossier["capability_coverage"]
        if row["predictor"] == "Housden"
    )
    assert housden["execution_state"] == "available_not_run"


def test_rejects_predictor_context_inconsistent_with_study(tmp_path):
    request = _request()
    request["predictors"]["crisprscan"] = {
        "nuclease": "SpCas9",
        "guide_expression": "u6",
        "developmental_context": "zebrafish_embryo",
        "guides": [],
    }

    with pytest.raises(ValueError, match="must match"):
        build_research_dossier(request, attachment_base_dir=tmp_path)


@pytest.mark.parametrize(
    (
        "species_profile",
        "strain",
        "genome_build",
        "assembly_accession",
    ),
    [
        ("mouse", "C57BL/6J", "GRCm39", "GCF_000001635.27"),
        ("rat", "BN/NHsdMcwi", "GRCr8", "GCF_036323735.1"),
        ("zebrafish", "Tuebingen", "GRCz12tu", "GCF_049306965.2"),
        (
            "fruit_fly",
            "ISO-1",
            "Release 6 plus ISO1 MT",
            "GCF_000001215.4",
        ),
        (
            "rhesus_macaque",
            "MMU2019108-1",
            "T2T-MMU8v2.0",
            "GCF_049350105.2",
        ),
        (
            "cynomolgus_macaque",
            "582-1",
            "T2T-MFA8v1.1",
            "GCF_037993035.2",
        ),
    ],
)
def test_builds_context_bound_dossier_for_every_registered_species(
    tmp_path,
    species_profile,
    strain,
    genome_build,
    assembly_accession,
):
    request = _request(
        study_id=f"{species_profile}-study",
        species_profile=species_profile,
        strain_or_breed=strain,
        genome_build=genome_build,
        assembly_accession=assembly_accession,
        delivery_context="declared_delivery",
        developmental_context="declared_developmental_context",
    )

    dossier = build_research_dossier(request, attachment_base_dir=tmp_path)

    assert dossier["study"]["species_profile"] == species_profile
    assert dossier["study"]["assembly_accession"] == assembly_accession
    assert verify_dossier_integrity(dossier).matches
    assert any(
        row["predictor"] == "CRISPRitz"
        and row["execution_state"] == "available_not_run"
        for row in dossier["capability_coverage"]
    )


def test_rejects_generic_monkey_profile(tmp_path):
    request = _request(species_profile="monkey")

    with pytest.raises(ValueError, match="exact registered key"):
        build_research_dossier(request, attachment_base_dir=tmp_path)
