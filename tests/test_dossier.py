import json

import pytest

from geneimpact.crispritz import CRISPRITZ_COMMIT
from geneimpact.dossier import build_research_dossier, verify_dossier_integrity


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
