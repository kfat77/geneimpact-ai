from dataclasses import asdict, replace
from hashlib import sha256
from io import BytesIO
import json
import sys

from openpyxl import Workbook
import pytest

import geneimpact.cli as cli
from geneimpact.cli import main
from geneimpact.cynomolgus_base_editing import (
    CynomolgusBaseEditingBlock,
    CynomolgusBaseEditingSource,
    evaluate_cynomolgus_base_editing_transfer,
    prepare_cynomolgus_base_editing_transfer_template,
)


TARGET_SEQUENCES = {
    "GENE1-S1": "ACGTACGTACGTACGTACGTTGG",
    "GENE2-S1": "TGCATGCATGCATGCATGCACGG",
}


def _workbook_bytes(workbook):
    stream = BytesIO()
    workbook.save(stream)
    return stream.getvalue()


def _source_workbooks():
    targets = Workbook()
    target_sheet = targets.active
    target_sheet.title = "sheet1"
    target_sheet.append(["Supplementary Data 1 Off-target sites"])
    target_sheet.append(
        [
            "Target site",
            "Gene",
            "Type",
            "Sequence",
            "Chromosome",
            "Position",
            "Direction",
        ]
    )
    target_sheet.append(
        [
            "GENE1-S1-On",
            "GENE1",
            "Exon",
            TARGET_SEQUENCES["GENE1-S1"],
            "Chr1",
            "10-32",
            "+",
        ]
    )
    target_sheet.append(
        [
            "GENE1-S1-OT1",
            "intergenic",
            "-",
            "ACGTACGTACGTACGTACGTCGG",
            "Chr2",
            "20-42",
            "-",
        ]
    )
    target_sheet.append(
        [
            "GENE2-S1-On",
            "GENE2",
            "Exon",
            TARGET_SEQUENCES["GENE2-S1"],
            "Chr3",
            "30-52",
            "+",
        ]
    )

    source_data = Workbook()
    source_sheet = source_data.active
    source_sheet.title = "Genotyping by Sanger Sequencing"
    source_sheet.append(["Source Data for genotyping by Sanger sequencing"])
    source_sheet.append(["GENE1-S1", None, "C4", None, "C6"])
    source_sheet.append(
        [
            "Embryo ID",
            "Number of Total clone",
            "Number of C>T clone",
            "Number of non C>T clone",
            "Number of C>T clone",
        ]
    )
    source_sheet.append([1, 10, "5 (50%)", "0 (0%)", "2 (20%)"])
    source_sheet.append([2, 20, "10 (50%)", "0 (0%)", "8 (40%)"])
    source_sheet.append([])
    source_sheet.append(["GENE2-S1", None, "A5"])
    source_sheet.append(
        [
            "Embryo ID",
            "Number of Total clone",
            "Number of A>G clone",
        ]
    )
    source_sheet.append([1, 8, "4 (50%)"])
    source_sheet.append([2, 12, "6 (50%)"])
    return _workbook_bytes(targets), _workbook_bytes(source_data)


def _source_profile(targets, source_data):
    return CynomolgusBaseEditingSource(
        source_id="synthetic-cynomolgus-base-editing",
        article_reference="https://example.test/article",
        target_sites_url="https://example.test/targets.xlsx",
        source_data_url="https://example.test/source-data.xlsx",
        target_sites_sha256=sha256(targets).hexdigest(),
        source_data_sha256=sha256(source_data).hexdigest(),
        source_genome_build="synthetic-build",
        expected_candidate_site_count=3,
        expected_on_target_site_count=2,
        expected_record_count=3,
        expected_context_count=2,
        expected_embryo_base_observation_count=6,
        expected_clone_denominator_count=80,
        blocks=(
            CynomolgusBaseEditingBlock(
                context_id="single_be3_gene1",
                editor="BE3",
                multiplex_guide_count=1,
                target_site_id="GENE1-S1",
                gene="GENE1",
                conversion="C_to_T",
                target_bases=("C4", "C6"),
                label_row=2,
                header_row=3,
                first_data_row=4,
                last_data_row=5,
                intended_count_columns=(3, 5),
            ),
            CynomolgusBaseEditingBlock(
                context_id="single_abe710_gene2",
                editor="ABE7.10",
                multiplex_guide_count=1,
                target_site_id="GENE2-S1",
                gene="GENE2",
                conversion="A_to_G",
                target_bases=("A5",),
                label_row=7,
                header_row=8,
                first_data_row=9,
                last_data_row=10,
                intended_count_columns=(3,),
            ),
        ),
    )


def _write_sources(tmp_path):
    targets, source_data = _source_workbooks()
    targets_path = tmp_path / "targets.xlsx"
    source_data_path = tmp_path / "source-data.xlsx"
    targets_path.write_bytes(targets)
    source_data_path.write_bytes(source_data)
    return (
        targets_path,
        source_data_path,
        _source_profile(targets, source_data),
    )


def test_prepares_sequence_redacted_cynomolgus_prediction_template(tmp_path):
    targets, source_data, source = _write_sources(tmp_path)

    template = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=source,
    )

    assert template["schema_version"] == (
        "geneimpact.cynomolgus_base_editing_transfer_predictions.v1"
    )
    assert len(template["records"]) == 3
    assert template["records"][0]["editor"] == "BE3"
    assert template["records"][0]["target_base"] == "C4"
    assert template["records"][0]["target_sequence_sha256"] == sha256(
        TARGET_SEQUENCES["GENE1-S1"].encode("ascii")
    ).hexdigest()
    assert template["records"][0]["predicted_score"] is None
    assert not any(sequence in repr(template) for sequence in TARGET_SEQUENCES.values())


def test_evaluates_only_within_context_and_hides_source_labels(tmp_path):
    targets, source_data, source = _write_sources(tmp_path)
    predictions = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=source,
    )
    predictions["prediction"].update(
        {
            "name": "synthetic-oracle",
            "version": "test-v1",
            "submitted_code_revision": "test-commit",
            "score_semantics": "expected_edit_fraction",
            "training_overlap_status": "declared_no_overlap",
            "evidence_reference": "https://example.test/model",
        }
    )
    expected = (0.5, 1 / 3, 0.5)
    for record, score in zip(predictions["records"], expected, strict=True):
        record["predicted_score"] = score

    report = evaluate_cynomolgus_base_editing_transfer(
        targets,
        source_data,
        predictions,
        source=source,
    )

    assert report.species_profile == "cynomolgus_macaque"
    assert report.evaluation_status == (
        "retrospective_external_transfer_benchmark"
    )
    assert report.predictive_adapter_available is False
    assert report.source_verification == "pinned_workbooks_verified"
    assert report.source_assembly_accession == "synthetic-build"
    assert report.target_assembly_accession == "GCF_037993035.2"
    assert report.liftover_status == "not_performed"
    assert report.publisher_target_sequence_record_verified is True
    assert report.target_sequence_verified_on_source_assembly is False
    assert report.target_sequence_verified_on_target is False
    assert report.geneimpact_version == "0.16.1"
    assert report.evaluator_code_revision.startswith("sha256:")
    assert report.evaluator_code_revision_verified is True
    assert report.evaluator_code_revision_status == "module_source_sha256"
    assert report.record_count == 3
    assert report.context_count == 2
    assert report.comparison_stratum_count == 2
    assert report.embryo_base_observation_count == 6
    assert report.clone_denominator_count == 80
    assert report.metrics.within_context_candidate_pair_count == 1
    assert report.metrics.within_context_eligible_pair_count == 1
    assert report.metrics.within_context_observation_tie_pair_count == 0
    assert report.metrics.within_context_prediction_tie_pair_count == 0
    assert report.metrics.within_context_pair_count == 1
    assert report.metrics.within_context_prediction_coverage == pytest.approx(
        1.0
    )
    assert report.metrics.within_context_pairwise_accuracy == pytest.approx(1.0)
    assert report.metrics.mean_absolute_error == pytest.approx(0.0)
    assert report.metrics.root_mean_squared_error == pytest.approx(0.0)
    assert report.independence_verified is False
    assert "records" not in asdict(report)
    assert "observed_fraction" not in repr(report)
    assert not any(
        sequence in repr(report) for sequence in TARGET_SEQUENCES.values()
    )


def test_does_not_compare_scores_across_editors_in_one_injection_context(
    tmp_path,
):
    targets, source_data, source = _write_sources(tmp_path)
    mixed_context_source = replace(
        source,
        expected_context_count=1,
        blocks=tuple(
            replace(block, context_id="mixed_editor_injection")
            for block in source.blocks
        ),
    )
    predictions = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=mixed_context_source,
    )
    predictions["prediction"].update(
        {
            "name": "mixed-editor-test",
            "version": "test-v1",
            "submitted_code_revision": "test-commit",
            "score_semantics": "ranking_score",
            "training_overlap_status": "declared_no_overlap",
            "evidence_reference": "https://example.test/model",
        }
    )
    for record, score in zip(
        predictions["records"],
        (0.5, 0.3, 0.9),
        strict=True,
    ):
        record["predicted_score"] = score

    report = evaluate_cynomolgus_base_editing_transfer(
        targets,
        source_data,
        predictions,
        source=mixed_context_source,
    )

    assert report.context_count == 1
    assert report.comparison_stratum_count == 2
    assert report.metrics.within_context_pair_count == 1


def test_prediction_ties_are_visible_and_score_as_half_credit(tmp_path):
    targets, source_data, source = _write_sources(tmp_path)
    predictions = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=source,
    )
    predictions["prediction"].update(
        {
            "name": "constant-ranking-model",
            "version": "test-v1",
            "submitted_code_revision": "test-commit",
            "score_semantics": "ranking_score",
            "training_overlap_status": "declared_no_overlap",
            "evidence_reference": "https://example.test/model",
        }
    )
    for record in predictions["records"]:
        record["predicted_score"] = 0.5

    report = evaluate_cynomolgus_base_editing_transfer(
        targets,
        source_data,
        predictions,
        source=source,
    )

    assert report.metrics.within_context_candidate_pair_count == 1
    assert report.metrics.within_context_observation_tie_pair_count == 0
    assert report.metrics.within_context_prediction_tie_pair_count == 1
    assert report.metrics.within_context_pair_count == 0
    assert report.metrics.within_context_prediction_coverage == pytest.approx(
        0.0
    )
    assert report.metrics.within_context_pairwise_accuracy == pytest.approx(0.5)
    assert report.metrics.within_context_weighted_concordant_score == pytest.approx(
        0.5
    )


def test_unknown_training_overlap_is_descriptive_not_external_validation(
    tmp_path,
):
    targets, source_data, source = _write_sources(tmp_path)
    predictions = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=source,
    )
    predictions["prediction"].update(
        {
            "name": "unknown-overlap-model",
            "version": "test-v1",
            "submitted_code_revision": "test-commit",
            "score_semantics": "ranking_score",
            "training_overlap_status": "unknown",
            "evidence_reference": "https://example.test/model",
        }
    )
    for index, record in enumerate(predictions["records"]):
        record["predicted_score"] = float(index)

    report = evaluate_cynomolgus_base_editing_transfer(
        targets,
        source_data,
        predictions,
        source=source,
    )

    assert report.evaluation_status == (
        "descriptive_evaluation_with_unverified_overlap"
    )
    assert report.use == "descriptive_only_unverified_overlap"


def test_ranking_scores_do_not_trigger_fraction_error_overflow(tmp_path):
    targets, source_data, source = _write_sources(tmp_path)
    predictions = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=source,
    )
    predictions["prediction"].update(
        {
            "name": "large-ranking-model",
            "version": "test-v1",
            "submitted_code_revision": "test-commit",
            "score_semantics": "ranking_score",
            "training_overlap_status": "declared_no_overlap",
            "evidence_reference": "https://example.test/model",
        }
    )
    for index, record in enumerate(predictions["records"]):
        record["predicted_score"] = 1e308 if index % 2 else -1e308

    report = evaluate_cynomolgus_base_editing_transfer(
        targets,
        source_data,
        predictions,
        source=source,
    )

    assert report.metrics.mean_absolute_error is None
    assert report.metrics.root_mean_squared_error is None


def test_rejects_tampered_publisher_workbook(tmp_path):
    targets, source_data, source = _write_sources(tmp_path)
    targets.write_bytes(targets.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="target-sites workbook SHA-256"):
        prepare_cynomolgus_base_editing_transfer_template(
            targets,
            source_data,
            source=source,
        )


def test_rejects_known_training_overlap(tmp_path):
    targets, source_data, source = _write_sources(tmp_path)
    predictions = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=source,
    )
    predictions["prediction"].update(
        {
            "name": "overlapping-model",
            "version": "test-v1",
            "submitted_code_revision": "test-commit",
            "training_overlap_status": "overlap_detected",
            "evidence_reference": "https://example.test/model",
        }
    )
    for record in predictions["records"]:
        record["predicted_score"] = 0.5

    with pytest.raises(ValueError, match="known training overlap"):
        evaluate_cynomolgus_base_editing_transfer(
            targets,
            source_data,
            predictions,
            source=source,
        )


def test_prepare_cli_writes_researcher_submission_template(
    tmp_path,
    monkeypatch,
    capsys,
):
    targets, source_data, source = _write_sources(tmp_path)
    output = tmp_path / "predictions.json"
    monkeypatch.setattr(
        cli,
        "ZHANG_2020_CYNOMOLGUS_BASE_EDITING_SOURCE",
        source,
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "geneimpact",
            "prepare-cynomolgus-base-editing-transfer",
            "--target-sites",
            str(targets),
            "--source-data",
            str(source_data),
            "--output",
            str(output),
        ],
    )

    main()

    template = json.loads(output.read_text(encoding="utf-8"))
    assert len(template["records"]) == 3
    assert template["records"][0]["predicted_score"] is None
    assert "written" in capsys.readouterr().out


def test_validate_cli_writes_bounded_transfer_report(
    tmp_path,
    monkeypatch,
    capsys,
):
    targets, source_data, source = _write_sources(tmp_path)
    predictions = prepare_cynomolgus_base_editing_transfer_template(
        targets,
        source_data,
        source=source,
    )
    predictions["prediction"].update(
        {
            "name": "synthetic-ranking-model",
            "version": "test-v1",
            "submitted_code_revision": "test-commit",
            "training_overlap_status": "unknown",
            "evidence_reference": "https://example.test/model",
        }
    )
    for index, record in enumerate(predictions["records"]):
        record["predicted_score"] = float(index)
    predictions_path = tmp_path / "predictions.json"
    predictions_path.write_text(json.dumps(predictions), encoding="utf-8")
    output = tmp_path / "report.json"
    monkeypatch.setattr(
        cli,
        "ZHANG_2020_CYNOMOLGUS_BASE_EDITING_SOURCE",
        source,
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "geneimpact",
            "validate-cynomolgus-base-editing-transfer",
            "--target-sites",
            str(targets),
            "--source-data",
            str(source_data),
            "--predictions",
            str(predictions_path),
            "--output",
            str(output),
        ],
    )

    main()

    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["record_count"] == 3
    assert report["predictive_adapter_available"] is False
    assert report["evaluation_status"] == (
        "descriptive_evaluation_with_unverified_overlap"
    )
    assert report["metrics"]["mean_absolute_error"] is None
    assert "written" in capsys.readouterr().out
