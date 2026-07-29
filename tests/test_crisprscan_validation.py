import json
from pathlib import Path

import pytest

from geneimpact.crisprscan_validation import evaluate_crisprscan_transfer


DATASET = (
    Path(__file__).parents[1]
    / "data"
    / "benchmarks"
    / "crisprscan-nhgri1-2022.json"
)


def _dataset():
    return json.loads(DATASET.read_text(encoding="utf-8"))


def test_reproduces_independent_transfer_metrics():
    report = evaluate_crisprscan_transfer(_dataset())

    assert report.record_count == 50
    assert report.gene_count == 14
    assert report.pearson_correlation == pytest.approx(0.2711476697010056)
    assert report.spearman_correlation == pytest.approx(0.2732474227858961)
    assert report.within_gene_pair_count == 64
    assert report.within_gene_concordant_pair_count == 36
    assert report.within_gene_pairwise_accuracy == pytest.approx(0.5625)
    assert report.domain_fit == "outside_declared_t7_expression_domain"
    assert report.training_exact_guide_overlap_count == 0
    assert report.training_reverse_complement_overlap_count == 0


def test_variant_aware_reported_scores_have_separate_metrics():
    report = evaluate_crisprscan_transfer(
        _dataset(),
        prediction_field="reported_crisprscan_nhgri_score",
    )

    assert report.pearson_correlation == pytest.approx(0.35179043379544656)
    assert report.spearman_correlation == pytest.approx(0.3105878356336742)
    assert report.within_gene_pair_count == 67
    assert report.within_gene_concordant_pair_count == 38
    assert report.within_gene_pairwise_accuracy == pytest.approx(
        0.5671641791044776
    )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("independent_lab", False, "independent_lab"),
        ("training_data_overlap", True, "training_data_overlap"),
        ("species_profile", "mouse", "zebrafish"),
        ("benchmark_scope", "prospective", "retrospective_external_transfer"),
        ("source_workbook_sha256", "bad", "SHA-256"),
        ("training_source_workbook_sha256", "bad", "SHA-256"),
        ("training_exact_guide_overlap_count", 1, "overlap"),
    ],
)
def test_rejects_ineligible_metadata(field, value, message):
    dataset = _dataset()
    dataset["metadata"][field] = value

    with pytest.raises(ValueError, match=message):
        evaluate_crisprscan_transfer(dataset)


def test_rejects_duplicate_guides():
    dataset = _dataset()
    dataset["records"][1]["guide_sha256"] = dataset["records"][0]["guide_sha256"]

    with pytest.raises(ValueError, match="duplicate guide"):
        evaluate_crisprscan_transfer(dataset)


def test_rejects_unsupported_prediction_field():
    with pytest.raises(ValueError, match="prediction_field"):
        evaluate_crisprscan_transfer(_dataset(), "observed_indel_fraction")
