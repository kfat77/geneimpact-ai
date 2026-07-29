import json

import pytest

from geneimpact.behive import BEHIVE_EFFICIENCY_COMMIT
from geneimpact.behive_validation import evaluate_behive_validation


def _document():
    return {
        "dataset": {
            "dataset_name": "independent-mes-be4",
            "source_reference": "doi:example",
            "species": "mouse",
            "cell_type": "mES",
            "base_editor": "BE4",
            "model_commit": BEHIVE_EFFICIENCY_COMMIT,
            "independent_lab": True,
            "training_sequence_overlap_count": 0,
            "training_sequence_overlap_method": "Exact SHA-256 comparison against training targets.",
        },
        "records": [
            {
                "sequence_sha256": f"{index:064x}",
                "base_editor": "BE4",
                "model_commit": BEHIVE_EFFICIENCY_COMMIT,
                "observed_fraction": observed,
                "predicted_logit_score": logit,
                "calibrated_fraction": calibrated,
            }
            for index, (observed, logit, calibrated) in enumerate(
                [(0.1, -1.0, 0.12), (0.4, 0.0, 0.35), (0.9, 1.0, 0.88)],
                start=1,
            )
        ],
    }


def _write(tmp_path, document):
    path = tmp_path / "validation.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def test_evaluates_independent_single_editor_dataset(tmp_path):
    report = evaluate_behive_validation(_write(tmp_path, _document()))

    assert report.metrics.record_count == 3
    assert report.metrics.pearson_logit_vs_observed > 0.98
    assert report.metrics.calibrated_mae == pytest.approx(0.03)
    assert report.training_sequence_overlap_count == 0
    assert len(report.input_sha256) == 64


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("independent_lab", False, "independent_lab"),
        ("training_sequence_overlap_count", 1, "must be zero"),
        ("species", "human", "mouse dataset"),
        ("cell_type", "HEK293T", "cell_type mES"),
    ],
)
def test_rejects_out_of_domain_or_non_independent_dataset(tmp_path, field, value, message):
    document = _document()
    document["dataset"][field] = value

    with pytest.raises(ValueError, match=message):
        evaluate_behive_validation(_write(tmp_path, document))


def test_rejects_duplicate_sequence_hash(tmp_path):
    document = _document()
    document["records"][1]["sequence_sha256"] = document["records"][0]["sequence_sha256"]

    with pytest.raises(ValueError, match="duplicates"):
        evaluate_behive_validation(_write(tmp_path, document))


def test_allows_uncalibrated_ranking_evaluation(tmp_path):
    document = _document()
    for record in document["records"]:
        record["calibrated_fraction"] = None

    report = evaluate_behive_validation(_write(tmp_path, document))

    assert report.metrics.pearson_logit_vs_observed > 0.98
    assert report.metrics.calibrated_mae is None
    assert report.metrics.calibrated_rmse is None
