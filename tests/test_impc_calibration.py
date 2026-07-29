import json

import pytest

from geneimpact.impc_calibration import evaluate_impc_calibration


def write(path, rows):
    path.write_text(
        "".join(json.dumps(row) + "\n" for row in rows),
        encoding="utf-8",
    )


def test_impc_calibration_reports_probability_metrics(tmp_path):
    calibration = tmp_path / "calibration.jsonl"
    test = tmp_path / "test.jsonl"
    output = tmp_path / "report.json"
    write(calibration, [
        {"gene_symbol": "A", "significant": True},
        {"gene_symbol": "A", "significant": False},
        {"gene_symbol": "B", "significant": False},
    ])
    write(test, [
        {"gene_symbol": "C", "significant": True},
        {"gene_symbol": "C", "significant": False},
    ])

    report = evaluate_impc_calibration(calibration, test, output_path=output)

    assert report.model.probability == pytest.approx(0.375)
    assert report.test.genes == 1
    assert report.test.brier_score > 0
    assert report.gene_overlap_checked
    assert output.exists()


def test_impc_calibration_rejects_gene_overlap(tmp_path):
    calibration = tmp_path / "calibration.jsonl"
    test = tmp_path / "test.jsonl"
    write(calibration, [{"gene_symbol": "A", "significant": False}])
    write(test, [{"gene_symbol": "A", "significant": True}])

    with pytest.raises(ValueError, match="gene leakage"):
        evaluate_impc_calibration(calibration, test)
