import json

import pytest

from geneimpact.baseline import evaluate_benchmark, evaluate_ranking, fit_phenotype_prior


def row(gene, phenotype):
    return {"gene_symbol": gene, "phenotype_id": phenotype}


def test_frequency_prior_and_grouped_recall():
    model = fit_phenotype_prior([
        row("A", "MP:1"),
        row("B", "MP:1"),
        row("C", "MP:2"),
    ])
    metrics = evaluate_ranking(
        model,
        [row("D", "MP:1"), row("D", "MP:3"), row("E", "MP:2")],
        k=1,
    )

    assert model.ranked_phenotype_ids[0] == "MP:1"
    assert metrics.genes == 2
    assert metrics.macro_recall_at_k == pytest.approx(0.25)
    assert metrics.gene_hit_rate_at_k == pytest.approx(0.5)


def test_benchmark_report_is_bound_to_manifest(tmp_path):
    (tmp_path / "manifest.json").write_text('{"version": 1}\n', encoding="utf-8")
    splits = {
        "train": [row("A", "MP:1"), row("B", "MP:2")],
        "validation": [row("C", "MP:1")],
        "test": [row("D", "MP:2")],
    }
    for split, records in splits.items():
        (tmp_path / f"{split}.jsonl").write_text(
            "".join(json.dumps(record) + "\n" for record in records),
            encoding="utf-8",
        )

    report = evaluate_benchmark(tmp_path, k=1)

    assert report.benchmark_manifest_sha256
    assert report.validation.genes == 1
    assert (tmp_path / "baseline-report.json").exists()
