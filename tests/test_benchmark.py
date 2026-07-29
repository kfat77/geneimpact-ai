import json

from geneimpact.benchmark import assign_gene_split, build_mgi_benchmark


def record(allele_id, gene, phenotype_ids, origin=None):
    return {
        "allele_id": allele_id,
        "marker_symbol": gene,
        "ensembl_gene_id": f"ENS_{gene}",
        "allele_attributes": ["Null/knockout"],
        "high_level_mp_ids": phenotype_ids,
        "origin_program": origin,
    }


def test_gene_grouping_is_deterministic():
    assert assign_gene_split("Prkdc") == assign_gene_split("Prkdc")
    assert assign_gene_split("Prkdc") == assign_gene_split("prkdc")


def test_benchmark_excludes_impc_origin_and_writes_manifest(tmp_path):
    source = tmp_path / "normalized.jsonl"
    rows = [
        record("MGI:1", "GeneA", ["MP:1", "MP:2"]),
        record("MGI:2", "GeneA", ["MP:3"]),
        record("MGI:3", "GeneB", ["MP:4"], "IMPC"),
        record("MGI:4", "GeneC", []),
    ]
    source.write_text(
        "".join(json.dumps(row) + "\n" for row in rows),
        encoding="utf-8",
    )

    manifest = build_mgi_benchmark(source, tmp_path / "benchmark")
    outputs = []
    for split in ("train", "validation", "test"):
        path = tmp_path / "benchmark" / f"{split}.jsonl"
        outputs.extend(json.loads(line) for line in path.read_text().splitlines())

    assert manifest.excluded_impc_alleles == 1
    assert manifest.output_associations == 3
    assert {row["gene_symbol"] for row in outputs} == {"GeneA"}
    assert len({row["split"] for row in outputs}) == 1
    assert (tmp_path / "benchmark" / "manifest.json").exists()
