import json

import pytest

from geneimpact.impc import ImpcClient
from geneimpact.impc_validation import build_impc_validation


def response(gene):
    return {
        "response": {
            "numFound": 2,
            "docs": [
                {
                    "marker_symbol": gene,
                    "significant": True,
                    "mp_term_id": "MP:1",
                    "procedure_name": "Hematology",
                    "parameter_name": "Cell count",
                    "sex": "female",
                    "zygosity": "homozygote",
                },
                {
                    "marker_symbol": gene,
                    "significant": False,
                    "procedure_name": "Hematology",
                    "parameter_name": "Other count",
                    "sex": "female",
                    "zygosity": "homozygote",
                },
            ],
        }
    }


def test_impc_validation_preserves_tested_labels(tmp_path):
    client = ImpcClient(reader=lambda url: response("GeneA"))
    output = tmp_path / "impc.jsonl"

    manifest = build_impc_validation(["GeneA"], output, client=client)
    records = [json.loads(line) for line in output.read_text().splitlines()]

    assert manifest.significant_documents == 1
    assert manifest.non_significant_documents == 1
    assert {record["significant"] for record in records} == {True, False}
    assert (tmp_path / "impc.jsonl.manifest.json").exists()


def test_impc_validation_bounds_gene_count(tmp_path):
    with pytest.raises(ValueError, match="at most 2"):
        build_impc_validation(
            ["A", "B", "C"],
            tmp_path / "output.jsonl",
            client=ImpcClient(reader=lambda _: response("A")),
            max_genes=2,
        )
