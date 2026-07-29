from urllib.parse import parse_qs, urlparse

import pytest

from geneimpact.impc import ImpcClient


def response():
    return {
        "response": {
            "numFound": 1,
            "docs": [{
                "marker_symbol": "Prkdc",
                "mp_term_id": "MP:0000221",
                "mp_term_name": "decreased leukocyte cell number",
                "top_level_mp_term_name": ["immune system phenotype"],
                "effect_size": -3.2,
                "p_value": 1e-12,
                "significant": True,
                "procedure_name": "Hematology",
                "parameter_name": "White blood cell count",
                "sex": "not_considered",
                "zygosity": "homozygote",
            }],
        }
    }


def test_impc_client_builds_bounded_significant_query():
    requested = []
    evidence = ImpcClient(reader=lambda url: requested.append(url) or response()).significant_gene_phenotypes(
        "Prkdc"
    )
    query = parse_qs(urlparse(requested[0]).query)

    assert query["q"] == ["marker_symbol:Prkdc AND significant:true"]
    assert query["rows"] == ["1000"]
    assert query["start"] == ["0"]
    assert evidence.num_found == 1
    assert evidence.results[0].mp_term_id == "MP:0000221"
    assert evidence.results[0].significant


def test_impc_client_rejects_unbounded_request():
    with pytest.raises(ValueError, match="between 1 and 1000"):
        ImpcClient(reader=lambda _: response()).significant_gene_phenotypes(
            "Prkdc", rows=1001
        )


def test_impc_client_can_fetch_tested_non_significant_results():
    requested = []
    ImpcClient(reader=lambda url: requested.append(url) or response()).gene_phenotypes(
        "Prkdc", significant=None
    )
    query = parse_qs(urlparse(requested[0]).query)

    assert query["q"] == ["marker_symbol:Prkdc"]


def test_impc_client_paginates_within_run_limit():
    requested = []

    def reader(url):
        requested.append(url)
        start = int(parse_qs(urlparse(url).query)["start"][0])
        docs = [
            {"marker_symbol": "Gene", "significant": False}
            for _ in range(2 if start == 0 else 1)
        ]
        return {"response": {"numFound": 3, "docs": docs}}

    evidence = ImpcClient(reader=reader).gene_phenotypes(
        "Gene", rows=2, max_documents=3
    )

    assert evidence.pages == 2
    assert len(evidence.results) == 3
    assert len(requested) == 2
