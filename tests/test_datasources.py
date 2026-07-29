from geneimpact.datasources import (
    EnsemblMetadataClient,
    NcbiDatasetsClient,
    check_ensembl_profile,
    check_ncbi_profile,
)
from geneimpact.species import MOUSE_PROFILE, ZEBRAFISH_PROFILE


def payload(assembly: str = "GRCm39"):
    return {
        "species": [{
            "name": "mus_musculus",
            "taxon_id": "10090",
            "assembly": assembly,
            "accession": "GCA_000001635.9",
            "release": 116,
        }]
    }


def test_ensembl_profile_matches_versioned_mouse_metadata():
    check = check_ensembl_profile(
        MOUSE_PROFILE, EnsemblMetadataClient(reader=lambda _: payload())
    )

    assert check.matches
    assert check.checked_release == "116"


def test_ensembl_profile_detects_source_drift():
    check = check_ensembl_profile(
        MOUSE_PROFILE, EnsemblMetadataClient(reader=lambda _: payload("future-build"))
    )

    assert not check.matches
    assert "assembly changed" in check.errors[0]


def ncbi_payload(
    *,
    accession: str = "GCF_049306965.2",
    assembly: str = "GRCz12tu",
    status: str = "current",
):
    return {
        "reports": [
            {
                "accession": accession,
                "organism": {
                    "organism_name": "Danio rerio",
                    "tax_id": 7955,
                },
                "assembly_info": {
                    "assembly_name": assembly,
                    "assembly_status": status,
                    "refseq_category": "reference genome",
                    "release_date": "2025-04-04",
                },
            }
        ]
    }


def test_ncbi_profile_matches_current_zebrafish_reference():
    check = check_ncbi_profile(
        ZEBRAFISH_PROFILE,
        NcbiDatasetsClient(reader=lambda _: ncbi_payload()),
    )

    assert check.matches
    assert check.source == "NCBI Datasets"
    assert check.checked_release == "2025-04-04"


def test_ncbi_profile_detects_superseded_assembly():
    check = check_ncbi_profile(
        ZEBRAFISH_PROFILE,
        NcbiDatasetsClient(reader=lambda _: ncbi_payload(status="suppressed")),
    )

    assert not check.matches
    assert "status changed" in check.errors[0]
