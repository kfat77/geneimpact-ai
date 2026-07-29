from geneimpact.datasources import EnsemblMetadataClient, check_ensembl_profile
from geneimpact.species import MOUSE_PROFILE


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
