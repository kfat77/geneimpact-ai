import pytest

from geneimpact.capabilities import (
    CapabilityStatus,
    capabilities_for_species,
    capability_matrix,
)
from geneimpact.species import PROFILES


def test_every_registered_species_has_available_reference_search():
    matrix = capability_matrix()

    assert set(matrix) == set(PROFILES)
    for capabilities in matrix.values():
        assert any(
            item.predictor == "CRISPRitz"
            and item.status is CapabilityStatus.AVAILABLE_REFERENCE_SEARCH
            for item in capabilities
        )


def test_behive_is_available_only_for_mouse_mes():
    mouse = capabilities_for_species("mouse")
    non_mouse = [
        item
        for profile in PROFILES
        if profile != "mouse"
        for item in capabilities_for_species(profile)
    ]

    assert sum(
        item.predictor.startswith("BE-Hive")
        and item.status is CapabilityStatus.AVAILABLE_DECLARED_DOMAIN
        for item in mouse
    ) == 2
    assert not any(item.predictor.startswith("BE-Hive") for item in non_mouse)


def test_crisprscan_is_candidate_only_in_zebrafish_domain():
    zebrafish = capabilities_for_species("zebrafish")
    fruit_fly = capabilities_for_species("fruit_fly")

    assert next(item for item in zebrafish if item.predictor == "CRISPRscan").status is (
        CapabilityStatus.VALIDATION_CANDIDATE
    )
    assert next(item for item in fruit_fly if item.predictor == "CRISPRscan").status is (
        CapabilityStatus.OUT_OF_DOMAIN_ONLY
    )


def test_unknown_profile_is_rejected():
    with pytest.raises(ValueError, match="unknown species"):
        capabilities_for_species("monkey")
