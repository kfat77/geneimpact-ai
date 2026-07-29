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


def test_indelphi_external_import_is_available_only_for_mouse_knockout():
    mouse = capabilities_for_species("mouse")
    indelphi = next(item for item in mouse if item.predictor == "inDelphi")

    assert indelphi.status is CapabilityStatus.AVAILABLE_DECLARED_DOMAIN
    assert indelphi.edit_classes == ("knockout",)
    assert all(
        not any(
            item.predictor == "inDelphi"
            for item in capabilities_for_species(profile)
        )
        for profile in PROFILES
        if profile != "mouse"
    )


def test_crisprscan_is_available_only_in_zebrafish_domain():
    zebrafish = capabilities_for_species("zebrafish")
    fruit_fly = capabilities_for_species("fruit_fly")

    assert next(item for item in zebrafish if item.predictor == "CRISPRscan").status is (
        CapabilityStatus.AVAILABLE_DECLARED_DOMAIN
    )
    assert next(item for item in fruit_fly if item.predictor == "CRISPRscan").status is (
        CapabilityStatus.OUT_OF_DOMAIN_ONLY
    )

    cas12a = next(
        item
        for item in fruit_fly
        if item.task == "in_vivo_cas12a_array_loh_evidence"
    )
    assert cas12a.status is CapabilityStatus.USABLE_BOUNDED_BENCHMARK
    assert cas12a.predictor == "Port 2026 Cas12a array LOH evidence"
    assert "not a predictor" in cas12a.note


def test_rat_transfer_benchmark_is_candidate_not_available_predictor():
    rat = capabilities_for_species("rat")
    benchmark = next(
        item
        for item in rat
        if item.task == "guide_activity_transfer_validation"
    )

    assert benchmark.status is CapabilityStatus.VALIDATION_CANDIDATE
    assert not any(
        item.status is CapabilityStatus.AVAILABLE_DECLARED_DOMAIN
        for item in rat
    )


def test_unknown_profile_is_rejected():
    with pytest.raises(ValueError, match="unknown species"):
        capabilities_for_species("monkey")
