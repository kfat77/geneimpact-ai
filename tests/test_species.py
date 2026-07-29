from geneimpact.provenance import StudyContext
import pytest

from geneimpact.species import validate_study_context


def context(species: str = "mouse", strain: str = "C57BL/6", build: str = "GRCm39"):
    return StudyContext(species, strain, build, "knockout", "snapshot")


def test_mouse_reference_context_is_supported():
    result = validate_study_context(context())

    assert result.supported
    assert result.errors == ()
    assert result.warnings == ()


def test_wrong_build_is_rejected_and_non_reference_strain_is_flagged():
    result = validate_study_context(context(strain="BALB/c", build="GRCm38"))

    assert not result.supported
    assert "GRCm39" in result.errors[0]
    assert "strain-specific validation" in result.warnings[0]


@pytest.mark.parametrize(
    ("species", "strain", "build", "profile_key"),
    [
        ("rat", "BN/NHsdMcwi", "GRCr8", "rat"),
        ("斑马鱼", "Tuebingen", "GRCz12tu", "zebrafish"),
        ("fruit fly", "ISO-1", "dm6", "fruit_fly"),
        ("恒河猴", "MMU2019108-1", "T2T-MMU8v2.0", "rhesus_macaque"),
        ("食蟹猴", "582-1", "T2T-MFA8v1.1", "cynomolgus_macaque"),
    ],
)
def test_registered_multispecies_aliases_are_supported(
    species, strain, build, profile_key
):
    result = validate_study_context(context(species, strain, build))

    assert result.supported
    assert result.profile_key == profile_key
    assert result.errors == ()
    assert result.warnings == ()


def test_generic_monkey_is_rejected_as_ambiguous():
    result = validate_study_context(context("猴", "unknown", "unknown"))

    assert not result.supported
    assert "ambiguous" in result.errors[0]
    assert "rhesus_macaque" in result.errors[0]
    assert "cynomolgus_macaque" in result.errors[0]
