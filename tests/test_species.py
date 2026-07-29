from geneimpact.provenance import StudyContext
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
