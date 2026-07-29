"""Declared species profiles and applicability checks."""

from __future__ import annotations

from dataclasses import dataclass

from .provenance import StudyContext


@dataclass(frozen=True)
class SpeciesProfile:
    """Versioned applicability boundary for a supported research species."""

    key: str
    scientific_name: str
    taxon_id: str
    genome_build: str
    assembly_accession: str
    reference_strain: str
    accepted_strain_names: tuple[str, ...]


@dataclass(frozen=True)
class SpeciesValidation:
    """Structured validation result retained in an assessment report."""

    supported: bool
    profile_key: str | None
    errors: tuple[str, ...]
    warnings: tuple[str, ...]


MOUSE_PROFILE = SpeciesProfile(
    key="mouse",
    scientific_name="mus_musculus",
    taxon_id="10090",
    genome_build="GRCm39",
    assembly_accession="GCA_000001635.9",
    reference_strain="C57BL/6J",
    accepted_strain_names=("C57BL/6", "C57BL/6J"),
)

PROFILES = {MOUSE_PROFILE.key: MOUSE_PROFILE}


def validate_study_context(context: StudyContext) -> SpeciesValidation:
    """Validate a study against a declared profile without silent fallback."""
    key = context.species.strip().casefold()
    profile = PROFILES.get(key)
    if profile is None:
        return SpeciesValidation(
            supported=False,
            profile_key=None,
            errors=(f"species {context.species!r} is not yet registered.",),
            warnings=(),
        )

    errors: list[str] = []
    warnings: list[str] = []
    if context.genome_build.casefold() != profile.genome_build.casefold():
        errors.append(
            f"genome build {context.genome_build!r} does not match registered "
            f"{profile.genome_build!r} ({profile.assembly_accession})."
        )
    accepted = {name.casefold() for name in profile.accepted_strain_names}
    if context.strain_or_breed.casefold() not in accepted:
        warnings.append(
            f"strain {context.strain_or_breed!r} differs from the registered "
            f"reference strain {profile.reference_strain!r}; strain-specific validation is required."
        )
    return SpeciesValidation(
        supported=not errors,
        profile_key=profile.key,
        errors=tuple(errors),
        warnings=tuple(warnings),
    )
