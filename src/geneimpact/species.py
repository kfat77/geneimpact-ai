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
    accepted_build_names: tuple[str, ...] = ()
    aliases: tuple[str, ...] = ()
    ensembl_genome_build: str | None = None
    ensembl_assembly_accession: str | None = None


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
    assembly_accession="GCF_000001635.27",
    reference_strain="C57BL/6J",
    accepted_strain_names=("C57BL/6", "C57BL/6J"),
    aliases=("mus musculus", "mus_musculus", "小鼠"),
    ensembl_genome_build="GRCm39",
    ensembl_assembly_accession="GCA_000001635.9",
)

RAT_PROFILE = SpeciesProfile(
    key="rat",
    scientific_name="rattus_norvegicus",
    taxon_id="10116",
    genome_build="GRCr8",
    assembly_accession="GCF_036323735.1",
    reference_strain="BN/NHsdMcwi",
    accepted_strain_names=("BN/NHsdMcwi", "BN"),
    aliases=("rattus norvegicus", "rattus_norvegicus", "大鼠"),
    ensembl_genome_build="GRCr8",
    ensembl_assembly_accession="GCA_036323735.1",
)

ZEBRAFISH_PROFILE = SpeciesProfile(
    key="zebrafish",
    scientific_name="danio_rerio",
    taxon_id="7955",
    genome_build="GRCz12tu",
    assembly_accession="GCF_049306965.2",
    reference_strain="Tuebingen",
    accepted_strain_names=("Tuebingen", "Tübingen", "TU"),
    aliases=("danio rerio", "danio_rerio", "斑马鱼"),
    ensembl_genome_build="GRCz11",
    ensembl_assembly_accession="GCA_000002035.4",
)

FRUIT_FLY_PROFILE = SpeciesProfile(
    key="fruit_fly",
    scientific_name="drosophila_melanogaster",
    taxon_id="7227",
    genome_build="Release 6 plus ISO1 MT",
    assembly_accession="GCF_000001215.4",
    reference_strain="ISO-1",
    accepted_strain_names=("ISO-1", "iso-1"),
    accepted_build_names=("BDGP6.54", "dm6"),
    aliases=(
        "fruit fly",
        "drosophila",
        "drosophila melanogaster",
        "drosophila_melanogaster",
        "果蝇",
    ),
    ensembl_genome_build="BDGP6.54",
    ensembl_assembly_accession="GCA_000001215.4",
)

RHESUS_MACAQUE_PROFILE = SpeciesProfile(
    key="rhesus_macaque",
    scientific_name="macaca_mulatta",
    taxon_id="9544",
    genome_build="T2T-MMU8v2.0",
    assembly_accession="GCF_049350105.2",
    reference_strain="MMU2019108-1",
    accepted_strain_names=("MMU2019108-1",),
    aliases=("rhesus macaque", "rhesus monkey", "macaca mulatta", "macaca_mulatta", "恒河猴"),
    ensembl_genome_build="Mmul_10",
    ensembl_assembly_accession="GCA_003339765.3",
)

CYNOMOLGUS_MACAQUE_PROFILE = SpeciesProfile(
    key="cynomolgus_macaque",
    scientific_name="macaca_fascicularis",
    taxon_id="9541",
    genome_build="T2T-MFA8v1.1",
    assembly_accession="GCF_037993035.2",
    reference_strain="582-1",
    accepted_strain_names=("582-1",),
    aliases=(
        "cynomolgus macaque",
        "crab-eating macaque",
        "macaca fascicularis",
        "macaca_fascicularis",
        "食蟹猴",
    ),
    ensembl_genome_build="Macaca_fascicularis_6.0",
    ensembl_assembly_accession="GCA_011100615.1",
)

PROFILES = {
    profile.key: profile
    for profile in (
        MOUSE_PROFILE,
        RAT_PROFILE,
        ZEBRAFISH_PROFILE,
        FRUIT_FLY_PROFILE,
        RHESUS_MACAQUE_PROFILE,
        CYNOMOLGUS_MACAQUE_PROFILE,
    )
}
PROFILE_ALIASES = {
    alias.casefold(): profile
    for profile in PROFILES.values()
    for alias in (profile.key, *profile.aliases)
}
AMBIGUOUS_ALIASES = {
    "monkey": ("rhesus_macaque", "cynomolgus_macaque"),
    "macaque": ("rhesus_macaque", "cynomolgus_macaque"),
    "猴": ("rhesus_macaque", "cynomolgus_macaque"),
    "猕猴": ("rhesus_macaque", "cynomolgus_macaque"),
}


def validate_study_context(context: StudyContext) -> SpeciesValidation:
    """Validate a study against a declared profile without silent fallback."""
    key = context.species.strip().casefold()
    if key in AMBIGUOUS_ALIASES:
        choices = ", ".join(AMBIGUOUS_ALIASES[key])
        return SpeciesValidation(
            supported=False,
            profile_key=None,
            errors=(
                f"species {context.species!r} is ambiguous; choose one registered species: {choices}.",
            ),
            warnings=(),
        )
    profile = PROFILE_ALIASES.get(key)
    if profile is None:
        return SpeciesValidation(
            supported=False,
            profile_key=None,
            errors=(f"species {context.species!r} is not yet registered.",),
            warnings=(),
        )

    errors: list[str] = []
    warnings: list[str] = []
    accepted_builds = {
        name.casefold()
        for name in (profile.genome_build, *profile.accepted_build_names)
    }
    if context.genome_build.casefold() not in accepted_builds:
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
