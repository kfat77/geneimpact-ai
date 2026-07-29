"""Evidence-status registry for species-specific prediction capabilities."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .species import PROFILES


class CapabilityStatus(str, Enum):
    AVAILABLE_DECLARED_DOMAIN = "available_declared_domain"
    VALIDATION_CANDIDATE = "validation_candidate"
    REFERENCE_SEARCH_CANDIDATE = "reference_search_candidate"
    OUT_OF_DOMAIN_ONLY = "out_of_domain_only"


@dataclass(frozen=True)
class PredictorCapability:
    species_profile: str
    predictor: str
    task: str
    edit_classes: tuple[str, ...]
    status: CapabilityStatus
    biological_domain: str
    evidence_reference: str
    note: str


_CRISPRITZ_REFERENCE = "https://pubmed.ncbi.nlm.nih.gov/31764961/"
_CRISPRSCAN_REFERENCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4589495/"
_INDELPHI_REFERENCE = "https://www.nature.com/articles/s41586-018-0686-x"


def capabilities_for_species(species_profile: str) -> tuple[PredictorCapability, ...]:
    """Return explicit available/candidate capabilities for one profile."""
    if species_profile not in PROFILES:
        raise ValueError(f"unknown species profile {species_profile!r}.")

    capabilities = [
        PredictorCapability(
            species_profile=species_profile,
            predictor="CRISPRitz",
            task="reference_genome_off_target_enumeration",
            edit_classes=("knockout", "base_editing", "prime_editing"),
            status=CapabilityStatus.REFERENCE_SEARCH_CANDIDATE,
            biological_domain=f"{species_profile} reference assembly",
            evidence_reference=_CRISPRITZ_REFERENCE,
            note=(
                "Genome-indexed mismatch/bulge search is species-configurable, "
                "but a version-locked GeneImpact adapter and species-specific empirical calibration are pending."
            ),
        )
    ]

    if species_profile == "mouse":
        capabilities.extend(
            (
                PredictorCapability(
                    species_profile="mouse",
                    predictor="BE-Hive efficiency",
                    task="base_editing_efficiency",
                    edit_classes=("base_editing",),
                    status=CapabilityStatus.AVAILABLE_DECLARED_DOMAIN,
                    biological_domain="mouse embryonic stem cells (mES)",
                    evidence_reference=(
                        "https://github.com/maxwshen/be_predict_efficiency/"
                        "tree/fbd495910d6c95b24081649015d6257a8badc9d7"
                    ),
                    note="Import adapter available; independent mES validation remains pending.",
                ),
                PredictorCapability(
                    species_profile="mouse",
                    predictor="BE-Hive bystander",
                    task="base_editing_bystander_outcomes",
                    edit_classes=("base_editing",),
                    status=CapabilityStatus.AVAILABLE_DECLARED_DOMAIN,
                    biological_domain="mouse embryonic stem cells (mES)",
                    evidence_reference=(
                        "https://github.com/maxwshen/be_predict_bystander/"
                        "tree/31aadd04a25c604857c7592b226ee987e9e20b31"
                    ),
                    note="Import adapter available; outcome frequencies are not phenotype probabilities.",
                ),
                PredictorCapability(
                    species_profile="mouse",
                    predictor="inDelphi",
                    task="repair_outcome",
                    edit_classes=("knockout",),
                    status=CapabilityStatus.VALIDATION_CANDIDATE,
                    biological_domain="mouse embryonic stem-cell repair data",
                    evidence_reference=_INDELPHI_REFERENCE,
                    note="Candidate for the next version-locked repair-outcome adapter.",
                ),
            )
        )
    elif species_profile == "zebrafish":
        capabilities.append(
            PredictorCapability(
                species_profile="zebrafish",
                predictor="CRISPRscan",
                task="guide_activity",
                edit_classes=("knockout",),
                status=CapabilityStatus.VALIDATION_CANDIDATE,
                biological_domain="in-vivo zebrafish embryo mutagenesis; predominantly TU background",
                evidence_reference=_CRISPRSCAN_REFERENCE,
                note=(
                    "Species-relevant candidate; adapter must account for reference-build "
                    "and strain differences before reporting applicable scores."
                ),
            )
        )
    elif species_profile in {"rat", "fruit_fly"}:
        capabilities.append(
            PredictorCapability(
                species_profile=species_profile,
                predictor="CRISPRscan",
                task="guide_activity",
                edit_classes=("knockout",),
                status=CapabilityStatus.OUT_OF_DOMAIN_ONLY,
                biological_domain="model trained on zebrafish in-vivo activity",
                evidence_reference=_CRISPRSCAN_REFERENCE,
                note=(
                    "The service exposes this genome, but the activity model is not treated "
                    "as species-validated for this profile."
                ),
            )
        )

    return tuple(capabilities)


def capability_matrix() -> dict[str, tuple[PredictorCapability, ...]]:
    """Return the full registered matrix in deterministic profile order."""
    return {
        profile_key: capabilities_for_species(profile_key)
        for profile_key in PROFILES
    }
