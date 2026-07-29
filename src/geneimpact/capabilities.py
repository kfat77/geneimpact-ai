"""Evidence-status registry for species-specific prediction capabilities."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .housden import (
    HOUSDEN_EDIT_CLASS,
    HOUSDEN_METHOD_REFERENCE,
    HOUSDEN_PREDICTOR,
    HOUSDEN_SPECIES_PROFILE,
    HOUSDEN_TRAINING_DOMAIN,
)
from .fruit_fly_cas12a import FRUIT_FLY_CAS12A_REFERENCE
from .cynomolgus_base_editing import CYNOMOLGUS_BASE_EDITING_REFERENCE
from .species import PROFILES


class CapabilityStatus(str, Enum):
    AVAILABLE_DECLARED_DOMAIN = "available_declared_domain"
    AVAILABLE_REFERENCE_SEARCH = "available_reference_search"
    USABLE_BOUNDED_BENCHMARK = "usable_bounded_benchmark"
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
_RAT_TRANSFER_REFERENCE = "https://doi.org/10.1038/s41592-018-0011-5"


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
            status=CapabilityStatus.AVAILABLE_REFERENCE_SEARCH,
            biological_domain=f"{species_profile} reference assembly",
            evidence_reference=_CRISPRITZ_REFERENCE,
            note=(
                "Version-locked audit import is available for genome-indexed mismatch/bulge "
                "search; species-specific empirical calibration remains pending."
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
                    status=CapabilityStatus.AVAILABLE_DECLARED_DOMAIN,
                    biological_domain=(
                        "mESC integrated-target repair model; external transfer "
                        "evidence in C57BL/6JJmsSlc blastocysts"
                    ),
                    evidence_reference=_INDELPHI_REFERENCE,
                    note=(
                        "Version-locked external-result import is available. The "
                        "upstream license is restricted, the model does not predict "
                        "editing efficiency, and prospective mESC validation is required."
                    ),
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
                status=CapabilityStatus.AVAILABLE_DECLARED_DOMAIN,
                biological_domain="in-vivo zebrafish embryo mutagenesis; predominantly TU background",
                evidence_reference=_CRISPRSCAN_REFERENCE,
                note=(
                    "Version-locked scoring is available for canonical SpCas9, T7 in-vitro-"
                    "transcribed guides in zebrafish embryos; an independent RNP transfer "
                    "benchmark is available but does not close the in-domain validation gap."
                ),
            )
        )
    elif species_profile == HOUSDEN_SPECIES_PROFILE:
        capabilities.extend(
            (
                PredictorCapability(
                    species_profile=HOUSDEN_SPECIES_PROFILE,
                    predictor=HOUSDEN_PREDICTOR,
                    task="guide_activity_ranking",
                    edit_classes=(HOUSDEN_EDIT_CLASS,),
                    status=CapabilityStatus.AVAILABLE_DECLARED_DOMAIN,
                    biological_domain=HOUSDEN_TRAINING_DOMAIN,
                    evidence_reference=HOUSDEN_METHOD_REFERENCE,
                    note=(
                        "Official FlyRNAi external-result import is available. "
                        "The live service is unversioned, the score is not a "
                        "probability, and in-vivo validity is not established."
                    ),
                ),
                PredictorCapability(
                    species_profile="fruit_fly",
                    predictor="CRISPRscan",
                    task="guide_activity",
                    edit_classes=("knockout",),
                    status=CapabilityStatus.OUT_OF_DOMAIN_ONLY,
                    biological_domain="model trained on zebrafish in-vivo activity",
                    evidence_reference=_CRISPRSCAN_REFERENCE,
                    note=(
                        "The service exposes this genome, but the activity model is "
                        "not treated as species-validated for this profile."
                    ),
                ),
                PredictorCapability(
                    species_profile="fruit_fly",
                    predictor="Port 2026 Cas12a array LOH evidence",
                    task="in_vivo_cas12a_array_loh_evidence",
                    edit_classes=("knockout",),
                    status=CapabilityStatus.USABLE_BOUNDED_BENCHMARK,
                    biological_domain=(
                        "larval wing imaginal-disc LOH screens using HD12aCFD "
                        "three- or four-guide arrays"
                    ),
                    evidence_reference=FRUIT_FLY_CAS12A_REFERENCE,
                    note=(
                        "Checksum-pinned array-level evidence audit is "
                        "available. This is not a predictor: observations "
                        "cannot be assigned to component guides or interpreted "
                        "as calibrated probabilities."
                    ),
                ),
            )
        )
    elif species_profile == "rat":
        capabilities.extend(
            (
                PredictorCapability(
                    species_profile=species_profile,
                    predictor="external SpCas9 guide-activity model",
                    task="guide_activity_transfer_validation",
                    edit_classes=("knockout",),
                    status=CapabilityStatus.VALIDATION_CANDIDATE,
                    biological_domain=(
                        "rat G0 animals and embryos; 14 uniquely mapped rn5 guides"
                    ),
                    evidence_reference=_RAT_TRANSFER_REFERENCE,
                    note=(
                        "A pinned external-transfer evaluator is available, but no rat "
                        "predictor is promoted. The selected, high-activity set is too "
                        "small for training or calibration, and the source uses rn5."
                    ),
                ),
                PredictorCapability(
                    species_profile=species_profile,
                    predictor="CRISPRscan",
                    task="guide_activity",
                    edit_classes=("knockout",),
                    status=CapabilityStatus.OUT_OF_DOMAIN_ONLY,
                    biological_domain="model trained on zebrafish in-vivo activity",
                    evidence_reference=_CRISPRSCAN_REFERENCE,
                    note=(
                        "The service exposes this genome, but the activity model is not "
                        "treated as species-validated for this profile."
                    ),
                ),
            )
        )
    elif species_profile == "cynomolgus_macaque":
        capabilities.append(
            PredictorCapability(
                species_profile=species_profile,
                predictor=(
                    "Zhang 2020 cynomolgus embryo base-editing benchmark"
                ),
                task="base_editing_embryo_transfer_validation",
                edit_classes=("base_editing",),
                status=CapabilityStatus.USABLE_BOUNDED_BENCHMARK,
                biological_domain=(
                    "BE3, ABE7.10, and SaKKH-BE3 mRNA/T7-sgRNA editing in "
                    "cynomolgus zygotes across 11 published target sites"
                ),
                evidence_reference=CYNOMOLGUS_BASE_EDITING_REFERENCE,
                note=(
                    "A checksum-pinned external-score evaluator is available. "
                    "This is not a predictor or species-level calibration; "
                    "ranking is compared only within shared editor and "
                    "multiplex injection contexts."
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
