"""Machine-readable evidence qualification for species-specific capabilities."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .species import PROFILES


class EvidenceUseStatus(str, Enum):
    """The strongest use justified by one evidence record."""

    USABLE_ADAPTER = "usable_adapter"
    TRANSFER_EVIDENCE_ONLY = "transfer_evidence_only"
    HAZARD_EVIDENCE_ONLY = "hazard_evidence_only"
    INSUFFICIENT_PUBLIC_DATA = "insufficient_public_data"


@dataclass(frozen=True)
class EvidenceQualification:
    record_id: str
    species_profile: str
    predictor_or_method: str
    task: str
    status: EvidenceUseStatus
    biological_domain: str
    strain_stock_or_population: str
    genome_build_or_target_sequence: str
    edit_system: str
    delivery_or_developmental_context: str
    labels_public: bool
    target_count: int | None
    sample_count: int | None
    independent_evidence: str
    license_status: str
    training_overlap_audited: bool
    source_reference: str
    limitations: str

    @property
    def eligible_for_predictive_capability(self) -> bool:
        """Return whether this record supports a bounded executable adapter."""
        return self.status is EvidenceUseStatus.USABLE_ADAPTER


@dataclass(frozen=True)
class SpeciesReadiness:
    species_profile: str
    scientific_name: str
    evidence_records: tuple[EvidenceQualification, ...]
    predictive_adapter_available: bool
    interpretation: str


_RECORDS = (
    EvidenceQualification(
        record_id="mouse-indelphi-mesc-v1",
        species_profile="mouse",
        predictor_or_method="inDelphi",
        task="repair_outcome",
        status=EvidenceUseStatus.USABLE_ADAPTER,
        biological_domain="mESC integrated-target model",
        strain_stock_or_population="C57BL/6JJmsSlc transfer study",
        genome_build_or_target_sequence="target-sequence model; dossier bound to GRCm39",
        edit_system="SpCas9 knockout",
        delivery_or_developmental_context="mESC model with mouse blastocyst transfer evidence",
        labels_public=True,
        target_count=14,
        sample_count=1182,
        independent_evidence="retrospective mouse-embryo outcome-pair comparison",
        license_status="restricted upstream model; external-result import only",
        training_overlap_audited=True,
        source_reference="https://doi.org/10.1038/s42003-026-09771-z",
        limitations="not editing efficiency, phenotype, welfare, or safety prediction",
    ),
    EvidenceQualification(
        record_id="zebrafish-crisprscan-transfer-v1",
        species_profile="zebrafish",
        predictor_or_method="CRISPRscan",
        task="guide_activity",
        status=EvidenceUseStatus.USABLE_ADAPTER,
        biological_domain="in-vivo zebrafish embryo mutagenesis",
        strain_stock_or_population="predominantly TU; independent RNP study",
        genome_build_or_target_sequence="sequence model; dossier bound to GRCz12tu",
        edit_system="SpCas9 knockout",
        delivery_or_developmental_context="T7 guide training; RNP transfer benchmark",
        labels_public=True,
        target_count=50,
        sample_count=50,
        independent_evidence="gene-disjoint training overlap audit and RNP transfer benchmark",
        license_status="adapter coefficients under upstream MIT notice",
        training_overlap_audited=True,
        source_reference="https://pmc.ncbi.nlm.nih.gov/articles/PMC4589495/",
        limitations="RNP transfer evidence does not establish in-domain calibration",
    ),
    EvidenceQualification(
        record_id="fruit-fly-housden-s2r-v1",
        species_profile="fruit_fly",
        predictor_or_method="Housden",
        task="guide_activity_ranking",
        status=EvidenceUseStatus.USABLE_ADAPTER,
        biological_domain="Drosophila S2R+ cell culture",
        strain_stock_or_population="S2R+ cell line; reference context ISO-1",
        genome_build_or_target_sequence="20-nt protospacer sequence model",
        edit_system="SpCas9 knockout",
        delivery_or_developmental_context="cell-culture sgRNA activity assay",
        labels_public=False,
        target_count=75,
        sample_count=75,
        independent_evidence="reported comparison with three prior Drosophila studies",
        license_status=(
            "official remote service; coefficients and restricted upstream code "
            "are not redistributed"
        ),
        training_overlap_audited=False,
        source_reference="https://pmc.ncbi.nlm.nih.gov/articles/PMC4642709/",
        limitations="in vivo and embryo/germline predictive validity are not established",
    ),
    EvidenceQualification(
        record_id="rat-embryo-protocol-gap-v1",
        species_profile="rat",
        predictor_or_method="rat embryo CRISPR protocol",
        task="guide_activity",
        status=EvidenceUseStatus.INSUFFICIENT_PUBLIC_DATA,
        biological_domain="rat embryo editing protocol",
        strain_stock_or_population="not qualified as a multi-guide prediction dataset",
        genome_build_or_target_sequence="not qualified",
        edit_system="CRISPR-Cas9 knockout",
        delivery_or_developmental_context="embryo protocol",
        labels_public=False,
        target_count=None,
        sample_count=None,
        independent_evidence="method evidence only",
        license_status="publication access does not create a reusable labeled benchmark",
        training_overlap_audited=False,
        source_reference="https://www.nature.com/articles/nprot.2014.171",
        limitations="no qualified public multi-target labels for calibration",
    ),
    EvidenceQualification(
        record_id="rhesus-cas9-trio-wgs-v1",
        species_profile="rhesus_macaque",
        predictor_or_method="Cas9-edited rhesus trio deep sequencing",
        task="unexpected_variant_hazard_observation",
        status=EvidenceUseStatus.HAZARD_EVIDENCE_ONLY,
        biological_domain="locus-specific edited rhesus monkeys",
        strain_stock_or_population="study animals; not a population-calibrated cohort",
        genome_build_or_target_sequence="study-specific targets",
        edit_system="Cas9",
        delivery_or_developmental_context="embryo editing followed by trio sequencing",
        labels_public=True,
        target_count=None,
        sample_count=None,
        independent_evidence="small locus-specific observational study",
        license_status="research article evidence; no model artifact",
        training_overlap_audited=False,
        source_reference="https://pmc.ncbi.nlm.nih.gov/articles/PMC6892871/",
        limitations="cannot estimate general guide, animal, or population-level risk",
    ),
    EvidenceQualification(
        record_id="cynomolgus-ccr5-wgs-v1",
        species_profile="cynomolgus_macaque",
        predictor_or_method="CCR5-edited blastomere whole-genome sequencing",
        task="structural_variant_and_off_target_hazard_observation",
        status=EvidenceUseStatus.HAZARD_EVIDENCE_ONLY,
        biological_domain="CCR5-edited Mauritian cynomolgus blastomeres",
        strain_stock_or_population="Mauritian cynomolgus macaque",
        genome_build_or_target_sequence="CCR5 and predicted off-target sites",
        edit_system="CRISPR-Cas9",
        delivery_or_developmental_context="blastomere editing and whole-genome sequencing",
        labels_public=True,
        target_count=1,
        sample_count=2,
        independent_evidence="two-embryo hazard observation",
        license_status="research article evidence; no model artifact",
        training_overlap_audited=False,
        source_reference="https://pmc.ncbi.nlm.nih.gov/articles/PMC9877282/",
        limitations="cannot calibrate a general predictor from two embryos and one target",
    ),
)


def readiness_for_species(species_profile: str) -> SpeciesReadiness:
    """Return qualified evidence without promoting observational hazards."""
    if species_profile not in PROFILES:
        raise ValueError(f"unknown species profile {species_profile!r}.")
    records = tuple(
        record for record in _RECORDS if record.species_profile == species_profile
    )
    available = any(record.eligible_for_predictive_capability for record in records)
    return SpeciesReadiness(
        species_profile=species_profile,
        scientific_name=PROFILES[species_profile].scientific_name,
        evidence_records=records,
        predictive_adapter_available=available,
        interpretation=(
            "Availability is task- and domain-specific. Hazard observations and "
            "method papers never become predictive capability by themselves."
        ),
    )


def readiness_matrix() -> dict[str, SpeciesReadiness]:
    """Return readiness for every registered species in deterministic order."""
    return {
        species_profile: readiness_for_species(species_profile)
        for species_profile in PROFILES
    }
