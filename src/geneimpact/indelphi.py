"""Strict audit import for externally executed inDelphi mESC predictions."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from math import isfinite, log
import re
from typing import Any, Mapping, Sequence

from .indelphi_validation import (
    IndelphiMouseTransferEvidence,
    indelphi_mouse_transfer_evidence,
)
from .species import MOUSE_PROFILE


INDELPHI_REPOSITORY = "https://github.com/maxwshen/inDelphi-model"
INDELPHI_COMMIT = "9ab67ca53ebb91e49aeb4530ec1e999ee9827ca1"
INDELPHI_MODEL_ARTIFACT_FAMILY = "model-sklearn-0.18.1"
INDELPHI_SKLEARN_VERSION = "0.18.1"
INDELPHI_REFERENCE = f"{INDELPHI_REPOSITORY}/tree/{INDELPHI_COMMIT}"
MAX_IMPORTED_OUTCOMES = 10_000
MAX_REPORTED_OUTCOMES = 25
MIN_FLANK_NT = 60
MAX_SEQUENCE_NT = 10_000
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9_.:-]{1,80}$")
_STAT_TOLERANCE = 1e-6


@dataclass(frozen=True)
class IndelphiOutcome:
    category: str
    length: int
    genotype_position: int | str | None
    inserted_base: str | None
    representation: str
    predicted_frequency_percent: float
    genotype_sha256: str | None


@dataclass(frozen=True)
class IndelphiPrediction:
    predictor: str
    predictor_version: str
    predictor_commit: str
    model_artifact_family: str
    task: str
    target_id: str
    species_profile: str
    genome_build: str
    assembly_accession: str
    sequence_source_strain_or_isolate: str
    cell_type: str
    biological_context: str
    nuclease: str
    delivery_context: str
    developmental_context: str
    target_context: str
    sequence_sha256: str
    cutsite: int
    upstream_context_nt: int
    downstream_context_nt: int
    frequency_semantics: str
    total_predicted_frequency_percent: float
    precision: float
    phi: float
    microhomology_strength_log_phi: float
    insertion_frequency_percent: float
    microhomology_deletion_frequency_percent: float
    microhomologyless_deletion_frequency_percent: float
    frameshift_frequency_percent: float | None
    frame_0_frequency_percent: float | None
    frame_1_frequency_percent: float | None
    frame_2_frequency_percent: float | None
    highest_outcome_frequency_percent: float
    expected_indel_length: float
    outcome_count: int
    reported_outcome_count: int
    reported_frequency_percent: float
    top_outcomes: tuple[IndelphiOutcome, ...]
    raw_output_sha256: str
    source_document_sha256: str | None
    model_bundle_sha256: str
    python_version: str
    sklearn_version: str
    evidence_reference: str
    external_validation: IndelphiMouseTransferEvidence
    source_license_status: str
    execution_status: str
    warnings: tuple[str, ...]


def normalize_indelphi(
    document: Mapping[str, Any],
    *,
    source_document_sha256: str | None = None,
) -> IndelphiPrediction:
    """Validate one external inDelphi result without redistributing its model."""
    request = _mapping(document, "request")
    execution = _mapping(document, "execution")
    raw_output = _mapping(document, "raw_output")
    stats = _mapping(raw_output, "stats")
    raw_outcomes = raw_output.get("outcomes")
    if (
        not isinstance(raw_outcomes, Sequence)
        or isinstance(raw_outcomes, (str, bytes))
        or not 1 <= len(raw_outcomes) <= MAX_IMPORTED_OUTCOMES
    ):
        raise ValueError(
            f"raw_output.outcomes must contain 1-{MAX_IMPORTED_OUTCOMES} objects."
        )
    if source_document_sha256 is not None and not _SHA256_PATTERN.fullmatch(
        source_document_sha256
    ):
        raise ValueError("source_document_sha256 must be a lowercase SHA-256 digest.")

    declaration = _validate_request(request)
    execution_record = _validate_execution(execution)
    if execution_record["cell_type"] != declaration["cell_type"]:
        raise ValueError("execution cell_type does not match the request.")

    outcomes = tuple(
        _normalize_outcome(raw, index)
        for index, raw in enumerate(raw_outcomes)
    )
    _reject_duplicate_outcomes(outcomes)
    calculated = _calculated_statistics(outcomes)
    reported = _validate_stats(stats, declaration, calculated)
    ordered = tuple(
        sorted(
            outcomes,
            key=lambda item: (
                -item.predicted_frequency_percent,
                item.category,
                item.length,
                str(item.genotype_position),
                str(item.inserted_base),
            ),
        )
    )
    top_outcomes = ordered[:MAX_REPORTED_OUTCOMES]
    frameshift_fields = (
        (
            reported["frameshift"],
            reported["frame_0"],
            reported["frame_1"],
            reported["frame_2"],
        )
        if declaration["target_context"] == "coding_sequence"
        else (None, None, None, None)
    )
    raw_digest = sha256(
        json.dumps(
            raw_output,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()

    return IndelphiPrediction(
        predictor="inDelphi",
        predictor_version=INDELPHI_COMMIT,
        predictor_commit=INDELPHI_COMMIT,
        model_artifact_family=INDELPHI_MODEL_ARTIFACT_FAMILY,
        task="cas9_repair_outcome_distribution",
        target_id=declaration["target_id"],
        species_profile="mouse",
        genome_build=declaration["genome_build"],
        assembly_accession=declaration["assembly_accession"],
        sequence_source_strain_or_isolate=declaration[
            "sequence_source_strain_or_isolate"
        ],
        cell_type="mESC",
        biological_context=(
            "mESC model trained on integrated 55-bp human genomic fragments"
        ),
        nuclease="SpCas9",
        delivery_context=declaration["delivery_context"],
        developmental_context=declaration["developmental_context"],
        target_context=declaration["target_context"],
        sequence_sha256=sha256(declaration["sequence"].encode("ascii")).hexdigest(),
        cutsite=declaration["cutsite"],
        upstream_context_nt=declaration["cutsite"],
        downstream_context_nt=(
            len(declaration["sequence"]) - declaration["cutsite"]
        ),
        frequency_semantics=(
            "conditional_percent_distribution_among_modeled_edited_products"
        ),
        total_predicted_frequency_percent=calculated["total"],
        precision=reported["precision"],
        phi=reported["phi"],
        microhomology_strength_log_phi=log(reported["phi"]),
        insertion_frequency_percent=calculated["insertion"],
        microhomology_deletion_frequency_percent=calculated["mh_deletion"],
        microhomologyless_deletion_frequency_percent=calculated[
            "mhless_deletion"
        ],
        frameshift_frequency_percent=frameshift_fields[0],
        frame_0_frequency_percent=frameshift_fields[1],
        frame_1_frequency_percent=frameshift_fields[2],
        frame_2_frequency_percent=frameshift_fields[3],
        highest_outcome_frequency_percent=calculated["highest_outcome"],
        expected_indel_length=calculated["expected_indel_length"],
        outcome_count=len(outcomes),
        reported_outcome_count=len(top_outcomes),
        reported_frequency_percent=sum(
            item.predicted_frequency_percent for item in top_outcomes
        ),
        top_outcomes=top_outcomes,
        raw_output_sha256=raw_digest,
        source_document_sha256=source_document_sha256,
        model_bundle_sha256=execution_record["model_bundle_sha256"],
        python_version=execution_record["python_version"],
        sklearn_version=INDELPHI_SKLEARN_VERSION,
        evidence_reference=INDELPHI_REFERENCE,
        external_validation=indelphi_mouse_transfer_evidence(),
        source_license_status=(
            "restricted_noncommercial_research_license; upstream code and models "
            "are not redistributed by GeneImpact AI"
        ),
        execution_status="externally_executed_output_validated_not_recomputed",
        warnings=(
            "Frequencies are conditional on modeled edited products; inDelphi does not predict cutting or editing efficiency.",
            "The model covers 1-bp insertions and deletions up to 59 bp and does not rule out large deletions, rearrangements, integrations, mosaicism, or other unmodeled outcomes.",
            "The mESC training system used integrated short human sequence fragments and does not reproduce endogenous chromatin or a mouse zygote.",
            "External mouse-embryo transfer evidence was moderate overall and heterogeneous by guide; prospective mESC and embryo validation remains required.",
            "Frameshift percentages are reported only for a declared coding-sequence target and remain conditional repair-product frequencies, not knockout or phenotype probabilities.",
            "This output does not establish animal phenotype, welfare impact, or safety.",
        ),
    )


def _validate_request(request: Mapping[str, Any]) -> dict[str, Any]:
    required = (
        "target_id",
        "species_profile",
        "genome_build",
        "assembly_accession",
        "sequence_source_strain_or_isolate",
        "sequence",
        "cutsite",
        "cell_type",
        "model_commit",
        "model_artifact_family",
        "nuclease",
        "delivery_context",
        "developmental_context",
        "target_context",
    )
    missing = [key for key in required if key not in request]
    if missing:
        raise ValueError("inDelphi request is missing fields: " + ", ".join(missing))
    if request["species_profile"] != "mouse":
        raise ValueError("the inDelphi mESC adapter accepts only species_profile mouse.")
    if request["cell_type"] != "mESC":
        raise ValueError("the inDelphi mouse adapter accepts only the mESC model.")
    if request["model_commit"] != INDELPHI_COMMIT:
        raise ValueError(
            f"adapter is verified only for inDelphi commit {INDELPHI_COMMIT}."
        )
    if request["model_artifact_family"] != INDELPHI_MODEL_ARTIFACT_FAMILY:
        raise ValueError(
            "model_artifact_family must be "
            f"{INDELPHI_MODEL_ARTIFACT_FAMILY!r}."
        )
    if request["nuclease"] != "SpCas9":
        raise ValueError("the adapter accepts only canonical SpCas9 repair outcomes.")
    if request["target_context"] not in {
        "coding_sequence",
        "noncoding_or_unknown",
    }:
        raise ValueError(
            "target_context must be coding_sequence or noncoding_or_unknown."
        )
    target_id = str(request["target_id"])
    if not _IDENTIFIER_PATTERN.fullmatch(target_id):
        raise ValueError("target_id must be a safe 1-80 character identifier.")
    profile = MOUSE_PROFILE
    accepted_builds = {
        value.casefold()
        for value in (profile.genome_build, *profile.accepted_build_names)
    }
    genome_build = str(request["genome_build"])
    if genome_build.casefold() not in accepted_builds:
        raise ValueError(
            f"genome_build must match registered {profile.genome_build!r}."
        )
    if request["assembly_accession"] != profile.assembly_accession:
        raise ValueError(
            "assembly_accession must match registered "
            f"{profile.assembly_accession!r}."
        )
    strain = _short_text(
        request["sequence_source_strain_or_isolate"],
        "sequence_source_strain_or_isolate",
    )
    delivery = _short_text(request["delivery_context"], "delivery_context")
    developmental = _short_text(
        request["developmental_context"],
        "developmental_context",
    )
    sequence = str(request["sequence"]).strip().upper()
    if (
        not 2 * MIN_FLANK_NT <= len(sequence) <= MAX_SEQUENCE_NT
        or any(base not in "ACGT" for base in sequence)
    ):
        raise ValueError(
            f"sequence must be {2 * MIN_FLANK_NT}-{MAX_SEQUENCE_NT} A/C/G/T bases."
        )
    cutsite = _integer(request["cutsite"], "cutsite")
    if cutsite < MIN_FLANK_NT or len(sequence) - cutsite < MIN_FLANK_NT:
        raise ValueError(
            f"cutsite requires at least {MIN_FLANK_NT} nt on each side."
        )
    return {
        "target_id": target_id,
        "genome_build": genome_build,
        "assembly_accession": profile.assembly_accession,
        "sequence_source_strain_or_isolate": strain,
        "sequence": sequence,
        "cutsite": cutsite,
        "cell_type": "mESC",
        "delivery_context": delivery,
        "developmental_context": developmental,
        "target_context": request["target_context"],
    }


def _validate_execution(execution: Mapping[str, Any]) -> dict[str, str]:
    required = (
        "repository_commit",
        "model_artifact_family",
        "model_bundle_sha256",
        "cell_type",
        "python_version",
        "sklearn_version",
    )
    missing = [key for key in required if key not in execution]
    if missing:
        raise ValueError(
            "inDelphi execution is missing fields: " + ", ".join(missing)
        )
    if execution["repository_commit"] != INDELPHI_COMMIT:
        raise ValueError("execution repository_commit is not version locked.")
    if execution["model_artifact_family"] != INDELPHI_MODEL_ARTIFACT_FAMILY:
        raise ValueError("execution model_artifact_family is not version locked.")
    if execution["sklearn_version"] != INDELPHI_SKLEARN_VERSION:
        raise ValueError(
            f"execution sklearn_version must be {INDELPHI_SKLEARN_VERSION}."
        )
    if execution["cell_type"] != "mESC":
        raise ValueError("execution cell_type must be mESC.")
    model_digest = str(execution["model_bundle_sha256"])
    if not _SHA256_PATTERN.fullmatch(model_digest):
        raise ValueError("model_bundle_sha256 must be a lowercase SHA-256 digest.")
    python_version = _short_text(execution["python_version"], "python_version")
    return {
        "model_bundle_sha256": model_digest,
        "python_version": python_version,
        "cell_type": "mESC",
    }


def _normalize_outcome(raw: Any, index: int) -> IndelphiOutcome:
    if not isinstance(raw, Mapping):
        raise ValueError(f"outcome {index} must be an object.")
    category = str(raw.get("Category", "")).strip()
    if category not in {"ins", "del"}:
        raise ValueError(f"outcome {index} Category must be ins or del.")
    length = _integer(raw.get("Length"), f"outcome {index} Length")
    if not 1 <= length <= 59:
        raise ValueError(f"outcome {index} Length must be between 1 and 59.")
    frequency = _number(
        raw.get("Predicted frequency"),
        f"outcome {index} Predicted frequency",
        minimum=0,
        maximum=100,
    )
    genotype_position: int | str | None
    inserted_base: str | None
    if category == "ins":
        if length != 1:
            raise ValueError(f"outcome {index} insertion Length must be 1.")
        inserted_base = str(raw.get("Inserted Bases", "")).strip().upper()
        if inserted_base not in {"A", "C", "G", "T"}:
            raise ValueError(
                f"outcome {index} Inserted Bases must be one DNA base."
            )
        genotype_position = None
        representation = "one_bp_insertion_genotype"
    else:
        if raw.get("Inserted Bases") not in {None, ""}:
            raise ValueError(
                f"outcome {index} deletion must not contain Inserted Bases."
            )
        inserted_base = None
        raw_position = raw.get("Genotype position")
        if raw_position == "e":
            genotype_position = "e"
            representation = "microhomologyless_deletion_length_group"
        else:
            genotype_position = _integer(
                raw_position,
                f"outcome {index} Genotype position",
            )
            if not 0 <= genotype_position <= length:
                raise ValueError(
                    f"outcome {index} Genotype position is outside its deletion."
                )
            representation = "microhomology_deletion_genotype"
    raw_genotype = raw.get("Genotype")
    genotype_digest = None
    if raw_genotype not in {None, ""}:
        genotype = str(raw_genotype).strip().upper()
        if (
            not 1 <= len(genotype) <= MAX_SEQUENCE_NT
            or any(base not in "ACGT-" for base in genotype)
        ):
            raise ValueError(f"outcome {index} Genotype is invalid.")
        genotype_digest = sha256(genotype.encode("ascii")).hexdigest()
    return IndelphiOutcome(
        category=category,
        length=length,
        genotype_position=genotype_position,
        inserted_base=inserted_base,
        representation=representation,
        predicted_frequency_percent=frequency,
        genotype_sha256=genotype_digest,
    )


def _reject_duplicate_outcomes(outcomes: Sequence[IndelphiOutcome]) -> None:
    keys = [
        (
            item.category,
            item.length,
            item.genotype_position,
            item.inserted_base,
        )
        for item in outcomes
    ]
    if len(keys) != len(set(keys)):
        raise ValueError("raw_output.outcomes contains duplicate outcome identities.")


def _calculated_statistics(
    outcomes: Sequence[IndelphiOutcome],
) -> dict[str, float]:
    total = sum(item.predicted_frequency_percent for item in outcomes)
    _assert_close(total, 100.0, "outcome frequencies must sum to 100")
    insertion = sum(
        item.predicted_frequency_percent
        for item in outcomes
        if item.category == "ins"
    )
    mh_deletion = sum(
        item.predicted_frequency_percent
        for item in outcomes
        if item.category == "del" and item.genotype_position != "e"
    )
    mhless_deletion = sum(
        item.predicted_frequency_percent
        for item in outcomes
        if item.category == "del" and item.genotype_position == "e"
    )
    frame = {0: 0.0, 1: 0.0, 2: 0.0}
    for item in outcomes:
        shift = item.length % 3 if item.category == "ins" else (-item.length) % 3
        frame[shift] += item.predicted_frequency_percent
    probabilities = [
        item.predicted_frequency_percent / 100
        for item in outcomes
        if item.predicted_frequency_percent > 0
    ]
    entropy = -sum(probability * log(probability) for probability in probabilities)
    precision = (
        1.0
        if len(outcomes) == 1
        else 1 - entropy / log(len(outcomes))
    )
    return {
        "total": total,
        "insertion": insertion,
        "mh_deletion": mh_deletion,
        "mhless_deletion": mhless_deletion,
        "frame_0": frame[0],
        "frame_1": frame[1],
        "frame_2": frame[2],
        "frameshift": frame[1] + frame[2],
        "precision": precision,
        "highest_outcome": max(
            item.predicted_frequency_percent for item in outcomes
        ),
        "highest_deletion": max(
            (
                item.predicted_frequency_percent
                for item in outcomes
                if item.category == "del"
            ),
            default=0.0,
        ),
        "highest_insertion": max(
            (
                item.predicted_frequency_percent
                for item in outcomes
                if item.category == "ins"
            ),
            default=0.0,
        ),
        "expected_indel_length": sum(
            item.length * item.predicted_frequency_percent / 100
            for item in outcomes
        ),
    }


def _validate_stats(
    stats: Mapping[str, Any],
    declaration: Mapping[str, Any],
    calculated: Mapping[str, float],
) -> dict[str, float]:
    if stats.get("Reference sequence") != declaration["sequence"]:
        raise ValueError("stats Reference sequence does not match the request.")
    if stats.get("Cutsite") != declaration["cutsite"]:
        raise ValueError("stats Cutsite does not match the request.")
    if stats.get("Celltype") not in {None, "mESC"}:
        raise ValueError("stats Celltype must be mESC when present.")
    names = {
        "Precision": ("precision", 0, 1),
        "Phi": ("phi", 0, None),
        "1-bp ins frequency": ("insertion", 0, 100),
        "MH del frequency": ("mh_deletion", 0, 100),
        "MHless del frequency": ("mhless_deletion", 0, 100),
        "Frameshift frequency": ("frameshift", 0, 100),
        "Frame +0 frequency": ("frame_0", 0, 100),
        "Frame +1 frequency": ("frame_1", 0, 100),
        "Frame +2 frequency": ("frame_2", 0, 100),
        "Highest outcome frequency": ("highest_outcome", 0, 100),
        "Highest del frequency": ("highest_deletion", 0, 100),
        "Highest ins frequency": ("highest_insertion", 0, 100),
        "Expected indel length": ("expected_indel_length", 0, 59),
    }
    result: dict[str, float] = {}
    for source_name, (key, minimum, maximum) in names.items():
        reported = _number(
            stats.get(source_name),
            f"stats {source_name}",
            minimum=minimum,
            maximum=maximum,
            minimum_exclusive=source_name == "Phi",
        )
        if key != "phi":
            _assert_close(
                reported,
                calculated[key],
                f"stats {source_name} does not match outcomes",
            )
        result[key] = reported
    return result


def _mapping(document: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = document.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object.")
    return value


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer.")
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer.") from error
    if isinstance(value, float) and not value.is_integer():
        raise ValueError(f"{label} must be an integer.")
    if isinstance(value, str) and value.strip() != str(number):
        raise ValueError(f"{label} must be an integer.")
    return number


def _number(
    value: Any,
    label: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    minimum_exclusive: bool = False,
) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric.") from error
    if not isfinite(number):
        raise ValueError(f"{label} must be finite.")
    if minimum is not None and (
        number < minimum or (minimum_exclusive and number == minimum)
    ):
        comparator = "greater than" if minimum_exclusive else "at least"
        raise ValueError(f"{label} must be {comparator} {minimum}.")
    if maximum is not None and number > maximum:
        raise ValueError(f"{label} must be at most {maximum}.")
    return number


def _short_text(value: Any, label: str) -> str:
    text = str(value).strip()
    if not text or len(text) > 160:
        raise ValueError(f"{label} must be a non-empty short string.")
    return text


def _assert_close(actual: float, expected: float, message: str) -> None:
    if abs(actual - expected) > _STAT_TOLERANCE:
        raise ValueError(f"{message}: expected {expected}, received {actual}.")
