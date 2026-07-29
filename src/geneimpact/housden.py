"""Audit normalization for externally obtained Drosophila Housden scores."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from hashlib import sha256
import math
import re
from typing import Any, Mapping


HOUSDEN_METHOD_REFERENCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4642709/"
HOUSDEN_SERVICE_URL = "https://www.flyrnai.org/evaluateCrispr/"
HOUSDEN_HELP_URL = "https://www.flyrnai.org/evaluateCrispr/help.jsp"
HOUSDEN_TRAINING_DOMAIN = "Drosophila S2R+ cell culture"
HOUSDEN_PUBLISHED_HIGH_EFFICIENCY_THRESHOLD = 7.5
HOUSDEN_CURRENT_RECOMMENDED_THRESHOLD = 5.0
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9_.:-]{1,80}$")
_PROTOSPACER_PATTERN = re.compile(r"^[ACGT]{20}$")
_ACCEPTED_BUILD_NAMES = {
    "release 6 plus iso1 mt",
    "bdgp6.54",
    "dm6",
}


@dataclass(frozen=True)
class HousdenPrediction:
    predictor: str
    guide_id: str
    species_profile: str
    genome_build: str
    assembly_accession: str
    sequence_source_strain_or_isolate: str
    sequence_sha256: str
    nuclease: str
    guide_expression: str
    developmental_context: str
    biological_domain: str
    housden_score: float
    score_semantics: str
    within_current_service_reference_range: bool
    published_high_efficiency_threshold: float
    current_recommended_threshold: float
    method_reference: str
    help_reference: str
    source: str
    source_url: str
    retrieved_at: str
    source_response_sha256: str
    source_document_sha256: str | None
    service_version_status: str
    warnings: tuple[str, ...]


def normalize_housden(
    document: Mapping[str, Any],
    *,
    source_document_sha256: str | None = None,
) -> HousdenPrediction:
    """Validate an official FlyRNAi score envelope and remove the raw sequence."""
    request = _mapping(document, "request")
    execution = _mapping(document, "execution")
    raw_output = _mapping(document, "raw_output")

    guide_id = str(request.get("guide_id", ""))
    if not _IDENTIFIER_PATTERN.fullmatch(guide_id):
        raise ValueError("request.guide_id must be a safe 1-80 character identifier.")
    if request.get("species_profile") != "fruit_fly":
        raise ValueError("Housden import is available only for species_profile fruit_fly.")
    genome_build = str(request.get("genome_build", ""))
    if genome_build.casefold() not in _ACCEPTED_BUILD_NAMES:
        raise ValueError("request.genome_build must match the registered fruit-fly build.")
    if request.get("assembly_accession") != "GCF_000001215.4":
        raise ValueError(
            "request.assembly_accession must be GCF_000001215.4."
        )
    strain = str(request.get("sequence_source_strain_or_isolate", ""))
    if strain.casefold() != "iso-1":
        raise ValueError(
            "request.sequence_source_strain_or_isolate must be ISO-1."
        )
    protospacer = str(request.get("protospacer", "")).upper()
    if not _PROTOSPACER_PATTERN.fullmatch(protospacer):
        raise ValueError("request.protospacer must contain exactly 20 A/C/G/T bases.")
    if request.get("nuclease") != "SpCas9":
        raise ValueError("request.nuclease must be SpCas9.")
    if request.get("guide_expression") != "u6_sgrna":
        raise ValueError("request.guide_expression must be u6_sgrna.")
    developmental_context = str(request.get("developmental_context", ""))
    if developmental_context != "drosophila_s2r_plus_cell_culture":
        raise ValueError(
            "Housden scores are accepted only in the declared Drosophila S2R+ "
            "cell-culture domain."
        )

    if execution.get("source") != "flyrnai_evaluate_crispr":
        raise ValueError("execution.source must be flyrnai_evaluate_crispr.")
    if execution.get("source_url") != HOUSDEN_SERVICE_URL:
        raise ValueError(
            f"execution.source_url must be the official service {HOUSDEN_SERVICE_URL!r}."
        )
    retrieved_at = _timestamp(execution.get("retrieved_at"))
    source_response_sha256 = _digest(
        execution.get("source_response_sha256"),
        "execution.source_response_sha256",
    )
    if source_document_sha256 is not None:
        source_document_sha256 = _digest(
            source_document_sha256,
            "source_document_sha256",
        )

    score = _score(raw_output.get("housden_score"))
    return HousdenPrediction(
        predictor="Housden",
        guide_id=guide_id,
        species_profile="fruit_fly",
        genome_build=genome_build,
        assembly_accession="GCF_000001215.4",
        sequence_source_strain_or_isolate=strain,
        sequence_sha256=sha256(protospacer.encode("ascii")).hexdigest(),
        nuclease="SpCas9",
        guide_expression="u6_sgrna",
        developmental_context=developmental_context,
        biological_domain=HOUSDEN_TRAINING_DOMAIN,
        housden_score=score,
        score_semantics="ranking_score_not_probability",
        within_current_service_reference_range=1.47 <= score <= 12.32,
        published_high_efficiency_threshold=(
            HOUSDEN_PUBLISHED_HIGH_EFFICIENCY_THRESHOLD
        ),
        current_recommended_threshold=HOUSDEN_CURRENT_RECOMMENDED_THRESHOLD,
        method_reference=HOUSDEN_METHOD_REFERENCE,
        help_reference=HOUSDEN_HELP_URL,
        source="flyrnai_evaluate_crispr",
        source_url=HOUSDEN_SERVICE_URL,
        retrieved_at=retrieved_at,
        source_response_sha256=source_response_sha256,
        source_document_sha256=source_document_sha256,
        service_version_status="live_service_unversioned",
        warnings=(
            "The live service does not expose a version; retain the response checksum.",
            "The Housden score is a ranking score, not an editing probability.",
            "Published and current service guidance use different score thresholds.",
            "The score does not evaluate off-target effects.",
            "Embryo, germline, in-vivo, phenotype, and safety validity are not established.",
        ),
    )


def _mapping(document: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = document.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object.")
    return value


def _digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or not _SHA256_PATTERN.fullmatch(value):
        raise ValueError(f"{label} must be a lowercase SHA-256 digest.")
    return value


def _timestamp(value: Any) -> str:
    if not isinstance(value, str):
        raise ValueError("execution.retrieved_at must be an ISO-8601 timestamp.")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(
            "execution.retrieved_at must be an ISO-8601 timestamp."
        ) from error
    if parsed.tzinfo is None:
        raise ValueError("execution.retrieved_at must include a timezone.")
    return value


def _score(value: Any) -> float:
    if isinstance(value, bool):
        raise ValueError("raw_output.housden_score must be numeric.")
    try:
        score = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("raw_output.housden_score must be numeric.") from error
    if not math.isfinite(score) or not 0.0 <= score <= 20.0:
        raise ValueError("raw_output.housden_score must be between 0 and 20.")
    return score
