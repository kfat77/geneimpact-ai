"""Strict normalization of externally executed BE-Hive bystander predictions."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from math import isfinite
from typing import Any, Mapping


BEHIVE_BYSTANDER_REPOSITORY = "https://github.com/maxwshen/be_predict_bystander"
BEHIVE_BYSTANDER_COMMIT = "31aadd04a25c604857c7592b226ee987e9e20b31"
BEHIVE_BYSTANDER_REFERENCE = (
    f"{BEHIVE_BYSTANDER_REPOSITORY}/tree/{BEHIVE_BYSTANDER_COMMIT}"
)
BEHIVE_BYSTANDER_MOUSE_EDITORS = frozenset(
    {
        "ABE",
        "ABE-CP1040",
        "ABE8",
        "AID",
        "BE4",
        "BE4-CP1028",
        "BE4-H47ES48A",
        "CDA",
        "CG-689",
        "CG-APOBEC1",
        "CG-EE",
        "CG-POLD2-APOBEC1-X",
        "CG-RBMX-eA3A-X",
        "CG-RBMX-eA3A-X-HF",
        "CG-X-689-X-RBMX",
        "CG-X-APOBEC1-X-HF",
        "CG-X-EE-X-X",
        "CG-eA3A-dead",
        "eA3A",
        "eA3A-T31A",
        "eA3A-T31AT44A",
        "evoAPOBEC",
    }
)
MAX_IMPORTED_OUTCOMES = 100_000
MAX_REPORTED_OUTCOMES = 25


@dataclass(frozen=True)
class BehiveBystanderOutcome:
    genotype_sha256: str
    edits: tuple[str, ...]
    predicted_frequency: float


@dataclass(frozen=True)
class BehiveBystanderPrediction:
    predictor: str
    predictor_version: str
    task: str
    base_editor: str
    cell_type: str
    species_scope: tuple[str, ...]
    biological_context: str
    sequence_sha256: str
    frequency_semantics: str
    total_predicted_probability: float
    unreported_probability_mass: float
    outcome_count: int
    reported_outcome_count: int
    reported_outcome_probability: float
    top_outcomes: tuple[BehiveBystanderOutcome, ...]
    raw_output_sha256: str
    evidence_reference: str
    source_license_status: str
    execution_status: str
    warnings: tuple[str, ...]


def normalize_behive_bystander(
    document: Mapping[str, Any],
) -> BehiveBystanderPrediction:
    """Validate a BE-Hive bystander result and retain bounded audit details."""
    request = _mapping(document, "request")
    raw_output = _mapping(document, "raw_output")
    stats = _mapping(raw_output, "stats")
    outcomes = raw_output.get("outcomes")
    if not isinstance(outcomes, list) or not outcomes:
        raise ValueError("raw_output.outcomes must be a non-empty list.")
    if len(outcomes) > MAX_IMPORTED_OUTCOMES:
        raise ValueError(f"raw_output.outcomes may contain at most {MAX_IMPORTED_OUTCOMES} rows.")

    sequence, editor = _validate_request(request)
    _validate_stats(stats, sequence, editor)
    total_probability = _bounded(
        stats.get("Total predicted probability"), "Total predicted probability"
    )
    if total_probability == 0:
        raise ValueError("Total predicted probability must be positive.")

    normalized: list[BehiveBystanderOutcome] = []
    genotype_hashes: set[str] = set()
    for index, raw in enumerate(outcomes):
        if not isinstance(raw, Mapping):
            raise ValueError(f"outcome {index} must be an object.")
        genotype = str(raw.get("Genotype", "")).strip().upper()
        if len(genotype) != 50 or any(base not in "ACGT" for base in genotype):
            raise ValueError(f"outcome {index} Genotype must be a 50-nt DNA sequence.")
        genotype_hash = sha256(genotype.encode("ascii")).hexdigest()
        if genotype_hash in genotype_hashes:
            raise ValueError(f"outcome {index} duplicates a genotype.")
        genotype_hashes.add(genotype_hash)
        frequency = _bounded(raw.get("Predicted frequency"), f"outcome {index} frequency")
        edits = tuple(
            f"{reference}{position}>{alternate}"
            for position, (reference, alternate) in enumerate(
                zip(sequence, genotype), start=-19
            )
            if reference != alternate
        )
        if not edits:
            raise ValueError(f"outcome {index} is an unedited genotype.")
        normalized.append(BehiveBystanderOutcome(genotype_hash, edits, frequency))

    frequency_sum = sum(item.predicted_frequency for item in normalized)
    if abs(frequency_sum - total_probability) > 1e-6:
        raise ValueError(
            "outcome frequencies do not sum to Total predicted probability."
        )
    normalized.sort(key=lambda item: item.predicted_frequency, reverse=True)
    top_outcomes = tuple(normalized[:MAX_REPORTED_OUTCOMES])
    reported_probability = sum(item.predicted_frequency for item in top_outcomes)
    raw_digest = sha256(
        json.dumps(raw_output, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()

    return BehiveBystanderPrediction(
        predictor="BE-Hive bystander",
        predictor_version=BEHIVE_BYSTANDER_COMMIT,
        task="base_editing_bystander_outcomes",
        base_editor=editor,
        cell_type="mES",
        species_scope=("mouse",),
        biological_context="mouse embryonic stem cells (mES)",
        sequence_sha256=sha256(sequence.encode("ascii")).hexdigest(),
        frequency_semantics="conditional_distribution_among_base_editing_outcomes",
        total_predicted_probability=total_probability,
        unreported_probability_mass=max(0.0, 1.0 - total_probability),
        outcome_count=len(normalized),
        reported_outcome_count=len(top_outcomes),
        reported_outcome_probability=reported_probability,
        top_outcomes=top_outcomes,
        raw_output_sha256=raw_digest,
        evidence_reference=BEHIVE_BYSTANDER_REFERENCE,
        source_license_status="not_declared_in_repository",
        execution_status="externally_reported_not_recomputed",
        warnings=(
            "Outcome frequencies are conditional editing-outcome predictions, not probabilities of animal phenotypes or harm.",
            "Only the highest-frequency outcomes are embedded in this bounded audit record; use raw_output_sha256 to bind the full external output.",
            "The model scope is mES culture and does not establish strain-specific, tissue, germline, mosaic, or whole-animal effects.",
        ),
    )


def _validate_request(request: Mapping[str, Any]) -> tuple[str, str]:
    required = ("sequence", "base_editor", "cell_type", "model_commit")
    missing = [key for key in required if key not in request]
    if missing:
        raise ValueError(f"BE-Hive request is missing required fields: {', '.join(missing)}")
    sequence = str(request["sequence"]).strip().upper()
    if len(sequence) != 50 or any(base not in "ACGT" for base in sequence):
        raise ValueError("BE-Hive bystander sequence must be exactly 50 A/C/G/T nucleotides.")
    if sequence[41:43] != "GG":
        raise ValueError("BE-Hive bystander sequence must contain an NGG PAM at positions 21-23.")
    if request["cell_type"] != "mES":
        raise ValueError("the mouse bystander adapter only accepts the mES model.")
    editor = str(request["base_editor"])
    if editor not in BEHIVE_BYSTANDER_MOUSE_EDITORS:
        raise ValueError(f"base editor {editor!r} is not declared for the mES bystander model.")
    if request["model_commit"] != BEHIVE_BYSTANDER_COMMIT:
        raise ValueError(
            f"this adapter is verified only for BE-Hive bystander commit {BEHIVE_BYSTANDER_COMMIT}."
        )
    return sequence, editor


def _validate_stats(stats: Mapping[str, Any], sequence: str, editor: str) -> None:
    if stats.get("50-nt target sequence") != sequence:
        raise ValueError("stats target sequence does not match the declared request.")
    if stats.get("Assumed protospacer sequence") != sequence[20:40]:
        raise ValueError("stats protospacer does not match positions 1-20 of the request.")
    if stats.get("Celltype") != "mES":
        raise ValueError("stats Celltype must be mES.")
    if stats.get("Base editor") != editor:
        raise ValueError("stats Base editor does not match the declared request.")


def _mapping(document: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = document.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object.")
    return value


def _bounded(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric.") from error
    if not isfinite(number) or not 0 <= number <= 1:
        raise ValueError(f"{label} must be finite and between 0 and 1.")
    return number
