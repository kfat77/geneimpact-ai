"""Leakage-aware validation of BE-Hive efficiency audit records."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from math import isfinite, sqrt
from pathlib import Path
import re
from typing import Any, Mapping

from .behive import BEHIVE_EFFICIENCY_COMMIT, BEHIVE_MOUSE_EDITORS


_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class BehiveValidationMetrics:
    record_count: int
    pearson_logit_vs_observed: float
    calibrated_mae: float | None
    calibrated_rmse: float | None


@dataclass(frozen=True)
class BehiveValidationReport:
    dataset_name: str
    source_reference: str
    species: str
    cell_type: str
    base_editor: str
    model_commit: str
    independent_lab: bool
    training_sequence_overlap_count: int
    training_sequence_overlap_method: str
    input_sha256: str
    metrics: BehiveValidationMetrics
    interpretation: str


def evaluate_behive_validation(path: Path) -> BehiveValidationReport:
    """Evaluate one independently declared, single-editor mES dataset."""
    raw_bytes = path.read_bytes()
    try:
        document = json.loads(raw_bytes)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid validation JSON: {error}") from error
    if not isinstance(document, Mapping):
        raise ValueError("validation document must be an object.")
    dataset = _mapping(document, "dataset")
    records = document.get("records")
    if not isinstance(records, list):
        raise ValueError("records must be a list.")

    declaration = _validate_declaration(dataset)
    observed: list[float] = []
    logits: list[float] = []
    calibrated: list[float | None] = []
    sequence_hashes: set[str] = set()
    for index, raw in enumerate(records):
        if not isinstance(raw, Mapping):
            raise ValueError(f"record {index} must be an object.")
        sequence_hash = str(raw.get("sequence_sha256", ""))
        if not _SHA256_PATTERN.fullmatch(sequence_hash):
            raise ValueError(f"record {index} has an invalid sequence_sha256.")
        if sequence_hash in sequence_hashes:
            raise ValueError(f"record {index} duplicates a sequence_sha256.")
        sequence_hashes.add(sequence_hash)
        if raw.get("base_editor") != declaration["base_editor"]:
            raise ValueError(f"record {index} does not match the declared base editor.")
        if raw.get("model_commit") != BEHIVE_EFFICIENCY_COMMIT:
            raise ValueError(f"record {index} does not use the verified model commit.")
        observed.append(_bounded(raw.get("observed_fraction"), f"record {index} observed_fraction"))
        logits.append(_finite(raw.get("predicted_logit_score"), f"record {index} predicted_logit_score"))
        calibrated.append(
            _optional_bounded(raw.get("calibrated_fraction"), f"record {index} calibrated_fraction")
        )

    if len(records) < 3:
        raise ValueError("an independent validation dataset must contain at least three records.")
    if any(value is None for value in calibrated) and any(value is not None for value in calibrated):
        raise ValueError("calibrated_fraction must be present for all records or none.")

    correlation = _pearson(logits, observed)
    calibrated_values = [value for value in calibrated if value is not None]
    mae = None
    rmse = None
    if calibrated_values:
        errors = [prediction - truth for prediction, truth in zip(calibrated_values, observed)]
        mae = sum(abs(error) for error in errors) / len(errors)
        rmse = sqrt(sum(error * error for error in errors) / len(errors))

    return BehiveValidationReport(
        dataset_name=declaration["dataset_name"],
        source_reference=declaration["source_reference"],
        species="mouse",
        cell_type="mES",
        base_editor=declaration["base_editor"],
        model_commit=BEHIVE_EFFICIENCY_COMMIT,
        independent_lab=True,
        training_sequence_overlap_count=0,
        training_sequence_overlap_method=declaration["training_sequence_overlap_method"],
        input_sha256=sha256(raw_bytes).hexdigest(),
        metrics=BehiveValidationMetrics(
            record_count=len(records),
            pearson_logit_vs_observed=correlation,
            calibrated_mae=mae,
            calibrated_rmse=rmse,
        ),
        interpretation=(
            "Performance applies only to the declared independent mES dataset and editor. "
            "It does not establish whole-animal predictive validity or safety."
        ),
    )


def _validate_declaration(dataset: Mapping[str, Any]) -> dict[str, str]:
    required = (
        "dataset_name",
        "source_reference",
        "species",
        "cell_type",
        "base_editor",
        "model_commit",
        "independent_lab",
        "training_sequence_overlap_count",
        "training_sequence_overlap_method",
    )
    missing = [key for key in required if key not in dataset]
    if missing:
        raise ValueError(f"dataset is missing required fields: {', '.join(missing)}")
    dataset_name = str(dataset["dataset_name"]).strip()
    source_reference = str(dataset["source_reference"]).strip()
    overlap_method = str(dataset["training_sequence_overlap_method"]).strip()
    if not dataset_name or not source_reference or not overlap_method:
        raise ValueError(
            "dataset_name, source_reference, and training_sequence_overlap_method are required."
        )
    if str(dataset["species"]).casefold() not in {"mouse", "mus musculus"}:
        raise ValueError("BE-Hive mES validation requires a mouse dataset.")
    if dataset["cell_type"] != "mES":
        raise ValueError("BE-Hive mES validation requires cell_type mES.")
    base_editor = str(dataset["base_editor"])
    if base_editor not in BEHIVE_MOUSE_EDITORS:
        raise ValueError("the declared base editor is not supported by the verified mES model.")
    if dataset["model_commit"] != BEHIVE_EFFICIENCY_COMMIT:
        raise ValueError("the validation dataset must use the verified BE-Hive model commit.")
    if dataset["independent_lab"] is not True:
        raise ValueError("independent_lab must be true for an independent validation report.")
    overlap_count = dataset["training_sequence_overlap_count"]
    if isinstance(overlap_count, bool) or not isinstance(overlap_count, int):
        raise ValueError("training_sequence_overlap_count must be an integer.")
    if overlap_count != 0:
        raise ValueError("training_sequence_overlap_count must be zero.")
    return {
        "dataset_name": dataset_name,
        "source_reference": source_reference,
        "base_editor": base_editor,
        "training_sequence_overlap_method": overlap_method,
    }


def _mapping(document: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = document.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object.")
    return value


def _finite(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric.") from error
    if not isfinite(number):
        raise ValueError(f"{label} must be finite.")
    return number


def _bounded(value: Any, label: str) -> float:
    number = _finite(value, label)
    if not 0 <= number <= 1:
        raise ValueError(f"{label} must be between 0 and 1.")
    return number


def _optional_bounded(value: Any, label: str) -> float | None:
    if value is None:
        return None
    return _bounded(value, label)


def _pearson(left: list[float], right: list[float]) -> float:
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    left_delta = [value - left_mean for value in left]
    right_delta = [value - right_mean for value in right]
    denominator = sqrt(
        sum(value * value for value in left_delta)
        * sum(value * value for value in right_delta)
    )
    if denominator == 0:
        raise ValueError("Pearson correlation requires variation in predictions and observations.")
    return sum(a * b for a, b in zip(left_delta, right_delta)) / denominator
