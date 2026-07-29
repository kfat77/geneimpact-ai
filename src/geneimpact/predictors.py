"""Normalized ingestion of validated external genome-edit predictor outputs.

Adapters deliberately consume *reported* model outputs rather than making
undocumented remote calls. This keeps provenance, model version, and declared
applicability visible to the researcher.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable


class PredictionTask(str, Enum):
    GUIDE_ACTIVITY = "guide_activity"
    OFF_TARGET = "off_target"
    REPAIR_OUTCOME = "repair_outcome"
    BASE_EDITING = "base_editing"
    PRIME_EDITING = "prime_editing"


class Applicability(str, Enum):
    DECLARED_MATCH = "declared_match"
    OUT_OF_SCOPE = "out_of_scope"


@dataclass(frozen=True)
class PredictorOutput:
    """A normalized, bounded result from an upstream prediction model."""

    predictor: str
    predictor_version: str
    task: PredictionTask
    concern_score: float
    confidence: float
    supported_species: tuple[str, ...]
    supported_edit_classes: tuple[str, ...]
    evidence_reference: str


@dataclass(frozen=True)
class IntegratedPrediction:
    """A single upstream result plus its applicability to the study context."""

    output: PredictorOutput
    applicability: Applicability
    note: str


def integrate_outputs(
    outputs: Iterable[PredictorOutput], species: str, edit_class: str
) -> tuple[IntegratedPrediction, ...]:
    """Label each output as applicable or out of scope without extrapolation."""
    integrated: list[IntegratedPrediction] = []
    for output in outputs:
        _validate(output)
        matches_species = species.casefold() in {item.casefold() for item in output.supported_species}
        matches_edit = edit_class.casefold() in {item.casefold() for item in output.supported_edit_classes}
        if matches_species and matches_edit:
            integrated.append(
                IntegratedPrediction(output, Applicability.DECLARED_MATCH, "Matches declared predictor scope.")
            )
        else:
            integrated.append(
                IntegratedPrediction(
                    output,
                    Applicability.OUT_OF_SCOPE,
                    "Does not match the predictor's declared species and edit-class scope; not used for scoring.",
                )
            )
    return tuple(integrated)


def _validate(output: PredictorOutput) -> None:
    if not output.predictor.strip() or not output.predictor_version.strip() or not output.evidence_reference.strip():
        raise ValueError("predictor, predictor_version, and evidence_reference are required.")
    if not output.supported_species or not output.supported_edit_classes:
        raise ValueError("predictor scope must declare supported species and edit classes.")
    for label, value in {"concern_score": output.concern_score, "confidence": output.confidence}.items():
        if not 0 <= value <= 1:
            raise ValueError(f"{label} must be between 0 and 1.")
