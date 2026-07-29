"""A file-based, auditable research workflow for edit-impact assessments."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, Mapping

from .behive import (
    BehiveEfficiencyPrediction,
    integrate_behive_efficiency,
    normalize_behive_efficiency,
)
from .behive_bystander import BehiveBystanderPrediction, normalize_behive_bystander
from .edit_assessment import EditEvidence, assess_edit
from .predictors import PredictionTask, PredictorOutput, integrate_outputs
from .provenance import StudyContext, create_record
from .species import validate_study_context


REQUIRED_CONTEXT = (
    "species",
    "strain_or_breed",
    "genome_build",
    "edit_class",
    "evidence_snapshot",
)
REQUIRED_EVIDENCE = (
    "on_target_uncertainty",
    "off_target_evidence",
    "network_impact_evidence",
    "welfare_relevance",
)
DEFAULT_MODEL_VERSION = "0.6.0"


def assess_request(
    request: Mapping[str, Any], model_version: str = DEFAULT_MODEL_VERSION
) -> dict[str, Any]:
    """Validate a standard request and return a JSON-serializable report.

    Inputs are evidence summaries produced under the user's approved research
    workflow. This function does not accept or generate edit-design details.
    """
    context_data = _mapping(request, "study_context")
    evidence_data = _mapping(request, "evidence")
    _required(context_data, REQUIRED_CONTEXT, "study_context")
    _required(evidence_data, REQUIRED_EVIDENCE, "evidence")

    context = StudyContext(**{key: context_data[key] for key in REQUIRED_CONTEXT})
    species_validation = validate_study_context(context)
    if species_validation.errors:
        raise ValueError("invalid study context: " + " ".join(species_validation.errors))
    evidence = EditEvidence(**{key: evidence_data[key] for key in REQUIRED_EVIDENCE})
    record = create_record(context, assess_edit(evidence), model_version)
    result = asdict(record)
    result["assessment"]["tier"] = record.assessment.tier.value
    result["species_validation"] = {
        "supported": species_validation.supported,
        "profile_key": species_validation.profile_key,
        "warnings": list(species_validation.warnings),
    }
    outputs = _predictor_outputs(request.get("predictor_outputs", []))
    integrated = integrate_outputs(outputs, context.species, context.edit_class)
    result["predictor_outputs"] = [
        {
            "predictor": item.output.predictor,
            "predictor_version": item.output.predictor_version,
            "task": item.output.task.value,
            "concern_score": item.output.concern_score,
            "confidence": item.output.confidence,
            "applicability": item.applicability.value,
            "note": item.note,
            "evidence_reference": item.output.evidence_reference,
        }
        for item in integrated
    ]
    behive_outputs = _behive_outputs(request.get("behive_efficiency_outputs", []))
    behive_integrated = integrate_behive_efficiency(
        behive_outputs, context.species, context.edit_class
    )
    result["model_predictions"] = []
    for item in behive_integrated:
        prediction = asdict(item.prediction)
        prediction["applicability"] = item.applicability
        prediction["applicability_note"] = item.note
        result["model_predictions"].append(prediction)
    bystander_outputs = _behive_bystander_outputs(
        request.get("behive_bystander_outputs", [])
    )
    bystander_matches = (
        context.species.casefold() in {"mouse", "mus musculus"}
        and context.edit_class.casefold().replace("-", "_").replace(" ", "_")
        in {"base_editing", "base_editor"}
    )
    for output in bystander_outputs:
        prediction = asdict(output)
        prediction["applicability"] = (
            "declared_match" if bystander_matches else "out_of_scope"
        )
        prediction["applicability_note"] = (
            "Matches the declared mouse/base-editing scope, limited to mES culture."
            if bystander_matches
            else "Does not match the declared mouse and base-editing scope; not used as applicable evidence."
        )
        result["model_predictions"].append(prediction)
    result["report_notice"] = (
        "Research decision-support only. This report does not establish safety, "
        "authorize an edit, or replace ethics, biosafety, veterinary, or experimental review."
    )
    return result


def _mapping(request: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = request.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object.")
    return value


def _required(data: Mapping[str, Any], fields: tuple[str, ...], section: str) -> None:
    missing = [field for field in fields if field not in data]
    if missing:
        raise ValueError(f"{section} is missing required fields: {', '.join(missing)}")


def _predictor_outputs(raw_outputs: Any) -> tuple[PredictorOutput, ...]:
    if not isinstance(raw_outputs, list):
        raise ValueError("predictor_outputs must be a list.")
    outputs: list[PredictorOutput] = []
    for raw in raw_outputs:
        if not isinstance(raw, Mapping):
            raise ValueError("each predictor output must be an object.")
        try:
            outputs.append(
                PredictorOutput(
                    predictor=str(raw["predictor"]),
                    predictor_version=str(raw["predictor_version"]),
                    task=PredictionTask(raw["task"]),
                    concern_score=float(raw["concern_score"]),
                    confidence=float(raw["confidence"]),
                    supported_species=tuple(raw["supported_species"]),
                    supported_edit_classes=tuple(raw["supported_edit_classes"]),
                    evidence_reference=str(raw["evidence_reference"]),
                )
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"invalid predictor output: {error}") from error
    return tuple(outputs)


def _behive_outputs(raw_outputs: Any) -> tuple[BehiveEfficiencyPrediction, ...]:
    if not isinstance(raw_outputs, list):
        raise ValueError("behive_efficiency_outputs must be a list.")
    outputs: list[BehiveEfficiencyPrediction] = []
    for raw in raw_outputs:
        if not isinstance(raw, Mapping):
            raise ValueError("each BE-Hive efficiency output must be an object.")
        try:
            outputs.append(normalize_behive_efficiency(raw))
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid BE-Hive efficiency output: {error}") from error
    return tuple(outputs)


def _behive_bystander_outputs(
    raw_outputs: Any,
) -> tuple[BehiveBystanderPrediction, ...]:
    if not isinstance(raw_outputs, list):
        raise ValueError("behive_bystander_outputs must be a list.")
    outputs: list[BehiveBystanderPrediction] = []
    for raw in raw_outputs:
        if not isinstance(raw, Mapping):
            raise ValueError("each BE-Hive bystander output must be an object.")
        try:
            outputs.append(normalize_behive_bystander(raw))
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid BE-Hive bystander output: {error}") from error
    return tuple(outputs)
