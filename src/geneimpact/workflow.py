"""A file-based, auditable research workflow for edit-impact assessments."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, Mapping

from .edit_assessment import EditEvidence, assess_edit
from .provenance import StudyContext, create_record


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


def assess_request(request: Mapping[str, Any], model_version: str = "0.2.0") -> dict[str, Any]:
    """Validate a standard request and return a JSON-serializable report.

    Inputs are evidence summaries produced under the user's approved research
    workflow. This function does not accept or generate edit-design details.
    """
    context_data = _mapping(request, "study_context")
    evidence_data = _mapping(request, "evidence")
    _required(context_data, REQUIRED_CONTEXT, "study_context")
    _required(evidence_data, REQUIRED_EVIDENCE, "evidence")

    context = StudyContext(**{key: context_data[key] for key in REQUIRED_CONTEXT})
    evidence = EditEvidence(**{key: evidence_data[key] for key in REQUIRED_EVIDENCE})
    record = create_record(context, assess_edit(evidence), model_version)
    result = asdict(record)
    result["assessment"]["tier"] = record.assessment.tier.value
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
