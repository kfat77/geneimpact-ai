"""Minimal immutable metadata required to audit a research assessment."""

from __future__ import annotations

from dataclasses import dataclass

from .edit_assessment import EditAssessment


@dataclass(frozen=True)
class StudyContext:
    """Declared applicability boundary for one assessment."""

    species: str
    strain_or_breed: str
    genome_build: str
    edit_class: str
    evidence_snapshot: str


@dataclass(frozen=True)
class AssessmentRecord:
    """A result that is meaningful only with its study context and code version."""

    context: StudyContext
    assessment: EditAssessment
    model_version: str


def create_record(
    context: StudyContext, assessment: EditAssessment, model_version: str
) -> AssessmentRecord:
    """Create an auditable record and reject missing applicability metadata."""
    values = {**vars(context), "model_version": model_version}
    missing = [name for name, value in values.items() if not value or not value.strip()]
    if missing:
        raise ValueError(f"missing required metadata: {', '.join(missing)}")
    return AssessmentRecord(context=context, assessment=assessment, model_version=model_version)
