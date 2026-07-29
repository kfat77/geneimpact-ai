"""Transparent research triage for proposed animal genome-edit consequences.

The module maps evidence and uncertainty to a review tier. It does not design
edits, declare an edit safe, or replace experimental and welfare review.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ReviewTier(str, Enum):
    STANDARD_REVIEW = "standard_review"
    ENHANCED_REVIEW = "enhanced_review"
    HIGH_CONCERN_REVIEW = "high_concern_review"


@dataclass(frozen=True)
class EditEvidence:
    """Bounded evidence inputs for a non-operational edit impact assessment."""

    on_target_uncertainty: float
    off_target_evidence: float
    network_impact_evidence: float
    welfare_relevance: float


@dataclass(frozen=True)
class EditAssessment:
    """Traceable prioritization result for human review."""

    concern_score: float
    tier: ReviewTier
    rationale: tuple[str, ...]


def assess_edit(evidence: EditEvidence) -> EditAssessment:
    """Return a conservative review tier from bounded evidence inputs.

    The largest signal controls the tier so a serious concern cannot be
    averaged away by lower scores elsewhere.
    """
    _validate(evidence)
    signals = {
        "on-target uncertainty": evidence.on_target_uncertainty,
        "candidate off-target evidence": evidence.off_target_evidence,
        "network-impact evidence": evidence.network_impact_evidence,
        "welfare relevance": evidence.welfare_relevance,
    }
    score = max(signals.values())
    tier = (
        ReviewTier.HIGH_CONCERN_REVIEW
        if score >= 0.7
        else ReviewTier.ENHANCED_REVIEW
        if score >= 0.4
        else ReviewTier.STANDARD_REVIEW
    )
    rationale = tuple(label for label, value in signals.items() if value == score)
    return EditAssessment(concern_score=score, tier=tier, rationale=rationale)


def _validate(evidence: EditEvidence) -> None:
    for label, value in vars(evidence).items():
        if not 0 <= value <= 1:
            raise ValueError(f"{label} must be between 0 and 1.")
