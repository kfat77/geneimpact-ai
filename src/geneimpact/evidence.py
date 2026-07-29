"""Evidence labels that prevent animal-edit concerns from being misreported."""

from __future__ import annotations

from enum import Enum


class EvidenceLevel(str, Enum):
    """Concern level permitted by a study's validation record."""

    EXPLORATORY = "exploratory"
    REPLICATED = "replicated"
    CAUSAL_SUPPORT = "causal_support"


def permitted_wording(level: EvidenceLevel) -> str:
    """Return the strongest wording allowed for an evidence level."""
    return {
        EvidenceLevel.EXPLORATORY: "candidate concern",
        EvidenceLevel.REPLICATED: "replicated concern",
        EvidenceLevel.CAUSAL_SUPPORT: "strongly supported concern",
    }[level]
