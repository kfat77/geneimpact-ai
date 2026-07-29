"""Evidence labels that prevent prioritization scores from being misreported."""

from __future__ import annotations

from enum import Enum


class EvidenceLevel(str, Enum):
    """Claim level permitted by a study's validation record."""

    EXPLORATORY = "exploratory"
    REPLICATED = "replicated"
    CAUSAL_SUPPORT = "causal_support"


def permitted_wording(level: EvidenceLevel) -> str:
    """Return the strongest wording allowed for an evidence level."""
    return {
        EvidenceLevel.EXPLORATORY: "candidate association",
        EvidenceLevel.REPLICATED: "replicated association",
        EvidenceLevel.CAUSAL_SUPPORT: "causal evidence consistent with the stated assumptions",
    }[level]
