"""Transparent prioritization of candidate gene-pair interactions.

This module deliberately ranks research hypotheses only. It must not be used
to infer an individual's disease risk or to guide clinical decisions.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Mapping


@dataclass(frozen=True, order=True)
class InteractionResult:
    """A candidate gene-pair interaction and its traceable priority score."""

    priority: float
    genes: tuple[str, str]
    evidence_weight: float
    rationale: str


def rank_interactions(
    gene_scores: Mapping[str, float],
    evidence: Mapping[frozenset[str], float],
) -> list[InteractionResult]:
    """Rank pairs using bounded gene signals and independently curated evidence.

    Scores must be in ``[0, 1]``. The priority is the geometric mean of the
    two gene signals, multiplied by the evidence weight. This intentionally
    makes missing evidence produce a zero priority rather than inventing an
    interaction.
    """
    _validate_scores(gene_scores, "gene score")
    _validate_scores(evidence, "evidence weight")

    results: list[InteractionResult] = []
    for first, second in combinations(sorted(gene_scores), 2):
        weight = evidence.get(frozenset((first, second)), 0.0)
        priority = (gene_scores[first] * gene_scores[second]) ** 0.5 * weight
        results.append(
            InteractionResult(
                priority=priority,
                genes=(first, second),
                evidence_weight=weight,
                rationale="Geometric mean of bounded gene signals × curated evidence weight.",
            )
        )
    return sorted(results, reverse=True)


def _validate_scores(scores: Mapping[object, float], label: str) -> None:
    for key, value in scores.items():
        if not 0 <= value <= 1:
            raise ValueError(f"{label} for {key!r} must be between 0 and 1.")
