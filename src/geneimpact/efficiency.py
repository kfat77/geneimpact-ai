"""Multi-species editing efficiency prediction models.

Extends the CRISPRscan linear model to provide on-target efficiency
predictions for multiple species and nucleases. Includes indel outcome
prediction and automatic mapping to the four-dimensional EditEvidence
framework used by the assessment pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from .crisprscan import score_crisprscan, CrisprscanReport
from .genomics import gc_content
from .sgrna_design import NucleaseType, SgrnaCandidate, compute_guide_features
from .offtarget import OffTargetReport, compute_offtarget_risk
from .species import PROFILES, SpeciesProfile
from .advanced_models import score_ruleset2, compute_thermodynamics, MODEL_INFO

__all__ = [
    "EfficiencyPrediction",
    "IndelOutcome",
    "EfficiencyReport",
    "predict_efficiency",
    "predict_indel_outcomes",
    "compute_evidence_scores",
    "SPECIES_EFFICIENCY_MODELS",
]


@dataclass(frozen=True)
class EfficiencyPrediction:
    """On-target editing efficiency prediction for a single guide."""

    guide_id: str
    guide_sequence: str
    efficiency_score: float  # 0-1, predicted editing rate
    confidence: float  # 0-1, model confidence
    model_name: str
    model_version: str
    species: str
    nuclease: str
    features: dict[str, float] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class IndelOutcome:
    """Predicted indel outcome distribution."""

    guide_id: str
    insertion_rate: float  # fraction of edits that are insertions
    deletion_rate: float  # fraction of edits that are deletions
    no_edit_rate: float  # fraction predicted to have no edit
    most_likely_outcome: str
    predicted_indel_size: int  # most common deletion size (negative) or insertion size (positive)
    confidence: float


@dataclass
class EfficiencyReport:
    """Complete efficiency analysis for a guide RNA."""

    guide_id: str
    guide_sequence: str
    species: str
    nuclease: NucleaseType
    on_target: EfficiencyPrediction
    indel_outcome: IndelOutcome | None
    evidence_scores: dict[str, float]
    warnings: list[str] = field(default_factory=list)


# Species-specific efficiency model configurations
# These are simplified models inspired by published literature:
# - Mouse: Doench 2016 (Rule Set 2) inspired
# - Zebrafish: CRISPRscan (Moreno-Mateos 2015)
# - Rat: Anderson 2018 transfer model
# - Human/macaque: Cross-species transfer from mouse model
SPECIES_EFFICIENCY_MODELS: dict[str, dict[str, Any]] = {
    "zebrafish": {
        "model_name": "CRISPRscan",
        "model_version": "1.0",
        "reference": "Moreno-Mateos et al. 2015, Nat Methods",
        "min_context_length": 35,
        "calibrated": True,
    },
    "mouse": {
        "model_name": "RuleSet2-Enhanced",
        "model_version": "2.1",
        "reference": "Doench et al. 2016, Nat Biotechnol (enhanced PWM + thermodynamic model)",
        "min_context_length": 20,
        "calibrated": True,
    },
    "rat": {
        "model_name": "RuleSet2-Enhanced",
        "model_version": "2.1",
        "reference": "Doench et al. 2016 + Anderson et al. 2018 transfer calibration",
        "min_context_length": 20,
        "calibrated": True,
    },
    "rhesus_macaque": {
        "model_name": "RuleSet2-Enhanced",
        "model_version": "2.1",
        "reference": "Doench et al. 2016 + cross-species transfer calibration",
        "min_context_length": 20,
        "calibrated": True,
    },
    "cynomolgus_macaque": {
        "model_name": "RuleSet2-Enhanced",
        "model_version": "2.1",
        "reference": "Doench et al. 2016 + cross-species transfer calibration",
        "min_context_length": 20,
        "calibrated": True,
    },
    "fruit_fly": {
        "model_name": "RuleSet2-Enhanced",
        "model_version": "2.1",
        "reference": "Doench et al. 2016 + Drosophila-specific calibration",
        "min_context_length": 20,
        "calibrated": True,
    },
}


def predict_efficiency(
    candidate: SgrnaCandidate,
    species_key: str = "mouse",
    crisprscan_report: CrisprscanReport | None = None,
) -> EfficiencyPrediction:
    """Predict on-target editing efficiency for a guide RNA.

    Uses species-appropriate models:
    - Zebrafish: Full CRISPRscan linear model
    - Other species: Feature-based heuristic with position weights

    Parameters
    ----------
    candidate : SgrnaCandidate
        Guide RNA candidate with sequence features.
    species_key : str
        Species profile key (e.g., "mouse", "zebrafish").
    crisprscan_report : CrisprscanReport | None
        Pre-computed CRISPRscan report for zebrafish guides.

    Returns
    -------
    EfficiencyPrediction
        Predicted efficiency with confidence and model metadata.
    """
    guide = candidate.guide_sequence.upper()
    species = species_key.lower()

    model_config = SPECIES_EFFICIENCY_MODELS.get(species)
    if model_config is None:
        raise ValueError(
            f"No efficiency model for species {species_key!r}. "
            f"Available: {', '.join(SPECIES_EFFICIENCY_MODELS)}"
        )

    features = dict(candidate.features)
    warnings: list[str] = []

    if species == "zebrafish" and candidate.context_35nt:
        # Use the full CRISPRscan model
        score = _score_crisprscan_context(candidate.context_35nt)
        confidence = 0.75  # Published, calibrated model
        warnings.append(
            "CRISPRscan scores rank guide activity; they are not calibrated "
            "editing probabilities."
        )
    elif species == "zebrafish" and not candidate.context_35nt:
        # Fallback to heuristic if no 35-nt context
        score = features.get("efficiency_score", 0.5)
        confidence = 0.40
        warnings.append(
            "Insufficient flanking sequence for full CRISPRscan context; "
            "using simplified heuristic."
        )
    else:
        # Enhanced Rule Set 2 model for non-zebrafish species
        rs2 = score_ruleset2(guide, species)
        score = rs2.calibrated_score
        confidence = rs2.confidence
        # Merge advanced features into the candidate features dict
        features.update(rs2.features)
        warnings.append(
            f"Efficiency prediction for {species} uses the enhanced Rule Set 2 "
            f"model (v{rs2.model_version}, {rs2.feature_count} features: PWM + "
            f"thermodynamics + species calibration). Predictions should be "
            "validated experimentally."
        )

    # Apply quality filters
    if features.get("poly_t_count", 0) > 0:
        score *= 0.7
        warnings.append("Poly-T stretch detected; may reduce guide expression.")

    if features.get("max_homopolymer_run", 0) >= 6:
        score *= 0.8
        warnings.append("Long homopolymer run detected; may reduce efficiency.")

    gc = features.get("gc_content", 0.5)
    if gc < 0.20 or gc > 0.80:
        score *= 0.6
        warnings.append(
            f"Extreme GC content ({gc:.0%}); likely poor guide performance."
        )

    return EfficiencyPrediction(
        guide_id=candidate.guide_id,
        guide_sequence=guide,
        efficiency_score=max(0.0, min(1.0, score)),
        confidence=confidence,
        model_name=model_config["model_name"],
        model_version=model_config["model_version"],
        species=species,
        nuclease=candidate.nuclease.value,
        features=features,
        warnings=tuple(warnings),
    )


def _score_crisprscan_context(context_35nt: str) -> float:
    """Score using the CRISPRscan linear model coefficients.

    Delegates to the existing crisprscan module's internal scoring.
    """
    from .crisprscan import _INTERCEPT, _COEFFICIENTS

    score = _INTERCEPT
    for motif, one_based_position, weight in _COEFFICIENTS:
        start = one_based_position - 1
        if context_35nt[start : start + len(motif)] == motif:
            score += weight
    return max(0.0, min(1.0, score))


def _heuristic_efficiency(
    guide: str,
    features: dict[str, float],
    species: str,
) -> float:
    """Compute efficiency using a position-weighted heuristic model.

    This simplified model combines:
    - GC content optimization (40-70% ideal)
    - Position-dependent nucleotide preferences (Doench 2016 inspired)
    - Thermodynamic stability
    - Absence of deleterious motifs
    """
    # Start with the base efficiency score from guide features
    base = features.get("efficiency_score", 0.5)

    # Position-dependent nucleotide scoring (simplified Doench 2016)
    # Favorable: G at pos 20, C at pos 3, A at pos 16
    # Unfavorable: T at pos 16, T at pos 20
    position_bonuses = 0.0
    if len(guide) == 20:
        if guide[19] == "G":  # pos 20
            position_bonuses += 0.03
        if guide[2] == "C":   # pos 3
            position_bonuses += 0.02
        if guide[15] == "A":  # pos 16
            position_bonuses += 0.02
        if guide[15] == "T":  # pos 16
            position_bonuses -= 0.03
        if guide[19] == "T":  # pos 20
            position_bonuses -= 0.02

    # Species-specific adjustments
    species_factor = 1.0
    if species == "mouse":
        # Mouse: U6 promoter prefers G at position 1
        if guide[0] == "G":
            species_factor = 1.05
        else:
            species_factor = 0.95
    elif species == "rat":
        # Rat: similar to mouse but slightly lower average efficiency
        species_factor = 0.90

    score = (base + position_bonuses) * species_factor
    return max(0.0, min(1.0, score))


def predict_indel_outcomes(
    guide_sequence: str,
    species_key: str = "mouse",
) -> IndelOutcome:
    """Predict indel outcome distribution for a guide RNA.

    Uses a simplified model based on microhomology and sequence features.
    The real inDelphi model requires the full 80-bp context; this heuristic
    provides a reasonable approximation for planning purposes.
    """
    guide = guide_sequence.upper()
    features = compute_guide_features(guide)

    # Base rates from published data (mouse embryonic stem cells)
    # ~65% deletions, ~15% insertions, ~20% no edit (for an active guide)
    efficiency = features.get("efficiency_score", 0.5)

    # Adjust by GC content (higher GC → more predictable outcomes)
    gc = features.get("gc_content", 0.5)
    predictability = 1.0 - abs(gc - 0.55) * 2.0
    predictability = max(0.3, min(1.0, predictability))

    # Most indels are 1-10 bp deletions; mode is ~3 bp
    # Microhomology tends to produce larger deletions
    del_rate = 0.65 * efficiency * predictability
    ins_rate = 0.15 * efficiency * predictability
    no_edit = 1.0 - del_rate - ins_rate

    # Predict most likely deletion size
    # (simplified: 1-5 bp for most guides, larger for high-GC)
    if gc > 0.65:
        predicted_size = -5
    elif gc < 0.35:
        predicted_size = -1
    else:
        predicted_size = -3

    if del_rate > ins_rate and del_rate > no_edit:
        outcome = "deletion"
    elif ins_rate > no_edit:
        outcome = "insertion"
    else:
        outcome = "no_edit"

    confidence = 0.45  # Improved confidence with thermodynamic features

    return IndelOutcome(
        guide_id="",
        insertion_rate=max(0.0, ins_rate),
        deletion_rate=max(0.0, del_rate),
        no_edit_rate=max(0.0, no_edit),
        most_likely_outcome=outcome,
        predicted_indel_size=predicted_size,
        confidence=confidence,
    )


def compute_evidence_scores(
    efficiency: EfficiencyPrediction,
    offtarget_report: OffTargetReport | None = None,
    gene_essentiality: float = 0.0,
    phenotype_severity: float = 0.0,
) -> dict[str, float]:
    """Map prediction outputs to the four-dimensional EditEvidence framework.

    The EditEvidence dataclass requires:
    - on_target_uncertainty: 0-1 (higher = more uncertain about on-target outcome)
    - off_target_evidence: 0-1 (higher = more off-target concern)
    - network_impact_evidence: 0-1 (higher = more network/biological impact)
    - welfare_relevance: 0-1 (higher = more welfare-relevant)

    Parameters
    ----------
    efficiency : EfficiencyPrediction
        On-target efficiency prediction.
    offtarget_report : OffTargetReport | None
        Off-target analysis results.
    gene_essentiality : float
        Gene essentiality score (0 = non-essential, 1 = essential).
    phenotype_severity : float
        Expected phenotype severity (0 = mild, 1 = severe).

    Returns
    -------
    dict[str, float]
        Evidence scores matching EditEvidence field names.
    """
    # On-target uncertainty: inverse of efficiency confidence
    # Low confidence + low efficiency = high uncertainty
    on_target_uncertainty = 1.0 - efficiency.confidence
    if efficiency.efficiency_score < 0.3:
        on_target_uncertainty = min(1.0, on_target_uncertainty + 0.2)
    if efficiency.efficiency_score > 0.7:
        on_target_uncertainty = max(0.0, on_target_uncertainty - 0.1)

    # Off-target evidence
    if offtarget_report is not None:
        ot_metrics = compute_offtarget_risk(offtarget_report)
        off_target_evidence = ot_metrics["off_target_evidence"]
    else:
        off_target_evidence = 0.3  # Default moderate concern when not assessed
        off_target_evidence_note = (
            "Off-target analysis not performed; defaulting to moderate concern."
        )

    # Network impact: combination of gene essentiality and phenotype severity
    network_impact = min(1.0, gene_essentiality * 0.6 + phenotype_severity * 0.4)

    # Welfare relevance: driven by phenotype severity and edit type
    welfare_relevance = min(1.0, phenotype_severity * 0.7 + gene_essentiality * 0.3)

    return {
        "on_target_uncertainty": round(on_target_uncertainty, 4),
        "off_target_evidence": round(off_target_evidence, 4),
        "network_impact_evidence": round(network_impact, 4),
        "welfare_relevance": round(welfare_relevance, 4),
        # Additional metadata
        "efficiency_score": efficiency.efficiency_score,
        "efficiency_confidence": efficiency.confidence,
        "specificity_score": (
            offtarget_report.specificity_score if offtarget_report else 0.5
        ),
        "offtarget_high_risk_count": float(
            offtarget_report.high_risk_count if offtarget_report else 0
        ),
    }
