"""End-to-end gene editing prediction pipeline.

Connects: sequence input -> sgRNA design -> efficiency prediction ->
off-target detection -> evidence scoring -> assessment -> report generation.
This module automates the full workflow that previously required manual
evidence score input.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping

from .edit_assessment import EditAssessment, EditEvidence, ReviewTier, assess_edit
from .efficiency import (
    EfficiencyPrediction,
    EfficiencyReport,
    compute_evidence_scores,
    predict_efficiency,
    predict_indel_outcomes,
)
from .genomics import FastaReader, gc_content, reverse_complement, validate_dna_sequence
from .offtarget import OffTargetReport, compute_offtarget_risk, find_offtargets
from .provenance import StudyContext, create_record
from .sgrna_design import NucleaseType, SgrnaCandidate, SgrnaDesignResult, design_sgrnas
from .species import PROFILES, SpeciesProfile, validate_study_context

__all__ = [
    "PipelineConfig",
    "GuideResult",
    "PipelineReport",
    "run_pipeline",
    "run_pipeline_from_fasta",
]


@dataclass
class PipelineConfig:
    """Configuration for the prediction pipeline."""

    species: str = "mouse"
    nuclease: NucleaseType = NucleaseType.SPCAS9
    guide_length: int = 20
    max_candidates: int = 50
    max_offtargets: int = 4
    search_reference: bool = True
    gene_essentiality: float = 0.0
    phenotype_severity: float = 0.0
    min_efficiency: float = 0.3
    min_specificity: float = 0.5
    top_k: int = 10


@dataclass
class GuideResult:
    """Complete analysis for a single guide RNA."""

    candidate: SgrnaCandidate
    efficiency: EfficiencyPrediction
    indel_outcome: Any | None = None
    offtarget_report: OffTargetReport | None = None
    evidence_scores: dict[str, float] = field(default_factory=dict)
    assessment: EditAssessment | None = None
    rank: int = 0
    recommendation: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dictionary."""
        result: dict[str, Any] = {
            "rank": self.rank,
            "recommendation": self.recommendation,
            "guide_id": self.candidate.guide_id,
            "guide_sequence": self.candidate.guide_sequence,
            "pam": self.candidate.pam,
            "strand": self.candidate.strand,
            "chrom": self.candidate.chrom,
            "start": self.candidate.start,
            "end": self.candidate.end,
            "gc_content": round(self.candidate.gc_content, 4),
            "nuclease": self.candidate.nuclease.value,
            "efficiency": {
                "score": round(self.efficiency.efficiency_score, 4),
                "confidence": round(self.efficiency.confidence, 4),
                "model": self.efficiency.model_name,
                "model_version": self.efficiency.model_version,
                "warnings": list(self.efficiency.warnings),
            },
            "evidence_scores": {
                k: round(v, 4) for k, v in self.evidence_scores.items()
            },
        }
        if self.indel_outcome is not None:
            result["indel_outcome"] = {
                "insertion_rate": round(self.indel_outcome.insertion_rate, 4),
                "deletion_rate": round(self.indel_outcome.deletion_rate, 4),
                "no_edit_rate": round(self.indel_outcome.no_edit_rate, 4),
                "most_likely": self.indel_outcome.most_likely_outcome,
                "predicted_size": self.indel_outcome.predicted_indel_size,
            }
        if self.offtarget_report is not None:
            result["offtarget"] = {
                "total_sites": self.offtarget_report.total_sites_scanned,
                "high_risk": self.offtarget_report.high_risk_count,
                "moderate_risk": self.offtarget_report.moderate_risk_count,
                "low_risk": self.offtarget_report.low_risk_count,
                "specificity_score": round(
                    self.offtarget_report.specificity_score, 4
                ),
                "top_sites": [
                    {
                        "chrom": ot.chrom,
                        "start": ot.start,
                        "end": ot.end,
                        "strand": ot.strand,
                        "mismatches": ot.mismatch_count,
                        "score": round(ot.score, 4),
                        "risk": ot.risk_level,
                        "sequence": ot.off_target_sequence,
                        "pam": ot.pam,
                    }
                    for ot in self.offtarget_report.off_targets[:10]
                ],
            }
        if self.assessment is not None:
            result["assessment"] = {
                "concern_score": round(self.assessment.concern_score, 4),
                "tier": self.assessment.tier.value,
                "rationale": list(self.assessment.rationale),
            }
        return result


@dataclass
class PipelineReport:
    """Complete pipeline report for all guides."""

    config: PipelineConfig
    study_context: StudyContext
    species_validation: Any
    guides: list[GuideResult]
    warnings: list[str] = field(default_factory=list)
    pipeline_version: str = "1.0.0"
    timestamp: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dictionary."""
        return {
            "pipeline_version": self.pipeline_version,
            "timestamp": self.timestamp,
            "config": {
                "species": self.config.species,
                "nuclease": self.config.nuclease.value,
                "guide_length": self.config.guide_length,
                "max_candidates": self.config.max_candidates,
                "max_offtargets": self.config.max_offtargets,
                "gene_essentiality": self.config.gene_essentiality,
                "phenotype_severity": self.config.phenotype_severity,
                "min_efficiency": self.config.min_efficiency,
                "min_specificity": self.config.min_specificity,
                "top_k": self.config.top_k,
            },
            "study_context": asdict(self.study_context),
            "species_validation": {
                "supported": self.species_validation.supported,
                "profile_key": self.species_validation.profile_key,
                "errors": list(self.species_validation.errors),
                "warnings": list(self.species_validation.warnings),
            },
            "guides": [g.to_dict() for g in self.guides],
            "warnings": self.warnings,
            "report_notice": (
                "Research decision-support only. Predictions are based on "
                "computational models and must be validated experimentally. "
                "This report does not establish safety, authorize an edit, "
                "or replace ethics, biosafety, veterinary, or experimental review."
            ),
        }

    def to_json(self, path: str | Path | None = None) -> str:
        """Serialize to JSON string or write to file."""
        text = json.dumps(self.to_dict(), indent=2, ensure_ascii=False) + "\n"
        if path is not None:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_text(text, encoding="utf-8")
        return text


def run_pipeline(
    sequence: str,
    config: PipelineConfig,
    study_context: StudyContext,
    reference_sequences: dict[str, str] | None = None,
    fasta_reader: FastaReader | None = None,
) -> PipelineReport:
    """Run the complete gene editing prediction pipeline.

    Parameters
    ----------
    sequence : str
        Target DNA sequence to design guides for.
    config : PipelineConfig
        Pipeline configuration.
    study_context : StudyContext
        Study context for species validation.
    reference_sequences : dict[str, str] | None
        Reference genome sequences for off-target search.
    fasta_reader : FastaReader | None
        FASTA reader for off-target search (alternative to reference_sequences).

    Returns
    -------
    PipelineReport
        Complete pipeline results.
    """
    # Validate species
    species_validation = validate_study_context(study_context)
    if species_validation.errors:
        raise ValueError(
            "Invalid study context: " + " ".join(species_validation.errors)
        )

    warnings: list[str] = []
    if species_validation.warnings:
        warnings.extend(species_validation.warnings)

    # Step 1: Design sgRNAs
    design_result = design_sgrnas(
        sequence=sequence,
        chrom=study_context.species,
        nuclease=config.nuclease,
        guide_length=config.guide_length,
        max_candidates=config.max_candidates,
    )
    warnings.extend(design_result.warnings)

    if not design_result.candidates:
        return PipelineReport(
            config=config,
            study_context=study_context,
            species_validation=species_validation,
            guides=[],
            warnings=warnings + ["No sgRNA candidates found in the target sequence."],
            timestamp=datetime.now(timezone.utc).isoformat(),
        )

    # Step 2: For each candidate, predict efficiency + off-targets
    guide_results: list[GuideResult] = []

    for candidate in design_result.candidates:
        # Predict on-target efficiency
        efficiency = predict_efficiency(
            candidate=candidate,
            species_key=study_context.species,
        )

        # Predict indel outcomes
        indel = predict_indel_outcomes(
            guide_sequence=candidate.guide_sequence,
            species_key=study_context.species,
        )

        # Off-target search
        ot_report: OffTargetReport | None = None
        if config.search_reference:
            try:
                ot_report = find_offtargets(
                    guide_sequence=candidate.guide_sequence,
                    reference_sequences=reference_sequences,
                    fasta_reader=fasta_reader,
                    nuclease=config.nuclease,
                    max_mismatches=config.max_offtargets,
                )
            except (ValueError, KeyError) as e:
                warnings.append(
                    f"Off-target search failed for {candidate.guide_id}: {e}"
                )
                ot_report = None

        # Compute evidence scores
        evidence_scores = compute_evidence_scores(
            efficiency=efficiency,
            offtarget_report=ot_report,
            gene_essentiality=config.gene_essentiality,
            phenotype_severity=config.phenotype_severity,
        )

        # Run assessment
        evidence = EditEvidence(
            on_target_uncertainty=evidence_scores["on_target_uncertainty"],
            off_target_evidence=evidence_scores["off_target_evidence"],
            network_impact_evidence=evidence_scores["network_impact_evidence"],
            welfare_relevance=evidence_scores["welfare_relevance"],
        )
        assessment = assess_edit(evidence)

        # Generate recommendation
        recommendation = _generate_recommendation(
            efficiency.efficiency_score,
            ot_report.specificity_score if ot_report else 0.5,
            assessment.tier,
            config,
        )

        guide_results.append(GuideResult(
            candidate=candidate,
            efficiency=efficiency,
            indel_outcome=indel,
            offtarget_report=ot_report,
            evidence_scores=evidence_scores,
            assessment=assessment,
            recommendation=recommendation,
        ))

    # Step 3: Rank guides
    guide_results.sort(
        key=lambda g: (
            -g.efficiency.efficiency_score,
            -(g.offtarget_report.specificity_score if g.offtarget_report else 0.5),
            g.assessment.concern_score,
        )
    )
    for i, g in enumerate(guide_results):
        g.rank = i + 1

    # Step 4: Filter to top_k
    if len(guide_results) > config.top_k:
        guide_results = guide_results[: config.top_k]

    return PipelineReport(
        config=config,
        study_context=study_context,
        species_validation=species_validation,
        guides=guide_results,
        warnings=warnings,
        timestamp=datetime.now(timezone.utc).isoformat(),
    )


def run_pipeline_from_fasta(
    fasta_path: str | Path,
    target_chrom: str,
    target_start: int,
    target_end: int,
    config: PipelineConfig,
    study_context: StudyContext,
    reference_fasta: str | Path | None = None,
) -> PipelineReport:
    """Run the pipeline using a FASTA file as input.

    Parameters
    ----------
    fasta_path : str | Path
        Path to the target sequence FASTA file.
    target_chrom : str
        Chromosome/sequence ID to extract the target from.
    target_start : int
        1-based start of the target region.
    target_end : int
        1-based end of the target region.
    config : PipelineConfig
        Pipeline configuration.
    study_context : StudyContext
        Study context.
    reference_fasta : str | Path | None
        Path to the reference genome FASTA for off-target search.
        If None, off-target search is skipped.
    """
    reader = FastaReader(fasta_path)
    target_seq = reader.fetch(target_chrom, target_start, target_end)

    ref_reader = None
    ref_seqs = None
    if reference_fasta is not None:
        ref_reader = FastaReader(reference_fasta)
    elif config.search_reference:
        # Use the target FASTA as both target and reference
        ref_seqs = {reader.sequence_ids[i]: reader[reader.sequence_ids[i]].sequence
                     for i in range(min(5, len(reader.sequence_ids)))}

    return run_pipeline(
        sequence=target_seq,
        config=config,
        study_context=study_context,
        reference_sequences=ref_seqs,
        fasta_reader=ref_reader,
    )


def _generate_recommendation(
    efficiency: float,
    specificity: float,
    tier: ReviewTier,
    config: PipelineConfig,
) -> str:
    """Generate a human-readable recommendation for a guide."""
    parts: list[str] = []

    if efficiency >= config.min_efficiency and specificity >= config.min_specificity:
        parts.append("Recommended candidate")
    elif efficiency >= config.min_efficiency:
        parts.append("Good efficiency but check off-target profile")
    elif specificity >= config.min_specificity:
        parts.append("Good specificity but efficiency may be suboptimal")
    else:
        parts.append("Use with caution")

    if tier == ReviewTier.HIGH_CONCERN_REVIEW:
        parts.append("requires enhanced review due to high concern signals")
    elif tier == ReviewTier.ENHANCED_REVIEW:
        parts.append("requires enhanced review")

    return "; ".join(parts) + "."
