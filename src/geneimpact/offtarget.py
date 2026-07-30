"""Off-target detection and scoring for CRISPR guide RNAs.

Implements a seed-region-weighted mismatch scoring algorithm to identify
and rank potential off-target sites in a reference genome or user-provided
sequence. Supports both Cas9 (NGG) and Cas12a (TTTV) PAM systems.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Sequence

from .genomics import reverse_complement
from .sgrna_design import NucleaseType, PAM_PATTERNS, PamPattern

__all__ = [
    "OffTargetSite",
    "OffTargetReport",
    "MismatchPattern",
    "find_offtargets",
    "score_offtarget",
    "compute_offtarget_risk",
]


# Position-dependent mismatch penalties for SpCas9 (Doench 2016 inspired)
# Higher penalty = more critical position (PAM-proximal seed region)
# Index 0 = PAM-distal, Index 19 = PAM-proximal
_SPcas9_POSITION_WEIGHTS: tuple[float, ...] = (
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5,
)


@dataclass(frozen=True)
class MismatchPattern:
    """Pattern of mismatches between guide and off-target."""

    positions: tuple[int, ...]  # 0-based mismatch positions
    count: int

    @property
    def has_seed_mismatch(self) -> bool:
        """True if any mismatch is in the seed region (PAM-proximal 8-12 nt)."""
        return any(pos >= 10 for pos in self.positions)

    @property
    def seed_mismatch_count(self) -> int:
        """Number of mismatches in the seed region."""
        return sum(1 for pos in self.positions if pos >= 10)


@dataclass(frozen=True)
class OffTargetSite:
    """A predicted off-target site."""

    chrom: str
    start: int  # 1-based start of the guide-binding region
    end: int  # 1-based end
    strand: str
    off_target_sequence: str  # 20-nt genomic sequence at this site
    pam: str
    mismatch_count: int
    mismatch_positions: tuple[int, ...]
    mismatch_pattern: MismatchPattern
    score: float  # 0-1, higher = more likely to be cut
    risk_level: str  # "high", "moderate", "low"


@dataclass
class OffTargetReport:
    """Complete off-target analysis for a guide RNA."""

    guide_sequence: str
    nuclease: NucleaseType
    total_sites_scanned: int
    off_targets: list[OffTargetSite]
    max_mismatches_searched: int
    warnings: list[str] = field(default_factory=list)

    @property
    def high_risk_count(self) -> int:
        return sum(1 for ot in self.off_targets if ot.risk_level == "high")

    @property
    def moderate_risk_count(self) -> int:
        return sum(1 for ot in self.off_targets if ot.risk_level == "moderate")

    @property
    def low_risk_count(self) -> int:
        return sum(1 for ot in self.off_targets if ot.risk_level == "low")

    @property
    def specificity_score(self) -> float:
        """Overall guide specificity (1 = perfectly specific, 0 = many off-targets)."""
        if not self.off_targets:
            return 1.0
        # Weighted by risk and mismatch count
        total_penalty = 0.0
        for ot in self.off_targets:
            weight = {"high": 1.0, "moderate": 0.5, "low": 0.2}[ot.risk_level]
            total_penalty += weight * (1.0 / max(1, ot.mismatch_count))
        return max(0.0, 1.0 - total_penalty / 10.0)


def find_offtargets(
    guide_sequence: str,
    reference_sequences: dict[str, str] | None = None,
    fasta_reader=None,
    nuclease: NucleaseType = NucleaseType.SPCAS9,
    max_mismatches: int = 4,
    pam_pattern: str | None = None,
    guide_length: int = 20,
) -> OffTargetReport:
    """Find off-target sites for a guide RNA in reference sequences.

    Parameters
    ----------
    guide_sequence : str
        20-nt guide RNA sequence (without PAM).
    reference_sequences : dict[str, str] | None
        Dictionary of {seq_id: sequence} to search. If None and fasta_reader
        is provided, all sequences from the reader will be searched.
    fasta_reader : FastaReader | None
        A FastaReader instance to load sequences from. Used when
        reference_sequences is None.
    nuclease : NucleaseType
        CRISPR nuclease type for PAM detection.
    max_mismatches : int
        Maximum number of mismatches to allow (1-4 recommended).
    pam_pattern : str | None
        Custom PAM pattern to override the nuclease default.

    Returns
    -------
    OffTargetReport
        Off-target analysis results.
    """
    guide = guide_sequence.upper()
    if len(guide) != guide_length:
        raise ValueError(f"guide must be {guide_length} nt, got {len(guide)}")

    if max_mismatches < 0 or max_mismatches > 6:
        raise ValueError(f"max_mismatches must be 0-6, got {max_mismatches}")

    # Determine PAM pattern
    if pam_pattern is not None:
        pam = PamPattern(
            pattern=pam_pattern, nuclease=nuclease,
            position=PAM_PATTERNS[nuclease].position,
            length=len(pam_pattern),
        )
    else:
        pam = PAM_PATTERNS[nuclease]

    warnings: list[str] = []

    # Gather sequences to search
    if reference_sequences is None:
        if fasta_reader is None:
            raise ValueError(
                "Either reference_sequences or fasta_reader must be provided"
            )
        reference_sequences = {}
        for seq_id in fasta_reader.sequence_ids:
            record = fasta_reader[seq_id]
            reference_sequences[seq_id] = record.sequence

    total_scanned = 0
    off_targets: list[OffTargetSite] = []

    for chrom, ref_seq in reference_sequences.items():
        ref = ref_seq.upper()
        # Scan plus strand
        for ot in _scan_strand(
            guide, ref, chrom, "+", pam, max_mismatches, guide_length, 1
        ):
            off_targets.append(ot)
            total_scanned += 1

        # Scan minus strand (check reverse complement)
        rc_ref = reverse_complement(ref)
        for ot in _scan_strand(
            guide, rc_ref, chrom, "-", pam, max_mismatches, guide_length, 1
        ):
            # Adjust coordinates back to plus strand
            adjusted_start = len(ref) - ot.end + 1
            adjusted_end = len(ref) - ot.start + 1
            off_targets.append(
                OffTargetSite(
                    chrom=ot.chrom,
                    start=adjusted_start,
                    end=adjusted_end,
                    strand=ot.strand,
                    off_target_sequence=ot.off_target_sequence,
                    pam=ot.pam,
                    mismatch_count=ot.mismatch_count,
                    mismatch_positions=ot.mismatch_positions,
                    mismatch_pattern=ot.mismatch_pattern,
                    score=ot.score,
                    risk_level=ot.risk_level,
                )
            )
            total_scanned += 1

    # Sort by score (descending) then by mismatch count (ascending)
    off_targets.sort(key=lambda ot: (-ot.score, ot.mismatch_count))

    # Limit results
    if len(off_targets) > 1000:
        warnings.append(
            f"Found {len(off_targets)} off-target sites; "
            "returning top 1000 by score."
        )
        off_targets = off_targets[:1000]

    return OffTargetReport(
        guide_sequence=guide,
        nuclease=nuclease,
        total_sites_scanned=total_scanned,
        off_targets=off_targets,
        max_mismatches_searched=max_mismatches,
        warnings=warnings,
    )


def _scan_strand(
    guide: str,
    ref: str,
    chrom: str,
    strand: str,
    pam: PamPattern,
    max_mismatches: int,
    guide_length: int,
    base_pos: int,
) -> Iterator[OffTargetSite]:
    """Scan one strand of a reference sequence for off-target sites."""
    total_window = guide_length + pam.length

    for i in range(len(ref) - total_window + 1):
        if pam.position == "3prime":
            target = ref[i : i + guide_length]
            pam_seq = ref[i + guide_length : i + total_window]
            guide_start = base_pos + i
            guide_end = guide_start + guide_length - 1
        elif pam.position == "5prime":
            pam_seq = ref[i : i + pam.length]
            target = ref[i + pam.length : i + total_window]
            guide_start = base_pos + i + pam.length
            guide_end = guide_start + guide_length - 1
        else:
            continue

        # Skip if N in target or PAM
        if "N" in target or "N" in pam_seq:
            continue

        # Check PAM match
        if not pam.matches(pam_seq):
            continue

        # Count mismatches
        mismatches: list[int] = []
        for j in range(guide_length):
            if guide[j] != target[j]:
                mismatches.append(j)

        mm_count = len(mismatches)
        if mm_count > max_mismatches:
            continue

        # Skip exact match (that's the on-target)
        if mm_count == 0:
            continue

        mm_pattern = MismatchPattern(
            positions=tuple(mismatches),
            count=mm_count,
        )

        score = score_offtarget(guide, target, mm_pattern, pam.nuclease)
        risk = _risk_level(mm_count, score, mm_pattern)

        yield OffTargetSite(
            chrom=chrom,
            start=guide_start,
            end=guide_end,
            strand=strand,
            off_target_sequence=target,
            pam=pam_seq,
            mismatch_count=mm_count,
            mismatch_positions=tuple(mismatches),
            mismatch_pattern=mm_pattern,
            score=score,
            risk_level=risk,
        )


def score_offtarget(
    guide: str,
    target: str,
    pattern: MismatchPattern,
    nuclease: NucleaseType,
) -> float:
    """Score an off-target site for cleavage likelihood.

    Returns a value in [0, 1] where 1 = as likely as on-target.
    Uses position-dependent mismatch penalties inspired by
    Hsu et al. (2013) and Doench et al. (2016).
    """
    # Perfect match (0 mismatches) = maximum score
    if pattern.count == 0:
        return 1.0

    if nuclease == NucleaseType.SPCAS9:
        weights = _SPcas9_POSITION_WEIGHTS
    else:
        # For other nucleases, use a uniform penalty
        weights = tuple(0.5 for _ in range(len(guide)))

    # Base penalty per mismatch (exponential decay with mismatch count)
    # 1 mismatch: ~0.7, 2: ~0.3, 3: ~0.1, 4: ~0.03
    base_penalty = {
        1: 0.30,
        2: 0.70,
        3: 0.90,
        4: 0.97,
        5: 0.99,
        6: 0.999,
    }
    base_score = 1.0 - base_penalty.get(pattern.count, 0.999)

    # Apply position-dependent weights
    position_penalty = 0.0
    for pos in pattern.positions:
        position_penalty += weights[pos] * 0.1

    # Seed mismatches are more detrimental
    if pattern.has_seed_mismatch:
        position_penalty += pattern.seed_mismatch_count * 0.15

    score = base_score * (1.0 - position_penalty)
    return max(0.0, min(1.0, score))


def _risk_level(
    mismatch_count: int,
    score: float,
    pattern: MismatchPattern,
) -> str:
    """Classify off-target risk level based on mismatch count and score."""
    if mismatch_count <= 1 and score > 0.5:
        return "high"
    if mismatch_count <= 2 and score > 0.2:
        return "moderate"
    if mismatch_count <= 3 and score > 0.05:
        return "moderate"
    return "low"


def compute_offtarget_risk(report: OffTargetReport) -> dict[str, float]:
    """Compute summary risk metrics from an off-target report.

    Returns a dictionary suitable for mapping to the EditEvidence
    off_target_evidence field (0-1 scale, higher = more concern).
    """
    if not report.off_targets:
        return {
            "off_target_evidence": 0.0,
            "high_risk_count": 0.0,
            "total_off_targets": 0.0,
            "specificity_score": 1.0,
            "max_offtarget_score": 0.0,
        }

    high_count = report.high_risk_count
    mod_count = report.moderate_risk_count
    total = len(report.off_targets)
    max_score = max(ot.score for ot in report.off_targets)
    specificity = report.specificity_score

    # Compute concern score (0-1, higher = more concerning)
    # Weighted by risk level
    concern = (
        high_count * 0.3
        + mod_count * 0.1
        + total * 0.01
    )
    concern += max_score * 0.3  # Add contribution from the worst site
    concern += (1.0 - specificity) * 0.3  # Add contribution from low specificity

    return {
        "off_target_evidence": min(1.0, concern),
        "high_risk_count": float(high_count),
        "total_off_targets": float(total),
        "specificity_score": specificity,
        "max_offtarget_score": max_score,
    }
