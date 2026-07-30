"""sgRNA design and PAM search engine for CRISPR systems.

Scans genomic sequences for guide RNA candidates with canonical PAM sites,
computes sequence features (GC content, position, thermodynamic properties),
and ranks candidates by predicted on-target efficiency.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterator, Sequence

from .genomics import gc_content, reverse_complement, validate_dna_sequence

__all__ = [
    "NucleaseType",
    "PamPattern",
    "SgrnaCandidate",
    "SgrnaDesignResult",
    "PAM_PATTERNS",
    "design_sgrnas",
    "search_pam_sites",
    "compute_guide_features",
]


class NucleaseType(str, Enum):
    """Supported CRISPR nucleases."""

    SPCAS9 = "SpCas9"
    SACAS9 = "SaCas9"
    CAS12A = "Cas12a"
    CASX = "CasX"


@dataclass(frozen=True)
class PamPattern:
    """A PAM recognition pattern with IUPAC ambiguity codes."""

    pattern: str
    nuclease: NucleaseType
    position: str  # "3prime" or "5prime" relative to guide
    length: int

    def matches(self, candidate: str) -> bool:
        """Check if a DNA string matches this PAM pattern."""
        iupac: dict[str, frozenset[str]] = {
            "A": frozenset("A"),
            "C": frozenset("C"),
            "G": frozenset("G"),
            "T": frozenset("T"),
            "R": frozenset("AG"),
            "Y": frozenset("CT"),
            "S": frozenset("GC"),
            "W": frozenset("AT"),
            "K": frozenset("GT"),
            "M": frozenset("AC"),
            "B": frozenset("CGT"),
            "D": frozenset("AGT"),
            "H": frozenset("ACT"),
            "V": frozenset("ACG"),
            "N": frozenset("ACGT"),
        }
        if len(candidate) != self.length:
            return False
        for pat_base, cand_base in zip(self.pattern.upper(), candidate.upper()):
            if cand_base not in iupac.get(pat_base, frozenset()):
                return False
        return True


# Registered PAM patterns for supported nucleases
PAM_PATTERNS: dict[NucleaseType, PamPattern] = {
    NucleaseType.SPCAS9: PamPattern(
        pattern="NGG", nuclease=NucleaseType.SPCAS9,
        position="3prime", length=3,
    ),
    NucleaseType.SACAS9: PamPattern(
        pattern="NNGRRT", nuclease=NucleaseType.SACAS9,
        position="3prime", length=6,
    ),
    NucleaseType.CAS12A: PamPattern(
        pattern="TTTV", nuclease=NucleaseType.CAS12A,
        position="5prime", length=4,
    ),
}


@dataclass(frozen=True)
class SgrnaCandidate:
    """A designed sgRNA candidate with sequence features."""

    guide_id: str
    guide_sequence: str  # 20-nt guide (without PAM)
    pam: str
    pam_strand: str  # "+" or "-"
    chrom: str
    start: int  # 1-based guide start on + strand
    end: int  # 1-based guide end on + strand
    strand: str  # guide is on + or - strand
    context_30nt: str  # 30-nt context for scoring
    context_35nt: str  # 35-nt context (CRISPRscan format) if available
    gc_content: float
    nuclease: NucleaseType
    features: dict[str, float] = field(default_factory=dict)

    @property
    def full_target(self) -> str:
        """Guide + PAM on the targeting strand."""
        return self.guide_sequence + self.pam if self.pam_strand == "+" else self.pam + self.guide_sequence


@dataclass
class SgrnaDesignResult:
    """Result of sgRNA design from a target sequence."""

    target_id: str
    nuclease: NucleaseType
    candidates: list[SgrnaCandidate]
    warnings: list[str] = field(default_factory=list)

    @property
    def count(self) -> int:
        return len(self.candidates)

    def top_candidates(self, n: int = 10) -> list[SgrnaCandidate]:
        """Return top N candidates sorted by efficiency score if available."""
        scored = sorted(
            self.candidates,
            key=lambda c: c.features.get("efficiency_score", 0.0),
            reverse=True,
        )
        return scored[:n]


def design_sgrnas(
    sequence: str,
    chrom: str = "target",
    start_position: int = 1,
    nuclease: NucleaseType = NucleaseType.SPCAS9,
    guide_length: int = 20,
    max_candidates: int = 500,
) -> SgrnaDesignResult:
    """Design sgRNA candidates from a target DNA sequence.

    Scans both strands for PAM sites and generates guide RNA candidates
    with computed sequence features.

    Parameters
    ----------
    sequence : str
        Target DNA sequence (5' to 3').
    chrom : str
        Chromosome or sequence identifier for coordinate reporting.
    start_position : int
        1-based start coordinate of the input sequence in the genome.
    nuclease : NucleaseType
        CRISPR nuclease to design guides for.
    guide_length : int
        Length of the guide RNA (default 20 for SpCas9).
    max_candidates : int
        Maximum number of candidates to return.

    Returns
    -------
    SgrnaDesignResult
        Design results with candidates and any warnings.
    """
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    validate_dna_sequence(seq, allowed="ACGTN")

    warnings: list[str] = []
    if "N" in seq:
        n_count = seq.count("N")
        warnings.append(
            f"Sequence contains {n_count} ambiguous bases (N); "
            "these positions will be skipped during PAM search."
        )

    if nuclease not in PAM_PATTERNS:
        raise ValueError(
            f"No PAM pattern registered for nuclease {nuclease.value}"
        )

    pam = PAM_PATTERNS[nuclease]
    candidates: list[SgrnaCandidate] = []

    # Search + strand
    for candidate in _scan_plus_strand(
        seq, chrom, start_position, pam, nuclease, guide_length
    ):
        candidates.append(candidate)
        if len(candidates) >= max_candidates:
            break

    # Search - strand
    if len(candidates) < max_candidates:
        for candidate in _scan_minus_strand(
            seq, chrom, start_position, pam, nuclease, guide_length
        ):
            candidates.append(candidate)
            if len(candidates) >= max_candidates:
                break

    # Compute features for each candidate
    for i, candidate in enumerate(candidates):
        features = compute_guide_features(candidate.guide_sequence)
        candidate.features.update(features)

    if not candidates:
        warnings.append(
            f"No {pam.pattern} PAM sites found for {nuclease.value} "
            f"in the provided sequence ({len(seq)} bp)."
        )

    return SgrnaDesignResult(
        target_id=chrom,
        nuclease=nuclease,
        candidates=candidates,
        warnings=warnings,
    )


def _scan_plus_strand(
    seq: str,
    chrom: str,
    base_pos: int,
    pam: PamPattern,
    nuclease: NucleaseType,
    guide_length: int,
) -> Iterator[SgrnaCandidate]:
    """Scan the plus strand for PAM sites and yield candidates."""
    for i in range(len(seq) - guide_length - pam.length + 1):
        # Check for N in the guide or PAM region
        region = seq[i : i + guide_length + pam.length]
        if "N" in region:
            continue

        if pam.position == "3prime":
            # Guide is upstream, PAM is downstream
            guide = seq[i : i + guide_length]
            pam_seq = seq[i + guide_length : i + guide_length + pam.length]
            if not pam.matches(pam_seq):
                continue
            guide_start = base_pos + i
            guide_end = guide_start + guide_length - 1

        elif pam.position == "5prime":
            # PAM is upstream, guide is downstream
            pam_seq = seq[i : i + pam.length]
            if not pam.matches(pam_seq):
                continue
            guide = seq[i + pam.length : i + pam.length + guide_length]
            guide_start = base_pos + i + pam.length
            guide_end = guide_start + guide_length - 1
        else:
            continue

        # Extract context windows
        ctx_start = max(0, i - 4)
        ctx_end = min(len(seq), i + guide_length + pam.length + 3)
        context_30 = seq[ctx_start:ctx_end]

        # 35-nt context for CRISPRscan (if available)
        cs_start = max(0, i - 6)
        cs_end = min(len(seq), i + guide_length + pam.length + 6)
        context_35 = seq[cs_start:cs_end]
        if len(context_35) < 35:
            context_35 = ""  # Not enough flanking sequence

        gc = gc_content(guide)
        guide_id = f"{chrom}:{guide_start}-{guide_end}(+)"

        yield SgrnaCandidate(
            guide_id=guide_id,
            guide_sequence=guide,
            pam=pam_seq,
            pam_strand="+",
            chrom=chrom,
            start=guide_start,
            end=guide_end,
            strand="+",
            context_30nt=context_30,
            context_35nt=context_35,
            gc_content=gc,
            nuclease=nuclease,
        )


def _scan_minus_strand(
    seq: str,
    chrom: str,
    base_pos: int,
    pam: PamPattern,
    nuclease: NucleaseType,
    guide_length: int,
) -> Iterator[SgrnaCandidate]:
    """Scan the minus strand by checking the reverse complement."""
    rc = reverse_complement(seq)
    for i in range(len(rc) - guide_length - pam.length + 1):
        region = rc[i : i + guide_length + pam.length]
        if "N" in region:
            continue

        if pam.position == "3prime":
            guide = rc[i : i + guide_length]
            pam_seq = rc[i + guide_length : i + guide_length + pam.length]
            if not pam.matches(pam_seq):
                continue
            # Map back to + strand coordinates
            guide_end_plus = base_pos + len(seq) - i - 1
            guide_start_plus = guide_end_plus - guide_length + 1

        elif pam.position == "5prime":
            pam_seq = rc[i : i + pam.length]
            if not pam.matches(pam_seq):
                continue
            guide = rc[i + pam.length : i + pam.length + guide_length]
            guide_end_plus = base_pos + len(seq) - i - pam.length - 1
            guide_start_plus = guide_end_plus - guide_length + 1
        else:
            continue

        ctx_start = max(0, i - 4)
        ctx_end = min(len(rc), i + guide_length + pam.length + 3)
        context_30 = rc[ctx_start:ctx_end]

        cs_start = max(0, i - 6)
        cs_end = min(len(rc), i + guide_length + pam.length + 6)
        context_35 = rc[cs_start:cs_end]
        if len(context_35) < 35:
            context_35 = ""

        gc = gc_content(guide)
        guide_id = f"{chrom}:{guide_start_plus}-{guide_end_plus}(-)"

        yield SgrnaCandidate(
            guide_id=guide_id,
            guide_sequence=guide,
            pam=pam_seq,
            pam_strand="-",
            chrom=chrom,
            start=guide_start_plus,
            end=guide_end_plus,
            strand="-",
            context_30nt=context_30,
            context_35nt=context_35,
            gc_content=gc,
            nuclease=nuclease,
        )


def compute_guide_features(guide_sequence: str) -> dict[str, float]:
    """Compute sequence-based features for a guide RNA.

    These features are used by efficiency prediction models and for
    filtering suboptimal guides (poly-T, extreme GC, etc.).
    """
    guide = guide_sequence.upper()
    if len(guide) != 20:
        raise ValueError(f"guide must be 20 nt, got {len(guide)}")

    features: dict[str, float] = {}

    # GC content
    features["gc_content"] = gc_content(guide)

    # Poly-T count (TTTT stretches cause Pol III transcription termination)
    features["poly_t_count"] = float(guide.count("TTTT"))

    # Homopolymer runs
    max_run = 1
    current_run = 1
    for i in range(1, len(guide)):
        if guide[i] == guide[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    features["max_homopolymer_run"] = float(max_run)

    # PAM-proximal seed region (positions 1-10 from PAM, i.e., guide 11-20)
    # For SpCas9, the seed is near the 3' end of the guide
    seed = guide[10:20]
    features["seed_gc_content"] = gc_content(seed)

    # PAM-distal region (positions 1-10 of guide)
    distal = guide[0:10]
    features["distal_gc_content"] = gc_content(distal)

    # Thermodynamic stability (simplified Tm estimation)
    # Wallace rule: Tm = 2*(A+T) + 4*(G+C)
    at = guide.count("A") + guide.count("T")
    gc = guide.count("G") + guide.count("C")
    features["wallace_tm"] = float(2 * at + 4 * gc)

    # Energy-based efficiency heuristic (simplified Rule Set 2-like)
    # High efficiency: 40-70% GC, no poly-T, moderate homopolymers
    gc_score = 1.0 - abs(features["gc_content"] - 0.55) * 3.0
    gc_score = max(0.0, min(1.0, gc_score))
    poly_t_penalty = features["poly_t_count"] * 0.15
    homopolymer_penalty = max(0, features["max_homopolymer_run"] - 4) * 0.10

    # Position-dependent weights (simplified from Doench 2016)
    # G at position 20 (PAM-proximal) is favorable
    if guide[19] == "G":
        features["g20_bonus"] = 0.05
    else:
        features["g20_bonus"] = 0.0

    efficiency = gc_score - poly_t_penalty - homopolymer_penalty + features["g20_bonus"]
    features["efficiency_score"] = max(0.0, min(1.0, efficiency))

    return features


def search_pam_sites(
    sequence: str,
    nuclease: NucleaseType = NucleaseType.SPCAS9,
    guide_length: int = 20,
) -> list[tuple[int, str, str]]:
    """Quick PAM site search returning (position, pam, strand) tuples.

    Useful for reconnaissance before full guide design.
    """
    pam = PAM_PATTERNS.get(nuclease)
    if pam is None:
        raise ValueError(f"No PAM pattern for {nuclease.value}")

    sites: list[tuple[int, str, str]] = []
    seq = sequence.upper()

    # Plus strand
    for i in range(len(seq) - pam.length + 1):
        candidate_pam = seq[i : i + pam.length]
        if pam.matches(candidate_pam):
            sites.append((i + 1, candidate_pam, "+"))

    # Minus strand
    rc = reverse_complement(seq)
    for i in range(len(rc) - pam.length + 1):
        candidate_pam = rc[i : i + pam.length]
        if pam.matches(candidate_pam):
            sites.append((len(seq) - i - pam.length + 1, candidate_pam, "-"))

    return sites
