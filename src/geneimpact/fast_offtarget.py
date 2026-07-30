"""Fast off-target detection using k-mer seed-and-extend algorithm.

Replaces the O(n) brute-force scan with a hash-based seed lookup
that reduces the search space to O(n/4^k) for exact seed matches,
then extends candidates with full mismatch verification.

Algorithm:
1. Build a hash index of k-mer seeds from the reference genome
2. For each guide, extract seed k-mers from the PAM-proximal region
3. Look up seed positions in the hash index (with optional seed mismatches)
4. Extend each candidate to full guide length and verify PAM + mismatch count
5. Score and rank confirmed off-target sites

For a 10-nt seed with SpCas9 (NGG PAM), the expected reduction is
~4^10 = ~1,000,000x fewer positions to check per guide.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Sequence

from .genomics import reverse_complement
from .sgrna_design import NucleaseType, PAM_PATTERNS, PamPattern
from .offtarget import (
    OffTargetSite,
    OffTargetReport,
    MismatchPattern,
    score_offtarget,
    _risk_level,
)

__all__ = [
    "SeedIndex",
    "fast_find_offtargets",
    "build_seed_index",
    "DEFAULT_SEED_LENGTH",
    "DEFAULT_SEED_MISMATCHES",
]

DEFAULT_SEED_LENGTH = 10
DEFAULT_SEED_MISMATCHES = 0  # Exact seed match by default; increase for sensitivity


@dataclass
class SeedIndex:
    """Hash index of k-mer seeds for fast off-target lookup.

    Maps each k-mer to a list of (chrom, position, strand) tuples
    where that k-mer appears in the reference genome.
    """

    kmer_length: int
    index: dict[str, list[tuple[str, int, str]]] = field(default_factory=dict)
    total_positions: int = 0
    reference_length: int = 0

    def lookup(self, kmer: str) -> list[tuple[str, int, str]]:
        """Return all positions where this k-mer appears."""
        return self.index.get(kmer.upper(), [])


def build_seed_index(
    reference_sequences: dict[str, str],
    seed_length: int = DEFAULT_SEED_LENGTH,
) -> SeedIndex:
    """Build a k-mer seed index from reference sequences.

    Parameters
    ----------
    reference_sequences : dict[str, str]
        Dictionary of {seq_id: sequence} to index.
    seed_length : int
        Length of k-mer seeds to index (default 10).

    Returns
    -------
    SeedIndex
        Hash index for fast seed lookup.
    """
    idx = SeedIndex(kmer_length=seed_length)
    idx.reference_length = sum(len(s) for s in reference_sequences.values())

    for chrom, seq in reference_sequences.items():
        ref = seq.upper()
        # Index plus strand
        for i in range(len(ref) - seed_length + 1):
            kmer = ref[i:i + seed_length]
            if "N" in kmer:
                continue
            if kmer not in idx.index:
                idx.index[kmer] = []
            idx.index[kmer].append((chrom, i, "+"))
            idx.total_positions += 1

        # Index minus strand (reverse complement)
        rc = reverse_complement(ref)
        for i in range(len(rc) - seed_length + 1):
            kmer = rc[i:i + seed_length]
            if "N" in kmer:
                continue
            if kmer not in idx.index:
                idx.index[kmer] = []
            idx.index[kmer].append((chrom, i, "-"))
            idx.total_positions += 1

    return idx


def fast_find_offtargets(
    guide_sequence: str,
    reference_sequences: dict[str, str] | None = None,
    seed_index: SeedIndex | None = None,
    nuclease: NucleaseType = NucleaseType.SPCAS9,
    max_mismatches: int = 4,
    seed_length: int = DEFAULT_SEED_LENGTH,
    seed_mismatches: int = DEFAULT_SEED_MISMATCHES,
    guide_length: int = 20,
) -> OffTargetReport:
    """Find off-target sites using k-mer seed-and-extend algorithm.

    This function is significantly faster than brute-force scanning for
    large reference sequences. For small sequences (<10 kb), the brute-force
    approach may be faster due to lower overhead.

    Parameters
    ----------
    guide_sequence : str
        20-nt guide RNA sequence (without PAM).
    reference_sequences : dict[str, str] | None
        Reference genome sequences. Required if seed_index is None.
    seed_index : SeedIndex | None
        Pre-built seed index. If None, one will be built from reference_sequences.
    nuclease : NucleaseType
        CRISPR nuclease type for PAM detection.
    max_mismatches : int
        Maximum total mismatches allowed (0-6).
    seed_length : int
        Length of seed k-mer for indexing (default 10).
    seed_mismatches : int
        Number of mismatches allowed in the seed region (0-2).
        Higher values increase sensitivity but reduce speed.

    Returns
    -------
    OffTargetReport
        Off-target analysis results (same format as find_offtargets).
    """
    guide = guide_sequence.upper()
    if len(guide) != guide_length:
        raise ValueError(f"guide must be {guide_length} nt, got {len(guide)}")

    if max_mismatches < 0 or max_mismatches > 6:
        raise ValueError(f"max_mismatches must be 0-6, got {max_mismatches}")

    # Build or use seed index
    if seed_index is None:
        if reference_sequences is None:
            raise ValueError(
                "Either reference_sequences or seed_index must be provided"
            )
        seed_index = build_seed_index(reference_sequences, seed_length)
    elif reference_sequences is None:
        raise ValueError(
            "reference_sequences required even when seed_index is provided "
            "(for PAM verification and sequence extraction)"
        )

    pam = PAM_PATTERNS[nuclease]
    warnings: list[str] = []

    # Extract seed k-mers from the guide
    # For SpCas9, the seed region is PAM-proximal (last 10-12 nt of guide)
    # We extract all possible seed_length-mers from the guide
    seed_kmers: set[str] = set()
    guide_seeds: list[tuple[str, int]] = []  # (kmer, offset_in_guide)

    for offset in range(guide_length - seed_length + 1):
        kmer = guide[offset:offset + seed_length]
        if "N" not in kmer:
            seed_kmers.add(kmer)
            guide_seeds.append((kmer, offset))

    # For inexact seed matching, generate k-mers with 1-2 substitutions
    if seed_mismatches > 0:
        expanded_seeds: set[str] = set()
        for kmer in seed_kmers:
            expanded_seeds.add(kmer)
            if seed_mismatches >= 1:
                for i in range(seed_length):
                    for base in "ACGT":
                        if base != kmer[i]:
                            mutated = kmer[:i] + base + kmer[i+1:]
                            expanded_seeds.add(mutated)
            if seed_mismatches >= 2:
                for i in range(seed_length):
                    for j in range(i + 1, seed_length):
                        for b1 in "ACGT":
                            if b1 == kmer[i]:
                                continue
                            for b2 in "ACGT":
                                if b2 == kmer[j]:
                                    continue
                                mutated = (
                                    kmer[:i] + b1 + kmer[i+1:j] + b2 + kmer[j+1:]
                                )
                                expanded_seeds.add(mutated)
        seed_kmers = expanded_seeds

    # Look up seed positions and collect candidate sites
    # Key: (chrom, position, strand) → set of (guide_seed_offset, mismatches_in_seed)
    candidates: dict[tuple[str, int, str], set[int]] = {}

    for kmer, guide_offset in guide_seeds:
        # For exact seed match
        if seed_mismatches == 0:
            positions = seed_index.lookup(kmer)
            for chrom, pos, strand in positions:
                key = (chrom, pos, strand)
                if key not in candidates:
                    candidates[key] = set()
                candidates[key].add(guide_offset)
        else:
            # For inexact seed match, look up all expanded seeds
            for expanded_kmer in _expand_seed_with_mismatches(kmer, seed_mismatches):
                positions = seed_index.lookup(expanded_kmer)
                for chrom, pos, strand in positions:
                    key = (chrom, pos, strand)
                    if key not in candidates:
                        candidates[key] = set()
                    candidates[key].add(guide_offset)

    # Verify candidates: check PAM, count full-guide mismatches
    off_targets: list[OffTargetSite] = []
    total_scanned = 0

    for (chrom, seed_pos, strand), guide_offsets in candidates.items():
        ref_seq = reference_sequences.get(chrom, "")
        if not ref_seq:
            continue
        ref = ref_seq.upper()

        for guide_offset in guide_offsets:
            # Determine the full guide position in the reference
            if strand == "+":
                guide_start_ref = seed_pos - guide_offset
            else:
                # On minus strand: seed_pos is in reverse-complement coordinates
                # Map back to plus strand
                guide_start_ref = len(ref) - seed_pos - seed_length - (guide_length - guide_offset - seed_length)

            # Check bounds
            if guide_start_ref < 0:
                continue

            if pam.position == "3prime":
                total_window_start = guide_start_ref
                total_window_end = guide_start_ref + guide_length + pam.length
            elif pam.position == "5prime":
                total_window_start = guide_start_ref - pam.length
                total_window_end = guide_start_ref + guide_length
            else:
                continue

            if total_window_start < 0 or total_window_end > len(ref):
                continue

            # Extract the target sequence and PAM
            if strand == "+":
                if pam.position == "3prime":
                    target = ref[guide_start_ref:guide_start_ref + guide_length]
                    pam_seq = ref[guide_start_ref + guide_length:guide_start_ref + guide_length + pam.length]
                    ref_start = guide_start_ref + 1  # 1-based
                    ref_end = ref_start + guide_length - 1
                else:
                    pam_seq = ref[total_window_start:total_window_start + pam.length]
                    target = ref[total_window_start + pam.length:total_window_start + pam.length + guide_length]
                    ref_start = total_window_start + pam.length + 1
                    ref_end = ref_start + guide_length - 1
            else:
                # Minus strand: use reverse complement
                rc = reverse_complement(ref)
                # Calculate position in RC coordinates
                rc_start = len(ref) - total_window_end
                rc_end = len(ref) - total_window_start
                rc_window = rc[rc_start:rc_end]

                if pam.position == "3prime":
                    target = rc_window[:guide_length]
                    pam_seq = rc_window[guide_length:guide_length + pam.length]
                else:
                    pam_seq = rc_window[:pam.length]
                    target = rc_window[pam.length:pam.length + guide_length]

                ref_start = total_window_start + 1
                ref_end = total_window_start + guide_length

            # Skip if N present
            if "N" in target or "N" in pam_seq:
                continue

            # Check PAM
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

            # Skip exact match (on-target)
            if mm_count == 0:
                continue

            total_scanned += 1

            mm_pattern = MismatchPattern(
                positions=tuple(mismatches),
                count=mm_count,
            )

            score = score_offtarget(guide, target, mm_pattern, nuclease)
            risk = _risk_level(mm_count, score, mm_pattern)

            off_targets.append(OffTargetSite(
                chrom=chrom,
                start=ref_start,
                end=ref_end,
                strand=strand,
                off_target_sequence=target,
                pam=pam_seq,
                mismatch_count=mm_count,
                mismatch_positions=tuple(mismatches),
                mismatch_pattern=mm_pattern,
                score=score,
                risk_level=risk,
            ))

    # Sort by score (descending) then by mismatch count (ascending)
    off_targets.sort(key=lambda ot: (-ot.score, ot.mismatch_count))

    # Limit results
    if len(off_targets) > 1000:
        warnings.append(
            f"Found {len(off_targets)} off-target sites; "
            "returning top 1000 by score."
        )
        off_targets = off_targets[:1000]

    if not off_targets and not warnings:
        warnings.append(
            "No off-target sites found. This may indicate high guide specificity "
            "or insufficient reference sequence coverage."
        )

    return OffTargetReport(
        guide_sequence=guide,
        nuclease=nuclease,
        total_sites_scanned=total_scanned,
        off_targets=off_targets,
        max_mismatches_searched=max_mismatches,
        warnings=warnings,
    )


def _expand_seed_with_mismatches(
    seed: str,
    max_mismatches: int,
) -> set[str]:
    """Generate all k-mers within max_mismatches of the seed."""
    result = {seed}
    if max_mismatches == 0:
        return result

    for i in range(len(seed)):
        for base in "ACGT":
            if base != seed[i]:
                mutated = seed[:i] + base + seed[i+1:]
                result.add(mutated)
                if max_mismatches >= 2:
                    for j in range(i + 1, len(seed)):
                        for b2 in "ACGT":
                            if b2 != seed[j]:
                                double_mut = (
                                    mutated[:j] + b2 + mutated[j+1:]
                                )
                                result.add(double_mut)
    return result
