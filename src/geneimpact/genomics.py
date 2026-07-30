"""Genomic sequence processing: FASTA I/O, reverse complement, sequence retrieval.

This module provides the foundation for all sequence-based gene editing
predictions. It handles FASTA file parsing, sequence manipulation, and
genomic coordinate resolution without external bioinformatics dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Sequence

__all__ = [
    "FastaRecord",
    "FastaReader",
    "GenomicInterval",
    "reverse_complement",
    "gc_content",
    "validate_dna_sequence",
    "extract_context",
]


_COMPLEMENT: dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "a": "t", "t": "a", "c": "g", "g": "c",
    "N": "N", "n": "n",
}

_VALID_BASES = frozenset("ACGTNacgtn")


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record with header and sequence."""

    header: str
    sequence: str

    @property
    def seq_id(self) -> str:
        """Return the first whitespace-delimited token of the header."""
        return self.header.split()[0] if self.header else ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    def subsequence(self, start: int, end: int) -> str:
        """Return 1-based inclusive subsequence [start, end]."""
        if start < 1 or end > self.length or start > end:
            raise ValueError(
                f"interval [{start}, {end}] out of range for sequence "
                f"of length {self.length}"
            )
        return self.sequence[start - 1 : end]


@dataclass(frozen=True)
class GenomicInterval:
    """A genomic interval with 1-based inclusive coordinates."""

    chrom: str
    start: int
    end: int
    strand: str = "+"

    def __post_init__(self) -> None:
        if self.start < 1:
            raise ValueError(f"start must be >= 1, got {self.start}")
        if self.end < self.start:
            raise ValueError(f"end ({self.end}) must be >= start ({self.start})")
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}")

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class FastaReader:
    """Memory-efficient FASTA file reader with chromosome indexing.

    Usage::

        reader = FastaReader("genome.fa")
        chr1 = reader["chr1"]
        seq = chr1.subsequence(1000, 1100)
    """

    def __init__(self, path: str | Path) -> None:
        self._path = Path(path)
        if not self._path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self._path}")
        self._index: dict[str, tuple[int, int]] = {}
        self._cache: dict[str, FastaRecord] = {}
        self._build_index()

    def _build_index(self) -> None:
        """Build a byte-offset index for random-access sequence retrieval."""
        with self._path.open("r", encoding="utf-8") as fh:
            current_id = ""
            current_start = 0
            current_end = 0
            offset = 0
            for line in fh:
                if line.startswith(">"):
                    if current_id:
                        self._index[current_id] = (current_start, current_end)
                    current_id = line[1:].split()[0] if len(line) > 1 else ""
                    current_start = offset + len(line)
                    current_end = current_start
                else:
                    current_end += len(line)
                offset += len(line)
            if current_id:
                self._index[current_id] = (current_start, current_end)

    def _load_record(self, seq_id: str) -> FastaRecord:
        """Load a single record by scanning from its byte offset."""
        if seq_id in self._cache:
            return self._cache[seq_id]
        if seq_id not in self._index:
            raise KeyError(f"sequence {seq_id!r} not found in FASTA index")
        start, end = self._index[seq_id]
        with self._path.open("r", encoding="utf-8") as fh:
            fh.seek(start)
            raw = fh.read(end - start)
        sequence = raw.replace("\n", "").replace("\r", "").upper()
        header = seq_id
        record = FastaRecord(header=header, sequence=sequence)
        self._cache[seq_id] = record
        return record

    def __getitem__(self, seq_id: str) -> FastaRecord:
        return self._load_record(seq_id)

    def __contains__(self, seq_id: str) -> bool:
        return seq_id in self._index

    @property
    def sequence_ids(self) -> tuple[str, ...]:
        return tuple(self._index.keys())

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch a 1-based inclusive genomic interval."""
        record = self._load_record(chrom)
        return record.subsequence(start, end)

    def close(self) -> None:
        self._cache.clear()


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    try:
        return "".join(_COMPLEMENT[base] for base in reversed(sequence))
    except KeyError as exc:
        raise ValueError(f"invalid base in sequence: {exc.args[0]!r}") from exc


def gc_content(sequence: str) -> float:
    """Return GC content as a fraction in [0, 1]."""
    if not sequence:
        return 0.0
    gc = sum(1 for base in sequence.upper() if base in ("G", "C"))
    return gc / len(sequence)


def validate_dna_sequence(sequence: str, allowed: str = "ACGT") -> None:
    """Raise ValueError if the sequence contains unexpected characters."""
    upper = sequence.upper()
    for i, base in enumerate(upper):
        if base not in allowed:
            raise ValueError(
                f"invalid base {base!r} at position {i + 1}; "
                f"allowed: {allowed}"
            )


def extract_context(
    sequence: str,
    guide_start: int,
    guide_length: int = 20,
    upstream: int = 4,
    downstream: int = 3,
) -> str:
    """Extract a context window around a guide RNA.

    For a guide at ``guide_start`` (1-based) of length ``guide_length``,
    returns ``upstream`` bases before the guide, the guide itself, and
    ``downstream`` bases after (typically the PAM).

    The returned context is used for CRISPRscan-style scoring.
    """
    seq = sequence.upper()
    start = guide_start - 1 - upstream
    end = guide_start - 1 + guide_length + downstream
    if start < 0:
        raise ValueError("not enough upstream sequence for context extraction")
    if end > len(seq):
        raise ValueError("not enough downstream sequence for context extraction")
    return seq[start:end]


@dataclass
class SequenceStats:
    """Summary statistics for a DNA sequence."""

    length: int
    gc_content: float
    at_content: float
    n_count: int
    base_counts: dict[str, int] = field(default_factory=dict)

    @classmethod
    def from_sequence(cls, sequence: str) -> "SequenceStats":
        upper = sequence.upper()
        counts: dict[str, int] = {}
        for base in "ACGTN":
            counts[base] = upper.count(base)
        gc = gc_content(upper)
        return cls(
            length=len(upper),
            gc_content=gc,
            at_content=1.0 - gc,
            n_count=counts.get("N", 0),
            base_counts=counts,
        )
