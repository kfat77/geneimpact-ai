"""Streaming normalization for public MGI phenotypic allele reports."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Iterator


MGI_COLUMN_COUNT = 13


@dataclass(frozen=True)
class MgiAlleleEvidence:
    """Normalized evidence from one MGI phenotypic allele row."""

    allele_id: str
    allele_symbol: str
    allele_name: str
    allele_type: str
    allele_attributes: tuple[str, ...]
    pubmed_id: str | None
    marker_id: str | None
    marker_symbol: str | None
    refseq_id: str | None
    ensembl_gene_id: str | None
    high_level_mp_ids: tuple[str, ...]
    synonyms: tuple[str, ...]
    marker_name: str | None
    origin_program: str | None

    @property
    def is_endonuclease_mediated(self) -> bool:
        return self.allele_type.casefold() == "endonuclease-mediated"


@dataclass(frozen=True)
class NormalizationSummary:
    """Counts and checksum for a normalized JSONL evidence file."""

    input_sha256: str
    output_sha256: str
    total_records: int
    genome_edited_records: int
    phenotype_annotated_records: int
    output_phenotype_annotated_records: int
    output_records: int


def parse_phenotypic_alleles(lines: Iterable[str]) -> Iterator[MgiAlleleEvidence]:
    """Parse MGI's 13-column report while ignoring comment lines."""
    rows = csv.reader((line for line in lines if line.strip() and not line.startswith("#")), delimiter="\t")
    for row_number, row in enumerate(rows, start=1):
        if len(row) != MGI_COLUMN_COUNT:
            raise ValueError(
                f"MGI data row {row_number} has {len(row)} columns; expected {MGI_COLUMN_COUNT}."
            )
        yield MgiAlleleEvidence(
            allele_id=row[0],
            allele_symbol=row[1],
            allele_name=row[2],
            allele_type=row[3],
            allele_attributes=_split(row[4], "|"),
            pubmed_id=row[5] or None,
            marker_id=row[6] or None,
            marker_symbol=row[7] or None,
            refseq_id=row[8] or None,
            ensembl_gene_id=row[9] or None,
            high_level_mp_ids=_split(row[10], ","),
            synonyms=_split(row[11], "|"),
            marker_name=row[12] or None,
            origin_program=_origin_program(row[1], row[2]),
        )


def normalize_phenotypic_alleles(
    input_path: Path,
    output_path: Path,
    *,
    genome_edited_only: bool = True,
) -> NormalizationSummary:
    """Write normalized JSONL records and return auditable source statistics."""
    input_hash = _sha256(input_path)
    total = genome_edited = phenotype_annotated = output_phenotype_annotated = output_count = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with input_path.open(encoding="utf-8") as source, output_path.open("w", encoding="utf-8", newline="\n") as target:
        for record in parse_phenotypic_alleles(source):
            total += 1
            if record.is_endonuclease_mediated:
                genome_edited += 1
            if record.high_level_mp_ids:
                phenotype_annotated += 1
            if genome_edited_only and not record.is_endonuclease_mediated:
                continue
            target.write(json.dumps(asdict(record), ensure_ascii=False) + "\n")
            output_count += 1
            if record.high_level_mp_ids:
                output_phenotype_annotated += 1
    summary = NormalizationSummary(
        input_sha256=input_hash,
        output_sha256=_sha256(output_path),
        total_records=total,
        genome_edited_records=genome_edited,
        phenotype_annotated_records=phenotype_annotated,
        output_phenotype_annotated_records=output_phenotype_annotated,
        output_records=output_count,
    )
    Path(f"{output_path}.manifest.json").write_text(
        json.dumps(asdict(summary), indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def _split(value: str, separator: str) -> tuple[str, ...]:
    return tuple(item for item in value.split(separator) if item)


def _origin_program(allele_symbol: str, allele_name: str) -> str | None:
    evidence = f"{allele_symbol} {allele_name}".casefold()
    if "(impc)" in evidence or "international mouse phenotyping consortium" in evidence:
        return "IMPC"
    return None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
