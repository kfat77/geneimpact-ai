"""Leakage-aware benchmark construction from normalized public evidence."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping


DEFAULT_SPLIT_SEED = "geneimpact-mouse-v1"


@dataclass(frozen=True)
class BenchmarkRecord:
    """One positive gene-edit/phenotype association for ranking evaluation."""

    record_id: str
    gene_symbol: str
    ensembl_gene_id: str | None
    allele_id: str
    edit_class: str
    allele_attributes: tuple[str, ...]
    phenotype_id: str
    source: str
    source_program: str | None
    split: str


@dataclass(frozen=True)
class BenchmarkManifest:
    """Reproducibility and leakage audit for a generated benchmark."""

    input_sha256: str
    split_seed: str
    split_strategy: str
    label_semantics: str
    include_impc_origin: bool
    source_records: int
    eligible_alleles: int
    excluded_impc_alleles: int
    output_associations: int
    split_associations: dict[str, int]
    split_genes: dict[str, int]
    output_sha256: dict[str, str]


def assign_gene_split(gene_symbol: str, seed: str = DEFAULT_SPLIT_SEED) -> str:
    """Assign all records for a gene to one deterministic split."""
    digest = hashlib.sha256(f"{seed}:{gene_symbol.casefold()}".encode()).digest()
    bucket = int.from_bytes(digest[:4], "big") % 10_000
    if bucket < 7_000:
        return "train"
    if bucket < 8_500:
        return "validation"
    return "test"


def build_mgi_benchmark(
    input_path: Path,
    output_dir: Path,
    *,
    include_impc_origin: bool = False,
    split_seed: str = DEFAULT_SPLIT_SEED,
) -> BenchmarkManifest:
    """Create grouped positive-association splits from normalized MGI JSONL."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {split: output_dir / f"{split}.jsonl" for split in ("train", "validation", "test")}
    handles = {split: path.open("w", encoding="utf-8", newline="\n") for split, path in paths.items()}
    source_records = eligible_alleles = excluded_impc = output_associations = 0
    split_associations = {split: 0 for split in paths}
    split_gene_sets = {split: set() for split in paths}

    try:
        with input_path.open(encoding="utf-8") as source:
            for line_number, line in enumerate(source, start=1):
                if not line.strip():
                    continue
                raw = _object(json.loads(line), line_number)
                source_records += 1
                gene_symbol = _optional_text(raw.get("marker_symbol"))
                phenotype_ids = raw.get("high_level_mp_ids", [])
                if not gene_symbol or not isinstance(phenotype_ids, list) or not phenotype_ids:
                    continue
                eligible_alleles += 1
                origin_program = _optional_text(raw.get("origin_program"))
                if origin_program == "IMPC" and not include_impc_origin:
                    excluded_impc += 1
                    continue
                split = assign_gene_split(gene_symbol, split_seed)
                split_gene_sets[split].add(gene_symbol)
                for phenotype_id in phenotype_ids:
                    record = BenchmarkRecord(
                        record_id=f"MGI:{raw['allele_id']}:{phenotype_id}",
                        gene_symbol=gene_symbol,
                        ensembl_gene_id=_optional_text(raw.get("ensembl_gene_id")),
                        allele_id=str(raw["allele_id"]),
                        edit_class="endonuclease-mediated",
                        allele_attributes=tuple(str(value) for value in raw.get("allele_attributes", [])),
                        phenotype_id=str(phenotype_id),
                        source="MGI",
                        source_program=origin_program,
                        split=split,
                    )
                    handles[split].write(json.dumps(asdict(record), ensure_ascii=False) + "\n")
                    output_associations += 1
                    split_associations[split] += 1
    finally:
        for handle in handles.values():
            handle.close()

    _assert_disjoint(split_gene_sets)
    output_hashes = {split: _sha256(path) for split, path in paths.items()}
    manifest = BenchmarkManifest(
        input_sha256=_sha256(input_path),
        split_seed=split_seed,
        split_strategy="SHA-256 grouped by gene symbol; 70% train, 15% validation, 15% test",
        label_semantics="Observed positive high-level Mammalian Phenotype association; missing is not negative",
        include_impc_origin=include_impc_origin,
        source_records=source_records,
        eligible_alleles=eligible_alleles,
        excluded_impc_alleles=excluded_impc,
        output_associations=output_associations,
        split_associations=split_associations,
        split_genes={split: len(genes) for split, genes in split_gene_sets.items()},
        output_sha256=output_hashes,
    )
    (output_dir / "manifest.json").write_text(
        json.dumps(asdict(manifest), indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def _object(value: Any, line_number: int) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"normalized MGI line {line_number} must be an object.")
    return value


def _optional_text(value: Any) -> str | None:
    return None if value is None or not str(value).strip() else str(value)


def _assert_disjoint(split_gene_sets: Mapping[str, set[str]]) -> None:
    train, validation, test = (
        split_gene_sets["train"],
        split_gene_sets["validation"],
        split_gene_sets["test"],
    )
    if train & validation or train & test or validation & test:
        raise AssertionError("gene leakage detected across benchmark splits.")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
