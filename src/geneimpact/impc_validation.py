"""Build independent, assay-level IMPC validation datasets."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from .impc import ImpcClient, ImpcGenePhenotype


@dataclass(frozen=True)
class ImpcValidationRecord:
    """One tested IMPC outcome with its statistical significance label."""

    gene_symbol: str
    outcome_key: str
    significant: bool
    mp_term_id: str | None
    mp_term_name: str | None
    top_level_mp_terms: tuple[str, ...]
    effect_size: float | None
    p_value: float | None
    procedure_name: str | None
    parameter_name: str | None
    sex: str | None
    zygosity: str | None
    source: str


@dataclass(frozen=True)
class ImpcValidationManifest:
    """Audit summary for a bounded multi-gene IMPC validation set."""

    requested_genes: tuple[str, ...]
    queried_genes: int
    genes_with_results: int
    documents: int
    significant_documents: int
    non_significant_documents: int
    label_semantics: str
    output_sha256: str


def build_impc_validation(
    genes: Iterable[str],
    output_path: Path,
    *,
    client: ImpcClient | None = None,
    max_genes: int = 50,
) -> ImpcValidationManifest:
    """Fetch bounded IMPC results and preserve assay-level labels."""
    unique_genes = tuple(dict.fromkeys(gene.strip() for gene in genes if gene.strip()))
    if not unique_genes:
        raise ValueError("at least one gene is required.")
    if len(unique_genes) > max_genes:
        raise ValueError(f"at most {max_genes} genes may be queried per run.")

    client = client or ImpcClient()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    genes_with_results = documents = significant = 0
    with output_path.open("w", encoding="utf-8", newline="\n") as target:
        for gene in unique_genes:
            evidence = client.gene_phenotypes(gene, significant=None)
            if evidence.results:
                genes_with_results += 1
            for result in evidence.results:
                record = _record(result)
                target.write(json.dumps(asdict(record), ensure_ascii=False) + "\n")
                documents += 1
                significant += int(record.significant)

    manifest = ImpcValidationManifest(
        requested_genes=unique_genes,
        queried_genes=len(unique_genes),
        genes_with_results=genes_with_results,
        documents=documents,
        significant_documents=significant,
        non_significant_documents=documents - significant,
        label_semantics=(
            "significant=true means an IMPC-tested procedure/parameter met the consortium's "
            "statistical significance criteria; false does not mean the gene has no phenotype"
        ),
        output_sha256=_sha256(output_path),
    )
    Path(f"{output_path}.manifest.json").write_text(
        json.dumps(asdict(manifest), indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def _record(result: ImpcGenePhenotype) -> ImpcValidationRecord:
    outcome_parts = (
        result.procedure_name or "unknown-procedure",
        result.parameter_name or "unknown-parameter",
        result.sex or "unknown-sex",
        result.zygosity or "unknown-zygosity",
    )
    return ImpcValidationRecord(
        gene_symbol=result.marker_symbol,
        outcome_key="|".join(outcome_parts),
        significant=result.significant,
        mp_term_id=result.mp_term_id,
        mp_term_name=result.mp_term_name,
        top_level_mp_terms=result.top_level_mp_terms,
        effect_size=result.effect_size,
        p_value=result.p_value,
        procedure_name=result.procedure_name,
        parameter_name=result.parameter_name,
        sex=result.sex,
        zygosity=result.zygosity,
        source="IMPC statistical-result",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
