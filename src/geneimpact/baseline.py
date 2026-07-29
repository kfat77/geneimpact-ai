"""Transparent baseline and ranking metrics for positive-association benchmarks."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping


@dataclass(frozen=True)
class PhenotypePriorModel:
    """Global training-set phenotype frequency ranking."""

    name: str
    ranked_phenotype_ids: tuple[str, ...]
    training_associations: int


@dataclass(frozen=True)
class RankingMetrics:
    """Macro metrics for unseen-gene positive-association ranking."""

    genes: int
    associations: int
    k: int
    macro_recall_at_k: float
    gene_hit_rate_at_k: float


@dataclass(frozen=True)
class BaselineReport:
    """Auditable baseline result tied to an exact benchmark manifest."""

    model: PhenotypePriorModel
    benchmark_manifest_sha256: str
    validation: RankingMetrics
    test: RankingMetrics
    calibration_status: str


def fit_phenotype_prior(records: Iterable[Mapping[str, Any]]) -> PhenotypePriorModel:
    """Rank phenotype IDs by positive-association frequency in training data."""
    counts: Counter[str] = Counter()
    total = 0
    for record in records:
        phenotype_id = str(record["phenotype_id"])
        counts[phenotype_id] += 1
        total += 1
    if not counts:
        raise ValueError("training benchmark contains no phenotype associations.")
    ranking = tuple(
        phenotype_id
        for phenotype_id, _ in sorted(
            counts.items(), key=lambda item: (-item[1], item[0])
        )
    )
    return PhenotypePriorModel(
        name="global-phenotype-frequency-prior-v1",
        ranked_phenotype_ids=ranking,
        training_associations=total,
    )


def evaluate_ranking(
    model: PhenotypePriorModel,
    records: Iterable[Mapping[str, Any]],
    *,
    k: int = 5,
) -> RankingMetrics:
    """Evaluate macro recall and any-hit rate for grouped unseen genes."""
    if k < 1:
        raise ValueError("k must be at least 1.")
    expected: dict[str, set[str]] = defaultdict(set)
    for record in records:
        expected[str(record["gene_symbol"])].add(str(record["phenotype_id"]))
    if not expected:
        raise ValueError("evaluation benchmark contains no gene associations.")
    predictions = set(model.ranked_phenotype_ids[:k])
    recalls = [
        len(phenotypes & predictions) / len(phenotypes)
        for phenotypes in expected.values()
    ]
    hits = [bool(phenotypes & predictions) for phenotypes in expected.values()]
    return RankingMetrics(
        genes=len(expected),
        associations=sum(len(phenotypes) for phenotypes in expected.values()),
        k=k,
        macro_recall_at_k=sum(recalls) / len(recalls),
        gene_hit_rate_at_k=sum(hits) / len(hits),
    )


def evaluate_benchmark(
    benchmark_dir: Path,
    *,
    k: int = 5,
) -> BaselineReport:
    """Fit on train and evaluate unchanged on validation and test splits."""
    manifest_path = benchmark_dir / "manifest.json"
    model = fit_phenotype_prior(_read_jsonl(benchmark_dir / "train.jsonl"))
    report = BaselineReport(
        model=model,
        benchmark_manifest_sha256=_sha256(manifest_path),
        validation=evaluate_ranking(
            model, _read_jsonl(benchmark_dir / "validation.jsonl"), k=k
        ),
        test=evaluate_ranking(model, _read_jsonl(benchmark_dir / "test.jsonl"), k=k),
        calibration_status=(
            "not_applicable: benchmark contains observed positive associations only; "
            "do not interpret ranks as probabilities"
        ),
    )
    (benchmark_dir / "baseline-report.json").write_text(
        json.dumps(asdict(report), indent=2) + "\n",
        encoding="utf-8",
    )
    return report


def _read_jsonl(path: Path) -> Iterable[Mapping[str, Any]]:
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, Mapping):
                raise ValueError(f"{path.name} line {line_number} must be an object.")
            yield value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
