"""Calibrated baseline for the bounded IMPC assay-level task."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from .calibration import brier_score, expected_calibration_error


@dataclass(frozen=True)
class ConstantSignificanceModel:
    """Jeffreys-smoothed significance prevalence estimated on calibration genes."""

    name: str
    probability: float
    calibration_documents: int
    calibration_positives: int


@dataclass(frozen=True)
class BinaryCalibrationMetrics:
    """Probability metrics for an untouched group of genes."""

    genes: int
    documents: int
    positives: int
    prevalence: float
    brier_score: float
    expected_calibration_error: float


@dataclass(frozen=True)
class ImpcCalibrationReport:
    """Audit-bound calibration baseline and task semantics."""

    model: ConstantSignificanceModel
    calibration_sha256: str
    test_sha256: str
    test: BinaryCalibrationMetrics
    gene_overlap_checked: bool
    task_semantics: str


def evaluate_impc_calibration(
    calibration_path: Path,
    test_path: Path,
    *,
    output_path: Path | None = None,
    bins: int = 10,
) -> ImpcCalibrationReport:
    """Fit a prevalence baseline on calibration genes and score untouched genes."""
    calibration_genes, calibration_labels = _read_dataset(calibration_path)
    test_genes, test_labels = _read_dataset(test_path)
    overlap = calibration_genes & test_genes
    if overlap:
        preview = ", ".join(sorted(overlap)[:5])
        raise ValueError(f"gene leakage between calibration and test datasets: {preview}")

    positives = sum(calibration_labels)
    probability = (positives + 0.5) / (len(calibration_labels) + 1)
    model = ConstantSignificanceModel(
        name="jeffreys-smoothed-impc-prevalence-v1",
        probability=probability,
        calibration_documents=len(calibration_labels),
        calibration_positives=positives,
    )
    predictions = [probability] * len(test_labels)
    test_positives = sum(test_labels)
    report = ImpcCalibrationReport(
        model=model,
        calibration_sha256=_sha256(calibration_path),
        test_sha256=_sha256(test_path),
        test=BinaryCalibrationMetrics(
            genes=len(test_genes),
            documents=len(test_labels),
            positives=test_positives,
            prevalence=test_positives / len(test_labels),
            brier_score=brier_score(predictions, test_labels),
            expected_calibration_error=expected_calibration_error(
                predictions, test_labels, bins=bins
            ),
        ),
        gene_overlap_checked=True,
        task_semantics=(
            "Predict whether a specific IMPC knockout procedure/parameter comparison is "
            "statistically significant; this is not a probability that an edit is safe"
        ),
    )
    destination = output_path or test_path.with_suffix(".calibration-report.json")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(asdict(report), indent=2) + "\n",
        encoding="utf-8",
    )
    return report


def _read_dataset(path: Path) -> tuple[set[str], list[int]]:
    genes: set[str] = set()
    labels: list[int] = []
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, Mapping):
                raise ValueError(f"{path.name} line {line_number} must be an object.")
            genes.add(str(value["gene_symbol"]))
            labels.append(int(bool(value["significant"])))
    if not labels:
        raise ValueError(f"{path.name} contains no IMPC results.")
    return genes, labels


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
