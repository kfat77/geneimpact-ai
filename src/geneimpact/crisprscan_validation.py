"""Retrospective independent transfer evaluation for CRISPRscan scores."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence


PREDICTION_FIELDS = (
    "reported_crisprscan_score",
    "reported_crisprscan_nhgri_score",
)
_SHA256_LENGTH = 64
_MIN_RECORDS = 20
_MIN_GENES = 5
_MAX_RECORDS = 10_000


@dataclass(frozen=True)
class CrisprscanTransferReport:
    predictor: str
    evaluation_status: str
    domain_fit: str
    prediction_field: str
    source_doi: str
    source_workbook_sha256: str
    training_source_workbook_sha256: str
    training_exact_guide_overlap_count: int
    training_reverse_complement_overlap_count: int
    study_strain: str
    study_assembly: str
    guide_expression: str
    delivery: str
    observed_assay: str
    record_count: int
    gene_count: int
    pearson_correlation: float
    pearson_ci95_lower: float
    pearson_ci95_upper: float
    spearman_correlation: float
    within_gene_pair_count: int
    within_gene_concordant_pair_count: int
    within_gene_pairwise_accuracy: float
    warnings: tuple[str, ...]


def evaluate_crisprscan_transfer(
    dataset: Mapping[str, Any],
    prediction_field: str = "reported_crisprscan_score",
) -> CrisprscanTransferReport:
    """Evaluate reported scores on an independent, explicitly transfer-domain dataset."""
    if prediction_field not in PREDICTION_FIELDS:
        raise ValueError(
            f"prediction_field must be one of: {', '.join(PREDICTION_FIELDS)}."
        )
    metadata, records = _validate_dataset(dataset, prediction_field)
    predictions = [record[prediction_field] for record in records]
    observations = [record["observed_indel_fraction"] for record in records]
    pearson = _pearson(predictions, observations)
    ci_lower, ci_upper = _pearson_ci95(pearson, len(records))
    spearman = _pearson(_average_ranks(predictions), _average_ranks(observations))
    pair_count, concordant_count = _within_gene_concordance(
        records, prediction_field
    )
    return CrisprscanTransferReport(
        predictor="CRISPRscan",
        evaluation_status="retrospective_external_transfer_benchmark",
        domain_fit="outside_declared_t7_expression_domain",
        prediction_field=prediction_field,
        source_doi=metadata["source_doi"],
        source_workbook_sha256=metadata["source_workbook_sha256"],
        training_source_workbook_sha256=metadata[
            "training_source_workbook_sha256"
        ],
        training_exact_guide_overlap_count=metadata[
            "training_exact_guide_overlap_count"
        ],
        training_reverse_complement_overlap_count=metadata[
            "training_reverse_complement_overlap_count"
        ],
        study_strain=metadata["study_strain"],
        study_assembly=metadata["study_assembly"],
        guide_expression=metadata["guide_expression"],
        delivery=metadata["delivery"],
        observed_assay=metadata["observed_assay"],
        record_count=len(records),
        gene_count=len({record["gene"] for record in records}),
        pearson_correlation=pearson,
        pearson_ci95_lower=ci_lower,
        pearson_ci95_upper=ci_upper,
        spearman_correlation=spearman,
        within_gene_pair_count=pair_count,
        within_gene_concordant_pair_count=concordant_count,
        within_gene_pairwise_accuracy=concordant_count / pair_count,
        warnings=(
            "This dataset is independent but uses crRNA:tracrRNA SpCas9 RNP delivery, not the adapter's declared T7 in-vitro-transcribed guide domain.",
            "The benchmark uses scores reported by the source study and does not independently validate GeneImpact's 35-nt sequence extraction.",
            "This is a retrospective analysis performed after inspecting the dataset; it is not a preregistered pass/fail test.",
            "Correlation and ranking agreement do not establish phenotype prediction, off-target safety, or animal-welfare acceptability.",
        ),
    )


def _validate_dataset(
    dataset: Mapping[str, Any], prediction_field: str
) -> tuple[Mapping[str, Any], tuple[dict[str, Any], ...]]:
    metadata = dataset.get("metadata")
    raw_records = dataset.get("records")
    if not isinstance(metadata, Mapping):
        raise ValueError("dataset metadata must be an object.")
    if (
        not isinstance(raw_records, Sequence)
        or isinstance(raw_records, (str, bytes))
        or not _MIN_RECORDS <= len(raw_records) <= _MAX_RECORDS
    ):
        raise ValueError(
            f"records must contain between {_MIN_RECORDS} and {_MAX_RECORDS} entries."
        )
    required_metadata = (
        "source_doi",
        "source_workbook_sha256",
        "independent_lab",
        "training_data_overlap",
        "training_source_workbook_sha256",
        "training_guide_count",
        "training_unique_guide_count",
        "training_exact_guide_overlap_count",
        "training_reverse_complement_overlap_count",
        "training_overlap_audit_method",
        "species_profile",
        "study_strain",
        "study_assembly",
        "guide_expression",
        "delivery",
        "observed_assay",
        "benchmark_scope",
    )
    missing_metadata = [key for key in required_metadata if key not in metadata]
    if missing_metadata:
        raise ValueError(
            "dataset metadata is missing fields: " + ", ".join(missing_metadata)
        )
    if metadata["independent_lab"] is not True:
        raise ValueError("independent_lab must be true.")
    if metadata["training_data_overlap"] is not False:
        raise ValueError("training_data_overlap must be false.")
    if not _is_sha256(metadata["training_source_workbook_sha256"]):
        raise ValueError(
            "training_source_workbook_sha256 must be a lowercase SHA-256 digest."
        )
    for field in (
        "training_guide_count",
        "training_unique_guide_count",
        "training_exact_guide_overlap_count",
        "training_reverse_complement_overlap_count",
    ):
        if isinstance(metadata[field], bool) or not isinstance(metadata[field], int):
            raise ValueError(f"{field} must be an integer.")
    if metadata["training_guide_count"] < _MIN_RECORDS:
        raise ValueError("training_guide_count is too small for a training overlap audit.")
    if not 0 < metadata["training_unique_guide_count"] <= metadata["training_guide_count"]:
        raise ValueError("training_unique_guide_count is inconsistent.")
    if (
        metadata["training_exact_guide_overlap_count"] != 0
        or metadata["training_reverse_complement_overlap_count"] != 0
    ):
        raise ValueError("training guide overlap counts must both be zero.")
    if not str(metadata["training_overlap_audit_method"]).strip():
        raise ValueError("training_overlap_audit_method is required.")
    if metadata["species_profile"] != "zebrafish":
        raise ValueError("CRISPRscan transfer evaluation requires zebrafish data.")
    if metadata["benchmark_scope"] != "retrospective_external_transfer":
        raise ValueError(
            "benchmark_scope must disclose retrospective_external_transfer."
        )
    digest = metadata["source_workbook_sha256"]
    if not _is_sha256(digest):
        raise ValueError("source_workbook_sha256 must be a lowercase SHA-256 digest.")
    records: list[dict[str, Any]] = []
    record_ids: set[str] = set()
    guide_hashes: set[str] = set()
    for index, item in enumerate(raw_records, start=1):
        if not isinstance(item, Mapping):
            raise ValueError(f"record {index} must be an object.")
        required_record = (
            "record_id",
            "gene",
            "guide_sha256",
            prediction_field,
            "observed_indel_fraction",
        )
        missing_record = [key for key in required_record if key not in item]
        if missing_record:
            raise ValueError(
                f"record {index} is missing fields: {', '.join(missing_record)}"
            )
        record_id = str(item["record_id"]).strip()
        gene = str(item["gene"]).strip()
        guide_hash = item["guide_sha256"]
        if not record_id or record_id in record_ids:
            raise ValueError(f"record {index} has a blank or duplicate record_id.")
        if not gene:
            raise ValueError(f"record {index} gene is required.")
        if not _is_sha256(guide_hash) or guide_hash in guide_hashes:
            raise ValueError(f"record {index} has an invalid or duplicate guide_sha256.")
        prediction = _fraction(item[prediction_field], prediction_field, index)
        observation = _fraction(
            item["observed_indel_fraction"], "observed_indel_fraction", index
        )
        record_ids.add(record_id)
        guide_hashes.add(guide_hash)
        records.append(
            {
                "record_id": record_id,
                "gene": gene,
                "guide_sha256": guide_hash,
                prediction_field: prediction,
                "observed_indel_fraction": observation,
            }
        )
    if len({record["gene"] for record in records}) < _MIN_GENES:
        raise ValueError(f"dataset must contain at least {_MIN_GENES} genes.")
    pair_count, _ = _within_gene_concordance(tuple(records), prediction_field)
    if pair_count == 0:
        raise ValueError("dataset has no comparable within-gene guide pairs.")
    return metadata, tuple(records)


def _fraction(value: Any, label: str, index: int) -> float:
    if isinstance(value, bool):
        raise ValueError(f"record {index} {label} must be numeric.")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"record {index} {label} must be numeric.") from error
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise ValueError(f"record {index} {label} must be between 0 and 1.")
    return number


def _pearson(x: Sequence[float], y: Sequence[float]) -> float:
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    centered = [(a - mean_x, b - mean_y) for a, b in zip(x, y, strict=True)]
    numerator = sum(a * b for a, b in centered)
    denominator = math.sqrt(
        sum(a * a for a, _ in centered) * sum(b * b for _, b in centered)
    )
    if denominator == 0:
        raise ValueError("correlation is undefined for a constant series.")
    return numerator / denominator


def _pearson_ci95(correlation: float, sample_size: int) -> tuple[float, float]:
    bounded = min(max(correlation, -0.999999999999), 0.999999999999)
    fisher_z = math.atanh(bounded)
    margin = 1.959963984540054 / math.sqrt(sample_size - 3)
    return math.tanh(fisher_z - margin), math.tanh(fisher_z + margin)


def _average_ranks(values: Sequence[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][1] == ordered[start][1]:
            end += 1
        average_rank = (start + 1 + end) / 2
        for ordered_index in range(start, end):
            ranks[ordered[ordered_index][0]] = average_rank
        start = end
    return ranks


def _within_gene_concordance(
    records: Sequence[Mapping[str, Any]], prediction_field: str
) -> tuple[int, int]:
    groups: dict[str, list[Mapping[str, Any]]] = {}
    for record in records:
        groups.setdefault(str(record["gene"]), []).append(record)
    comparable = 0
    concordant = 0
    for group in groups.values():
        for left_index, left in enumerate(group):
            for right in group[left_index + 1 :]:
                prediction_delta = (
                    left[prediction_field] - right[prediction_field]
                )
                observation_delta = (
                    left["observed_indel_fraction"]
                    - right["observed_indel_fraction"]
                )
                if prediction_delta == 0 or observation_delta == 0:
                    continue
                comparable += 1
                concordant += int(
                    math.copysign(1, prediction_delta)
                    == math.copysign(1, observation_delta)
                )
    return comparable, concordant


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == _SHA256_LENGTH
        and all(character in "0123456789abcdef" for character in value)
    )
