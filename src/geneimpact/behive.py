"""Strict normalization of externally executed BE-Hive efficiency predictions.

The original BE-Hive environment is intentionally not bundled here. Its
published implementation uses an older Python stack and its repository does
not declare a standard software license. This adapter validates the documented
input domain and turns a reported result into an auditable, non-design output.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from math import exp, isfinite
import re
from typing import Any, Mapping


BEHIVE_EFFICIENCY_REPOSITORY = "https://github.com/maxwshen/be_predict_efficiency"
BEHIVE_EFFICIENCY_COMMIT = "fbd495910d6c95b24081649015d6257a8badc9d7"
BEHIVE_EFFICIENCY_REFERENCE = (
    f"{BEHIVE_EFFICIENCY_REPOSITORY}/tree/{BEHIVE_EFFICIENCY_COMMIT}"
)
BEHIVE_MOUSE_CELL_TYPE = "mES"
BEHIVE_MOUSE_CONTEXT = "mouse embryonic stem cells (mES)"
BEHIVE_MOUSE_EDITORS = frozenset(
    {
        "ABE",
        "ABE-CP1040",
        "ABE8",
        "AID",
        "BE4",
        "BE4-CP1028",
        "BE4-H47ES48A",
        "CDA",
        "CG-689",
        "CG-APOBEC1",
        "CG-EE",
        "CG-POLD2-APOBEC1-X",
        "CG-RBMX-eA3A-X",
        "CG-RBMX-eA3A-X-HF",
        "CG-X-689-X-RBMX",
        "CG-X-APOBEC1-X-HF",
        "CG-X-EE-X-X",
        "CG-eA3A-dead",
        "eA3A",
        "eA3A_T31A",
        "eA3A_T31AT44A",
        "evoAPOBEC",
    }
)
_COMMIT_PATTERN = re.compile(r"^[0-9a-f]{40}$")
_LOGIT_KEYS = ("Predicted logit score", "predicted_logit_score")
_FRACTION_KEYS = (
    "Predicted fraction of sequenced reads with base editing activity",
    "predicted_fraction",
)


@dataclass(frozen=True)
class BehiveEfficiencyRequest:
    """Declared inputs and provenance for one external BE-Hive execution."""

    sequence: str
    base_editor: str
    cell_type: str
    model_commit: str
    calibration_mean: float | None = None
    calibration_std: float | None = None


@dataclass(frozen=True)
class BehiveEfficiencyPrediction:
    """Normalized efficiency result with explicit biological scope."""

    predictor: str
    predictor_version: str
    task: str
    base_editor: str
    cell_type: str
    species_scope: tuple[str, ...]
    biological_context: str
    sequence_sha256: str
    predicted_logit_score: float
    calibrated_fraction: float | None
    calibration_mean: float | None
    calibration_std: float | None
    evidence_reference: str
    source_license_status: str
    execution_status: str
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class BehiveApplicability:
    """A BE-Hive result plus its fit to an assessment context."""

    prediction: BehiveEfficiencyPrediction
    applicability: str
    note: str


def normalize_behive_efficiency(
    document: Mapping[str, Any],
) -> BehiveEfficiencyPrediction:
    """Validate and normalize a reported BE-Hive efficiency result.

    ``document`` must contain ``request`` and ``raw_output`` objects. The raw
    target sequence is hashed and is not retained in the returned record.
    """
    request_data = _mapping(document, "request")
    raw_output = _mapping(document, "raw_output")
    request = _parse_request(request_data)
    sequence = _validate_request(request)
    logit_score = _finite_number(_first(raw_output, _LOGIT_KEYS), "predicted logit score")
    raw_fraction = _optional_number(_first(raw_output, _FRACTION_KEYS, required=False))

    warnings = [
        "The declared output was produced outside GeneImpact AI and was not independently recomputed.",
        "The model scope is mES culture; it does not establish strain-specific, tissue, germline, or whole-animal effects.",
    ]
    calibrated_fraction: float | None = None
    if request.calibration_mean is None and request.calibration_std is None:
        if raw_fraction is not None:
            warnings.append(
                "A reported edited fraction was omitted because calibration_mean and calibration_std were not declared."
            )
    elif request.calibration_mean is None or request.calibration_std is None:
        raise ValueError("calibration_mean and calibration_std must be declared together.")
    else:
        mean = _finite_number(request.calibration_mean, "calibration_mean")
        std = _finite_number(request.calibration_std, "calibration_std")
        if std <= 0:
            raise ValueError("calibration_std must be positive.")
        calibrated_fraction = _expit(logit_score * std + mean)
        if raw_fraction is not None and abs(raw_fraction - calibrated_fraction) > 1e-6:
            raise ValueError(
                "reported edited fraction does not match the declared calibration parameters."
            )

    return BehiveEfficiencyPrediction(
        predictor="BE-Hive efficiency",
        predictor_version=request.model_commit,
        task="base_editing_efficiency",
        base_editor=request.base_editor,
        cell_type=request.cell_type,
        species_scope=("mouse",),
        biological_context=BEHIVE_MOUSE_CONTEXT,
        sequence_sha256=sha256(sequence.encode("ascii")).hexdigest(),
        predicted_logit_score=logit_score,
        calibrated_fraction=calibrated_fraction,
        calibration_mean=request.calibration_mean,
        calibration_std=request.calibration_std,
        evidence_reference=BEHIVE_EFFICIENCY_REFERENCE,
        source_license_status="not_declared_in_repository",
        execution_status="externally_reported_not_recomputed",
        warnings=tuple(warnings),
    )


def integrate_behive_efficiency(
    predictions: tuple[BehiveEfficiencyPrediction, ...],
    species: str,
    edit_class: str,
) -> tuple[BehiveApplicability, ...]:
    """Label imported predictions without extrapolating beyond declared scope."""
    species_matches = species.casefold() in {"mouse", "mus musculus"}
    edit_matches = edit_class.casefold().replace("-", "_").replace(" ", "_") in {
        "base_editing",
        "base_editor",
    }
    applicability = "declared_match" if species_matches and edit_matches else "out_of_scope"
    note = (
        "Matches the declared mouse/base-editing scope, limited to mES culture."
        if applicability == "declared_match"
        else "Does not match the declared mouse and base-editing scope; not used as applicable evidence."
    )
    return tuple(BehiveApplicability(item, applicability, note) for item in predictions)


def _parse_request(data: Mapping[str, Any]) -> BehiveEfficiencyRequest:
    required = ("sequence", "base_editor", "cell_type", "model_commit")
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"BE-Hive request is missing required fields: {', '.join(missing)}")
    return BehiveEfficiencyRequest(
        sequence=str(data["sequence"]),
        base_editor=str(data["base_editor"]),
        cell_type=str(data["cell_type"]),
        model_commit=str(data["model_commit"]),
        calibration_mean=data.get("calibration_mean"),
        calibration_std=data.get("calibration_std"),
    )


def _validate_request(request: BehiveEfficiencyRequest) -> str:
    sequence = request.sequence.strip().upper()
    if len(sequence) != 50:
        raise ValueError("BE-Hive efficiency sequence must contain exactly 50 nucleotides.")
    if any(base not in "ACGT" for base in sequence):
        raise ValueError("BE-Hive efficiency sequence may contain only A, C, G, and T.")
    if sequence[41:43] != "GG":
        raise ValueError("BE-Hive efficiency sequence must contain an NGG PAM at positions 21-23.")
    if request.cell_type != BEHIVE_MOUSE_CELL_TYPE:
        raise ValueError("the mouse adapter only accepts the BE-Hive mES cell-type model.")
    if request.base_editor not in BEHIVE_MOUSE_EDITORS:
        raise ValueError(f"base editor {request.base_editor!r} is not declared for the mES efficiency model.")
    if not _COMMIT_PATTERN.fullmatch(request.model_commit):
        raise ValueError("model_commit must be a full 40-character lowercase Git commit.")
    if request.model_commit != BEHIVE_EFFICIENCY_COMMIT:
        raise ValueError(
            f"this adapter is verified only for BE-Hive efficiency commit {BEHIVE_EFFICIENCY_COMMIT}."
        )
    return sequence


def _mapping(document: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = document.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object.")
    return value


def _first(
    values: Mapping[str, Any], keys: tuple[str, ...], *, required: bool = True
) -> Any:
    for key in keys:
        if key in values:
            return values[key]
    if required:
        raise ValueError(f"raw_output must include one of: {', '.join(keys)}")
    return None


def _finite_number(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric.") from error
    if not isfinite(number):
        raise ValueError(f"{label} must be finite.")
    return number


def _optional_number(value: Any) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("reported edited fraction must be numeric when present.") from error
    if not isfinite(number):
        return None
    if not 0 <= number <= 1:
        raise ValueError("reported edited fraction must be between 0 and 1.")
    return number


def _expit(value: float) -> float:
    if value >= 0:
        return 1 / (1 + exp(-value))
    exponential = exp(value)
    return exponential / (1 + exponential)
