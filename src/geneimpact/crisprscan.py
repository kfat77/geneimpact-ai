"""Domain-gated implementation of the published CRISPRscan linear score."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from hashlib import sha256
import json
import re
from typing import Any, Mapping

from .species import ZEBRAFISH_PROFILE


CRISPRSCAN_METHOD_REFERENCE = "https://doi.org/10.1038/nmeth.3543"
CRISPRSCORE_VERSION = "1.15.3"
CRISPRSCORE_COMMIT = "cbd6f9f60dc7fb50d14b90485b9561d582caf21e"
CRISPRSCORE_REFERENCE = (
    "https://github.com/crisprVerse/crisprScore/"
    f"tree/{CRISPRSCORE_COMMIT}"
)
CRISPRSCAN_COEFFICIENTS_SHA256 = (
    "6e3f1bbfd58e5426651a15cfd0db6ac2094e0a93158dc51639b5929fc9ced5a4"
)
MAX_GUIDES = 10_000
_GUIDE_ID_PATTERN = re.compile(r"^[A-Za-z0-9._:-]{1,128}$")
_SEQUENCE_PATTERN = re.compile(r"^[ACGT]{35}$")

# Published CRISPRscan coefficients, transcribed from the model supplement via
# crisprScore 1.15.3. Positions are one-based and include the 35-nt context.
_INTERCEPT = 0.183930943629
_COEFFICIENTS: tuple[tuple[str, int, float], ...] = (
    ("AA", 19, -0.0973770966031),
    ("TT", 18, -0.0944240749849),
    ("TT", 13, -0.0861877104467),
    ("CT", 26, -0.0842648932315),
    ("GC", 25, -0.0734536086796),
    ("T", 21, -0.0687304966631),
    ("TG", 23, -0.0663880753745),
    ("AG", 23, -0.0543384556991),
    ("G", 30, -0.0463159143121),
    ("A", 4, -0.0421535206047),
    ("AG", 34, -0.0419359084607),
    ("GA", 34, -0.0377977074249),
    ("A", 18, -0.0338204320002),
    ("C", 25, -0.031648352971),
    ("C", 31, -0.030715556059),
    ("G", 1, -0.0296937089124),
    ("C", 16, -0.0216386089249),
    ("A", 14, -0.0184872292452),
    ("A", 11, -0.0182872920807),
    ("T", 34, -0.0176476916909),
    ("AA", 10, -0.0169054146176),
    ("A", 19, -0.0155764988875),
    ("G", 34, -0.0141671225608),
    ("C", 30, -0.0131827325942),
    ("GA", 31, -0.01227989),
    ("T", 24, -0.0119961724142),
    ("A", 15, -0.0105952955376),
    ("G", 4, -0.00544886928417),
    ("GG", 9, -0.00157798989048),
    ("T", 23, -0.00142224326996),
    ("C", 15, -0.000477727288834),
    ("C", 26, -0.000368972859222),
    ("T", 27, -0.000280844938195),
    ("A", 31, 0.00158974977347),
    ("GT", 18, 0.00239174406352),
    ("C", 9, 0.0024492239247),
    ("GA", 20, 0.00974079889824),
    ("A", 25, 0.0105064050328),
    ("A", 12, 0.0116332345075),
    ("A", 32, 0.0124352308169),
    ("T", 22, 0.0132240349237),
    ("C", 20, 0.0150895135396),
    ("G", 17, 0.0154937801148),
    ("G", 18, 0.016457815932),
    ("T", 30, 0.0172631618141),
    ("A", 13, 0.017628924254),
    ("G", 19, 0.0179168436986),
    ("A", 27, 0.0191268154262),
    ("G", 11, 0.0209290394428),
    ("TG", 3, 0.0229499956202),
    ("GC", 3, 0.0246817853421),
    ("G", 14, 0.0251167144259),
    ("GG", 10, 0.0268021578544),
    ("G", 12, 0.0275911379406),
    ("G", 32, 0.0307124896871),
    ("A", 22, 0.0319309087749),
    ("G", 20, 0.0339570080799),
    ("C", 21, 0.0342629206803),
    ("TT", 17, 0.0349288099016),
    ("T", 13, 0.0354451707029),
    ("G", 26, 0.0361466486748),
    ("A", 24, 0.0374664775295),
    ("C", 22, 0.0376316196793),
    ("G", 16, 0.0379709424474),
    ("GG", 12, 0.041883008595),
    ("TG", 18, 0.0459089908401),
    ("TG", 31, 0.04813681225),
    ("A", 35, 0.0485962591677),
    ("G", 15, 0.0511297166878),
    ("C", 24, 0.0529723136618),
    ("TG", 15, 0.0533728222344),
    ("GT", 11, 0.0536784362118),
    ("GC", 9, 0.0541714023256),
    ("CA", 30, 0.0577598506643),
    ("GT", 24, 0.0609521142769),
    ("G", 13, 0.061360904679),
    ("CA", 24, 0.0622193697604),
    ("AG", 10, 0.0637170933285),
    ("G", 10, 0.0677391822097),
    ("C", 13, 0.0694959439929),
    ("GT", 31, 0.0734253503704),
    ("GG", 13, 0.0743558475817),
    ("C", 27, 0.0799339215437),
    ("G", 27, 0.085151051597),
    ("CC", 21, 0.0889196009093),
    ("CC", 23, 0.0950722864525),
    ("G", 22, 0.101144380251),
    ("G", 24, 0.105488324852),
    ("GT", 23, 0.106718562559),
    ("GG", 25, 0.111559440706),
    ("G", 9, 0.114600681211),
)


@dataclass(frozen=True)
class CrisprscanGuideScore:
    guide_id: str
    context_sha256: str
    score: float
    published_threshold_label: str


@dataclass(frozen=True)
class CrisprscanReport:
    predictor: str
    task: str
    species_profile: str
    genome_build: str
    assembly_accession: str
    reference_strain_or_isolate: str
    nuclease: str
    guide_expression: str
    developmental_context: str
    model_training_assembly: str
    model_training_strain: str
    method_reference: str
    implementation_reference: str
    implementation_version: str
    implementation_commit: str
    coefficients_sha256: str
    request_sha256: str
    guide_count: int
    predictions: tuple[CrisprscanGuideScore, ...]
    warnings: tuple[str, ...]


def score_crisprscan(request: Mapping[str, Any]) -> CrisprscanReport:
    """Score canonical SpCas9 contexts in the declared zebrafish embryo domain."""
    declaration, guides = _validate_request(request)
    predictions = []
    for guide_id, context in guides:
        score = _score_context(context)
        predictions.append(
            CrisprscanGuideScore(
                guide_id=guide_id,
                context_sha256=sha256(context.encode("ascii")).hexdigest(),
                score=score,
                published_threshold_label=_threshold_label(score),
            )
        )
    request_digest = sha256(
        json.dumps(
            request,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()
    return CrisprscanReport(
        predictor="CRISPRscan",
        task="spcas9_guide_activity_ranking",
        species_profile="zebrafish",
        genome_build=declaration["genome_build"],
        assembly_accession=ZEBRAFISH_PROFILE.assembly_accession,
        reference_strain_or_isolate=ZEBRAFISH_PROFILE.reference_strain,
        nuclease="SpCas9",
        guide_expression="t7_in_vitro_transcription",
        developmental_context="zebrafish_embryo",
        model_training_assembly="Zv9",
        model_training_strain="TU",
        method_reference=CRISPRSCAN_METHOD_REFERENCE,
        implementation_reference=CRISPRSCORE_REFERENCE,
        implementation_version=CRISPRSCORE_VERSION,
        implementation_commit=CRISPRSCORE_COMMIT,
        coefficients_sha256=CRISPRSCAN_COEFFICIENTS_SHA256,
        request_sha256=request_digest,
        guide_count=len(predictions),
        predictions=tuple(predictions),
        warnings=(
            "CRISPRscan scores rank expected guide activity; they are not calibrated editing probabilities.",
            "The model was trained using in-vitro-transcribed guides in zebrafish embryos and is not generalized to other species or delivery contexts.",
            "The model was trained on Zv9/TU data; each 35-nt context must be rechecked against the declared current assembly and study-animal variants.",
            "On-target activity does not establish phenotype, off-target risk, edit safety, or animal-welfare acceptability.",
        ),
    )


def _validate_request(
    request: Mapping[str, Any],
) -> tuple[dict[str, str], tuple[tuple[str, str], ...]]:
    required = (
        "species_profile",
        "genome_build",
        "assembly_accession",
        "reference_strain_or_isolate",
        "nuclease",
        "guide_expression",
        "developmental_context",
        "guides",
    )
    missing = [key for key in required if key not in request]
    if missing:
        raise ValueError(f"CRISPRscan request is missing fields: {', '.join(missing)}")
    if request["species_profile"] != "zebrafish":
        raise ValueError("CRISPRscan is enabled only for the declared zebrafish domain.")
    accepted_builds = {
        value.casefold()
        for value in (
            ZEBRAFISH_PROFILE.genome_build,
            *ZEBRAFISH_PROFILE.accepted_build_names,
        )
    }
    genome_build = str(request["genome_build"]).strip()
    if genome_build.casefold() not in accepted_builds:
        raise ValueError(
            f"genome_build must match registered {ZEBRAFISH_PROFILE.genome_build!r}."
        )
    if request["assembly_accession"] != ZEBRAFISH_PROFILE.assembly_accession:
        raise ValueError(
            "assembly_accession must match registered "
            f"{ZEBRAFISH_PROFILE.assembly_accession!r}."
        )
    if (
        str(request["reference_strain_or_isolate"]).strip().casefold()
        != ZEBRAFISH_PROFILE.reference_strain.casefold()
    ):
        raise ValueError(
            "reference_strain_or_isolate must match registered "
            f"{ZEBRAFISH_PROFILE.reference_strain!r}."
        )
    exact_context = {
        "nuclease": "SpCas9",
        "guide_expression": "t7_in_vitro_transcription",
        "developmental_context": "zebrafish_embryo",
    }
    for field, expected in exact_context.items():
        if request[field] != expected:
            raise ValueError(
                f"{field} must be {expected!r} for the declared CRISPRscan domain."
            )
    raw_guides = request["guides"]
    if (
        not isinstance(raw_guides, Sequence)
        or isinstance(raw_guides, (str, bytes))
        or not 1 <= len(raw_guides) <= MAX_GUIDES
    ):
        raise ValueError(f"guides must contain between 1 and {MAX_GUIDES} entries.")
    guides: list[tuple[str, str]] = []
    guide_ids: set[str] = set()
    for index, item in enumerate(raw_guides, start=1):
        if not isinstance(item, Mapping):
            raise ValueError(f"guide {index} must be an object.")
        guide_id = str(item.get("guide_id", ""))
        if not _GUIDE_ID_PATTERN.fullmatch(guide_id):
            raise ValueError(
                f"guide {index} guide_id must use 1-128 safe identifier characters."
            )
        if guide_id in guide_ids:
            raise ValueError(f"duplicate guide_id {guide_id!r}.")
        guide_ids.add(guide_id)
        context = str(item.get("context_35nt", "")).upper()
        if not _SEQUENCE_PATTERN.fullmatch(context):
            raise ValueError(
                f"guide {guide_id!r} context_35nt must contain exactly 35 A/C/G/T bases."
            )
        if context[27:29] != "GG":
            raise ValueError(
                f"guide {guide_id!r} requires a canonical NGG PAM at positions 27-29."
            )
        guides.append((guide_id, context))
    return {"genome_build": genome_build}, tuple(guides)


def _score_context(context: str) -> float:
    score = _INTERCEPT
    for motif, one_based_position, weight in _COEFFICIENTS:
        start = one_based_position - 1
        if context[start : start + len(motif)] == motif:
            score += weight
    return score


def _threshold_label(score: float) -> str:
    if score > 0.70:
        return "above_published_highly_efficient_threshold"
    if score > 0.55:
        return "above_published_efficient_threshold"
    return "below_published_efficient_threshold"
