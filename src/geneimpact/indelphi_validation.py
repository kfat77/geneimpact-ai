"""Pinned external mouse-embryo transfer evidence for inDelphi."""

from __future__ import annotations

from dataclasses import dataclass


INDELPHI_MOUSE_VALIDATION_DOI = "10.1038/s42003-026-09771-z"
INDELPHI_MOUSE_VALIDATION_REFERENCE = (
    "https://www.nature.com/articles/s42003-026-09771-z"
)
INDELPHI_MOUSE_SUPPLEMENT_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs42003-026-09771-z/MediaObjects/"
    "42003_2026_9771_MOESM4_ESM.xlsx"
)
INDELPHI_MOUSE_SUPPLEMENT_SHA256 = (
    "1dee32aa2c4145454b1b73fb4dc10bc46abf83ba4e4a3be068ddc4403a2c8156"
)


@dataclass(frozen=True)
class IndelphiMouseTransferEvidence:
    evidence_status: str
    study_reference: str
    doi: str
    source_supplement_url: str
    source_supplement_sha256: str
    source_license: str
    source_rows_redistributed: bool
    species: str
    strain: str
    gene: str
    developmental_context: str
    delivery_context: str
    guide_count: int
    compared_outcome_pairs: int
    overall_pearson_r: float
    very_strong_guide_count: int
    strong_guide_count: int
    moderate_guide_count: int
    weak_guide_count: int
    interpretation: str
    limitations: tuple[str, ...]


def indelphi_mouse_transfer_evidence() -> IndelphiMouseTransferEvidence:
    """Return independently checked aggregate evidence without copying source rows."""
    return IndelphiMouseTransferEvidence(
        evidence_status="retrospective_external_mouse_embryo_transfer",
        study_reference=INDELPHI_MOUSE_VALIDATION_REFERENCE,
        doi=INDELPHI_MOUSE_VALIDATION_DOI,
        source_supplement_url=INDELPHI_MOUSE_SUPPLEMENT_URL,
        source_supplement_sha256=INDELPHI_MOUSE_SUPPLEMENT_SHA256,
        source_license="CC BY-NC-ND 4.0",
        source_rows_redistributed=False,
        species="Mus musculus",
        strain="C57BL/6JJmsSlc",
        gene="Tyr",
        developmental_context="pooled blastocysts",
        delivery_context="HiFi-Cas9 RNP electroporation into one-cell embryos",
        guide_count=14,
        compared_outcome_pairs=1182,
        overall_pearson_r=0.6375914615783601,
        very_strong_guide_count=1,
        strong_guide_count=5,
        moderate_guide_count=3,
        weak_guide_count=5,
        interpretation=(
            "The aggregate mutation-frequency association was moderate and varied "
            "substantially by guide; this supports prioritization, not direct embryo "
            "outcome assurance."
        ),
        limitations=(
            "Wild-type reads were excluded and remaining edited reads were normalized to 100%.",
            "The primary transfer comparison used 14 Tyr guides and pooled blastocysts.",
            "The study found experimental mESC profiles more concordant with blastocysts than inDelphi alone.",
            "This is retrospective external transfer evidence, not prospective validation for a new locus, strain, delivery protocol, or laboratory.",
        ),
    )
