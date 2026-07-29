"""Reproducible snapshots of public mouse evidence reports."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable
from urllib.request import Request, urlopen


@dataclass(frozen=True)
class ReportSpec:
    """A fixed public report endpoint and its local filename."""

    key: str
    source_url: str
    filename: str
    description: str


@dataclass(frozen=True)
class SnapshotManifest:
    """Audit metadata for one downloaded source file."""

    source: str
    report_key: str
    source_url: str
    retrieved_at: str
    filename: str
    sha256: str
    byte_count: int


MGI_REPORTS = {
    "all-phenotypes": ReportSpec(
        key="all-phenotypes",
        source_url="https://www.informatics.jax.org/downloads/reports/ALL_Phenotype.rpt",
        filename="ALL_Phenotype.rpt",
        description="All mouse genotypes and Mammalian Phenotype annotations.",
    ),
    "phenotypic-alleles": ReportSpec(
        key="phenotypic-alleles",
        source_url="https://www.informatics.jax.org/downloads/reports/MGI_PhenotypicAllele.rpt",
        filename="MGI_PhenotypicAllele.rpt",
        description="Mouse alleles with marker and high-level phenotype annotations.",
    ),
}


def create_mgi_snapshot(
    report_key: str,
    output_dir: Path,
    downloader: Callable[[str], bytes] | None = None,
) -> SnapshotManifest:
    """Download a fixed MGI report and write a checksum-bearing manifest."""
    try:
        spec = MGI_REPORTS[report_key]
    except KeyError as error:
        raise ValueError(f"unknown MGI report {report_key!r}.") from error

    content = (downloader or _download)(spec.source_url)
    if not content:
        raise ValueError(f"MGI report {report_key!r} was empty.")

    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / spec.filename
    report_path.write_bytes(content)
    manifest = SnapshotManifest(
        source="Mouse Genome Informatics",
        report_key=spec.key,
        source_url=spec.source_url,
        retrieved_at=datetime.now(timezone.utc).isoformat(),
        filename=spec.filename,
        sha256=hashlib.sha256(content).hexdigest(),
        byte_count=len(content),
    )
    manifest_path = output_dir / f"{spec.filename}.manifest.json"
    manifest_path.write_text(
        json.dumps(asdict(manifest), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def _download(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "GeneImpact-AI/0.3"})
    with urlopen(request, timeout=60) as response:
        return response.read()
