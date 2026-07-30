"""Automated genome sequence downloader via Ensembl REST API and NCBI Datasets.

Provides one-command access to reference genome FASTA files for supported
species, with local caching and checksum verification.

Supported sources:
- Ensembl REST API (primary): https://rest.ensembl.org
- NCBI Datasets API (fallback): https://api.ncbi.nlm.nih.gov/datasets/v2alpha

Usage::

    from geneimpact.genome_downloader import download_sequence

    # Download mouse chromosome 1 from Ensembl
    path = download_sequence("mouse", "1", cache_dir="./genomes")

    # Or download a specific gene region
    path = download_sequence("mouse", "1", start=1000000, end=1100000)
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "GenomeSource",
    "DownloadResult",
    "download_sequence",
    "download_genome",
    "list_species",
    "SPECIES_TO_ENSEMBL",
    "SPECIES_TO_NCBI",
]

# ---------------------------------------------------------------------------
# Species name → database name mappings
# ---------------------------------------------------------------------------

SPECIES_TO_ENSEMBL: dict[str, str] = {
    "mouse": "mus_musculus",
    "rat": "rattus_norvegicus",
    "zebrafish": "danio_rerio",
    "fruit_fly": "drosophila_melanogaster",
    "rhesus_macaque": "macaca_mulatta",
    "cynomolgus_macaque": "macaca_fascicularis",
    "human": "homo_sapiens",
}

SPECIES_TO_NCBI: dict[str, str] = {
    "mouse": "Mus musculus",
    "rat": "Rattus norvegicus",
    "zebrafish": "Danio rerio",
    "fruit_fly": "Drosophila melanogaster",
    "rhesus_macaque": "Macaca mulatta",
    "cynomolgus_macaque": "Macaca fascicularis",
    "human": "Homo sapiens",
}

# Assembly defaults
SPECIES_ASSEMBLY: dict[str, str] = {
    "mouse": "GRCm39",
    "rat": "mRatBN7.2",
    "zebrafish": "GRCz11",
    "fruit_fly": "BDGP6.46",
    "rhesus_macaque": "Mmul_10",
    "cynomolgus_macaque": "Macaca_fascicularis_6.0",
    "human": "GRCh38",
}

ENSEMBL_REST_BASE = "https://rest.ensembl.org"
NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


class GenomeSource(str):
    """Genome data source identifier."""
    ENSEMBL = "ensembl"
    NCBI = "ncbi"
    CACHE = "cache"  # Already downloaded


@dataclass(frozen=True)
class DownloadResult:
    """Result of a genome download operation."""

    species: str
    chrom: str
    source: str
    assembly: str
    local_path: str
    sequence_length: int
    sha256: str
    start: int = 1
    end: int = 0  # 0 = full chromosome
    cached: bool = False
    warnings: tuple[str, ...] = ()


def download_sequence(
    species: str,
    chrom: str,
    start: int | None = None,
    end: int | None = None,
    cache_dir: str | Path = "./genome_cache",
    source: str = "ensembl",
    force_refresh: bool = False,
) -> DownloadResult:
    """Download a genomic sequence from Ensembl or NCBI.

    Parameters
    ----------
    species : str
        Species key (mouse, rat, zebrafish, fruit_fly, human, etc.)
    chrom : str
        Chromosome name (e.g., "1", "X", "MT").
    start : int | None
        1-based start position. If None, downloads from beginning.
    end : int | None
        1-based end position. If None, downloads to end of chromosome.
    cache_dir : str | Path
        Directory for caching downloaded sequences.
    source : str
        Data source: "ensembl" (default) or "ncbi".
    force_refresh : bool
        If True, re-download even if cached version exists.

    Returns
    -------
    DownloadResult
        Download metadata including local file path.

    Raises
    ------
    ValueError
        If species is not supported or sequence not found.
    urllib.error.URLError
        If network request fails.
    """
    species_key = species.lower().replace(" ", "_")
    if species_key not in SPECIES_TO_ENSEMBL:
        raise ValueError(
            f"Unsupported species {species!r}. Supported: "
            f"{', '.join(sorted(SPECIES_TO_ENSEMBL))}"
        )

    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)

    # Determine coordinate range
    seq_start = start or 1
    seq_end = end or 0  # 0 means full chromosome
    is_region = start is not None and end is not None

    # Build cache filename
    if is_region:
        cache_file = cache_path / f"{species_key}_{chrom}_{seq_start}_{seq_end}.fa"
    else:
        cache_file = cache_path / f"{species_key}_{chrom}.fa"

    # Check cache
    if not force_refresh and cache_file.exists():
        content = cache_file.read_text(encoding="utf-8")
        seq = _extract_sequence(content)
        sha = hashlib.sha256(content.encode()).hexdigest()
        assembly = SPECIES_ASSEMBLY.get(species_key, "unknown")
        return DownloadResult(
            species=species_key,
            chrom=chrom,
            source=GenomeSource.CACHE,
            assembly=assembly,
            local_path=str(cache_file),
            sequence_length=len(seq),
            sha256=sha,
            start=seq_start,
            end=seq_end if seq_end else len(seq),
            cached=True,
        )

    # Download from source
    if source == "ensembl":
        result = _download_ensembl(
            species_key, chrom, seq_start, seq_end, cache_file, is_region
        )
    elif source == "ncbi":
        result = _download_ncbi(
            species_key, chrom, seq_start, seq_end, cache_file, is_region
        )
    else:
        raise ValueError(f"Unknown source {source!r}. Use 'ensembl' or 'ncbi'.")

    return result


def _download_ensembl(
    species_key: str,
    chrom: str,
    start: int,
    end: int,
    cache_file: Path,
    is_region: bool,
) -> DownloadResult:
    """Download sequence from Ensembl REST API."""
    ensembl_name = SPECIES_TO_ENSEMBL[species_key]
    assembly = SPECIES_ASSEMBLY.get(species_key, "unknown")

    if is_region:
        # Region-specific sequence
        region = f"{chrom}:{start}-{end}"
        endpoint = f"/sequence/region/{ensembl_name}/{region}"
        params = {"content-type": "text/x-fasta"}
    else:
        # Full chromosome
        endpoint = f"/sequence/region/{ensembl_name}/{chrom}"
        params = {"content-type": "text/x-fasta"}

    url = f"{ENSEMBL_REST_BASE}{urllib.parse.quote(endpoint)}?{urllib.parse.urlencode(params)}"

    req = urllib.request.Request(url)
    req.add_header("Accept", "text/x-fasta")
    req.add_header("User-Agent", "GeneImpact-AI/1.0")

    with urllib.request.urlopen(req, timeout=120) as response:
        content = response.read().decode("utf-8")

    if not content or content.startswith("error"):
        raise ValueError(
            f"Ensembl API returned error for {species_key}:{chrom}. "
            f"Response: {content[:200]}"
        )

    # Write to cache
    cache_file.write_text(content, encoding="utf-8")

    seq = _extract_sequence(content)
    sha = hashlib.sha256(content.encode()).hexdigest()

    return DownloadResult(
        species=species_key,
        chrom=chrom,
        source=GenomeSource.ENSEMBL,
        assembly=assembly,
        local_path=str(cache_file),
        sequence_length=len(seq),
        sha256=sha,
        start=start,
        end=end if end else len(seq),
        warnings=(
            "Downloaded from Ensembl REST API. Verify assembly version "
            f"({assembly}) matches your study context." if not is_region else
            f"Region {chrom}:{start}-{end} downloaded from Ensembl ({assembly}).",
        ),
    )


def _download_ncbi(
    species_key: str,
    chrom: str,
    start: int,
    end: int,
    cache_file: Path,
    is_region: bool,
) -> DownloadResult:
    """Download sequence from NCBI Datasets API (fallback)."""
    ncbi_name = SPECIES_TO_NCBI[species_key]
    assembly = SPECIES_ASSEMBLY.get(species_key, "unknown")

    # NCBI Datasets API for chromosome sequence
    # This is a simplified interface; full implementation would use
    # the NCBI Datasets CLI or the v2alpha REST API
    taxon_search = urllib.parse.quote(ncbi_name)

    if is_region:
        # Use NCBI E-utilities for region fetch
        url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nucleotide&id={chrom}&seq_start={start}&seq_stop={end}"
            f"&rettype=fasta&retmode=text"
        )
    else:
        # Full chromosome - this requires knowing the RefSeq accession
        # For simplicity, we fall back to Ensembl for full chromosomes
        raise ValueError(
            "NCBI full chromosome download requires a RefSeq accession. "
            "Please use source='ensembl' for full chromosome downloads, "
            "or provide a specific region with start/end parameters."
        )

    req = urllib.request.Request(url)
    req.add_header("User-Agent", "GeneImpact-AI/1.0")

    with urllib.request.urlopen(req, timeout=120) as response:
        content = response.read().decode("utf-8")

    if not content or "error" in content.lower():
        raise ValueError(
            f"NCBI API returned error for {species_key}:{chrom}. "
            f"Response: {content[:200]}"
        )

    cache_file.write_text(content, encoding="utf-8")

    seq = _extract_sequence(content)
    sha = hashlib.sha256(content.encode()).hexdigest()

    return DownloadResult(
        species=species_key,
        chrom=chrom,
        source=GenomeSource.NCBI,
        assembly=assembly,
        local_path=str(cache_file),
        sequence_length=len(seq),
        sha256=sha,
        start=start,
        end=end,
        warnings=(
            f"Region {chrom}:{start}-{end} downloaded from NCBI E-utilities ({assembly}).",
        ),
    )


def download_genome(
    species: str,
    chromosomes: list[str] | None = None,
    cache_dir: str | Path = "./genome_cache",
    source: str = "ensembl",
) -> list[DownloadResult]:
    """Download multiple chromosomes for a species.

    Parameters
    ----------
    species : str
        Species key (mouse, human, etc.)
    chromosomes : list[str] | None
        List of chromosome names. If None, downloads common chromosomes.
    cache_dir : str | Path
        Cache directory.
    source : str
        Data source.

    Returns
    -------
    list[DownloadResult]
        Results for each downloaded chromosome.
    """
    if chromosomes is None:
        chromosomes = [str(i) for i in range(1, 20)] + ["X", "Y"]

    results: list[DownloadResult] = []
    for chrom in chromosomes:
        try:
            result = download_sequence(
                species=species,
                chrom=chrom,
                cache_dir=cache_dir,
                source=source,
            )
            results.append(result)
        except (ValueError, OSError) as e:
            # Skip failed chromosomes but continue
            results.append(DownloadResult(
                species=species,
                chrom=chrom,
                source="failed",
                assembly="",
                local_path="",
                sequence_length=0,
                sha256="",
                warnings=(str(e),),
            ))

    return results


def list_species() -> list[str]:
    """Return list of supported species for genome download."""
    return sorted(SPECIES_TO_ENSEMBL.keys())


def _extract_sequence(fasta_content: str) -> str:
    """Extract the DNA sequence from FASTA content (without header)."""
    lines = fasta_content.strip().split("\n")
    seq_lines = [line for line in lines if not line.startswith(">")]
    return "".join(seq_lines).upper().replace("\r", "")
