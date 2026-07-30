"""Tests for the genome downloader module."""

import pytest
from geneimpact.genome_downloader import (
    download_sequence,
    list_species,
    SPECIES_TO_ENSEMBL,
    SPECIES_TO_NCBI,
    SPECIES_ASSEMBLY,
    _extract_sequence,
)


class TestSpeciesMapping:
    def test_all_supported_species(self):
        species = list_species()
        assert "mouse" in species
        assert "rat" in species
        assert "zebrafish" in species
        assert "human" in species

    def test_ensembl_mapping(self):
        assert SPECIES_TO_ENSEMBL["mouse"] == "mus_musculus"
        assert SPECIES_TO_ENSEMBL["zebrafish"] == "danio_rerio"

    def test_assembly_mapping(self):
        assert SPECIES_ASSEMBLY["mouse"] == "GRCm39"
        assert SPECIES_ASSEMBLY["human"] == "GRCh38"

    def test_ncbi_mapping(self):
        assert SPECIES_TO_NCBI["mouse"] == "Mus musculus"


class TestExtractSequence:
    def test_basic_extraction(self):
        fasta = ">chr1\nATCGATCGATCG\nATCGATCGATCG\n"
        seq = _extract_sequence(fasta)
        assert seq == "ATCGATCGATCGATCGATCGATCG"

    def test_no_header(self):
        fasta = "ATCGATCG"
        assert _extract_sequence(fasta) == "ATCGATCG"

    def test_uppercase(self):
        fasta = ">chr1\natcgatcg\n"
        assert _extract_sequence(fasta) == "ATCGATCG"


class TestDownloadSequence:
    def test_unsupported_species_raises(self):
        with pytest.raises(ValueError, match="Unsupported species"):
            download_sequence("alien", "1")

    def test_cache_hit(self, tmp_path):
        """Pre-existing cache file should be returned without download."""
        cache_file = tmp_path / "mouse_1.fa"
        cache_file.write_text(">chr1\nATCGATCGATCG\n", encoding="utf-8")
        result = download_sequence(
            "mouse", "1", cache_dir=str(tmp_path)
        )
        assert result.cached is True
        assert result.sequence_length == 12
        assert result.source == "cache"

    def test_unknown_source_raises(self):
        with pytest.raises(ValueError, match="Unknown source"):
            download_sequence("mouse", "1", source="invalid_source")
