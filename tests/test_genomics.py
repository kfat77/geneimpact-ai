"""Tests for genomic sequence processing (genomics module)."""

from __future__ import annotations

import pytest
from pathlib import Path

from geneimpact.genomics import (
    FastaReader,
    FastaRecord,
    GenomicInterval,
    reverse_complement,
    gc_content,
    validate_dna_sequence,
    extract_context,
    SequenceStats,
)


class TestReverseComplement:
    def test_basic_complement(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        assert reverse_complement("GAATTC") == "GAATTC"

    def test_lowercase(self):
        assert reverse_complement("atcg") == "cgat"

    def test_empty(self):
        assert reverse_complement("") == ""

    def test_invalid_base_raises(self):
        with pytest.raises(ValueError, match="invalid base"):
            reverse_complement("ATCX")


class TestGcContent:
    def test_all_gc(self):
        assert gc_content("GCGC") == 1.0

    def test_no_gc(self):
        assert gc_content("ATAT") == 0.0

    def test_half(self):
        assert gc_content("ATGC") == 0.5

    def test_empty(self):
        assert gc_content("") == 0.0

    def test_mixed_case(self):
        assert gc_content("AtGc") == 0.5


class TestValidateDnaSequence:
    def test_valid_sequence(self):
        validate_dna_sequence("ATCGATCG")

    def test_invalid_base(self):
        with pytest.raises(ValueError, match="invalid base"):
            validate_dna_sequence("ATCX")

    def test_allows_n(self):
        validate_dna_sequence("ATCGNATC", allowed="ACGTN")


class TestFastaRecord:
    def test_seq_id(self):
        rec = FastaRecord(header="chr1 description here", sequence="ATCG")
        assert rec.seq_id == "chr1"

    def test_length(self):
        rec = FastaRecord(header="chr1", sequence="ATCG")
        assert rec.length == 4

    def test_subsequence(self):
        rec = FastaRecord(header="chr1", sequence="ATCGATCG")
        assert rec.subsequence(2, 5) == "TCGA"

    def test_subsequence_out_of_range(self):
        rec = FastaRecord(header="chr1", sequence="ATCG")
        with pytest.raises(ValueError, match="out of range"):
            rec.subsequence(1, 10)


class TestGenomicInterval:
    def test_valid(self):
        iv = GenomicInterval(chrom="chr1", start=100, end=200)
        assert iv.length == 101

    def test_negative_start(self):
        with pytest.raises(ValueError):
            GenomicInterval(chrom="chr1", start=0, end=10)

    def test_end_before_start(self):
        with pytest.raises(ValueError):
            GenomicInterval(chrom="chr1", start=100, end=50)

    def test_invalid_strand(self):
        with pytest.raises(ValueError):
            GenomicInterval(chrom="chr1", start=1, end=10, strand="*")


class TestFastaReader:
    @pytest.fixture
    def fasta_file(self, tmp_path):
        path = tmp_path / "test.fa"
        path.write_text(
            ">chr1\nATCGATCGATCGATCGATCGATCG\n"
            ">chr2\nGGGGCCCCAAAATTTT\n",
            encoding="utf-8",
        )
        return path

    def test_load_and_access(self, fasta_file):
        reader = FastaReader(fasta_file)
        assert "chr1" in reader
        assert "chr2" in reader
        assert "chr3" not in reader

    def test_fetch(self, fasta_file):
        reader = FastaReader(fasta_file)
        assert reader.fetch("chr1", 1, 4) == "ATCG"

    def test_sequence_ids(self, fasta_file):
        reader = FastaReader(fasta_file)
        assert set(reader.sequence_ids) == {"chr1", "chr2"}

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            FastaReader("nonexistent.fa")


class TestExtractContext:
    def test_basic_extraction(self):
        seq = "AAAATCGATCGATCGATCGATGGAAAA"
        # guide at position 5, length 20, upstream 4, downstream 3
        ctx = extract_context(seq, guide_start=5, guide_length=20, upstream=4, downstream=3)
        assert len(ctx) == 27

    def test_not_enough_upstream(self):
        seq = "ATCGATCG"
        with pytest.raises(ValueError, match="upstream"):
            extract_context(seq, guide_start=1, upstream=4)


class TestSequenceStats:
    def test_basic(self):
        stats = SequenceStats.from_sequence("ATGCATGC")
        assert stats.length == 8
        assert stats.gc_content == 0.5
        assert stats.n_count == 0
        assert stats.base_counts["A"] == 2
