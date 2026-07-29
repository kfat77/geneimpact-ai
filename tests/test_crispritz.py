import json

import pytest

from geneimpact.crispritz import CRISPRITZ_COMMIT, import_crispritz_targets


HEADER = (
    "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
    "Direction\tMismatches\tBulge Size\tTotal\n"
)


def _metadata(**overrides):
    metadata = {
        "species_profile": "rat",
        "genome_build": "GRCr8",
        "assembly_accession": "GCF_036323735.1",
        "reference_strain_or_isolate": "BN/NHsdMcwi",
        "edit_class": "knockout",
        "pam_definition": "NNNNNNNNNNNNNNNNNNNNNGG 3",
        "max_mismatches": 4,
        "max_dna_bulge": 1,
        "max_rna_bulge": 1,
        "crispritz_commit": CRISPRITZ_COMMIT,
        "reference_fasta_sha256": "a" * 64,
        "variant_aware": False,
        "variant_snapshot_sha256": None,
    }
    metadata.update(overrides)
    return metadata


def _targets(tmp_path, body=None):
    body = body or (
        "X\tAAAAAAAAAAAAAAAAAAAANNN\tAAAAAAAAAAAAAAAAAAAATGG\tchr1\t10\t10\t+\t0\t0\t0\n"
        "DNA\tAAAAAAAAAAAAAAAAAAAANNN\tAAAAAAAAAAAAAAAAAA-ATGG\tchr2\t20\t20\t-\t1\t1\t2\n"
    )
    path = tmp_path / "example.targets.txt"
    path.write_text(HEADER + body, encoding="utf-8")
    return path


def test_imports_bounded_reference_search_audit(tmp_path):
    report = import_crispritz_targets(_metadata(), _targets(tmp_path))

    assert report.species_profile == "rat"
    assert report.candidate_site_count == 2
    assert report.exact_sequence_match_count == 1
    assert report.observed_guide_count == 1
    assert report.counts_by_total_difference[1].total_differences == 2
    assert report.top_candidate_hits[0].total_differences == 0
    assert len(report.targets_file_sha256) == 64


@pytest.mark.parametrize(
    ("species_profile", "genome_build", "assembly_accession", "reference_context"),
    [
        ("mouse", "GRCm39", "GCF_000001635.27", "C57BL/6J"),
        ("rat", "GRCr8", "GCF_036323735.1", "BN/NHsdMcwi"),
        ("zebrafish", "GRCz12tu", "GCF_049306965.2", "Tuebingen"),
        ("fruit_fly", "Release 6 plus ISO1 MT", "GCF_000001215.4", "ISO-1"),
        ("rhesus_macaque", "T2T-MMU8v2.0", "GCF_049350105.2", "MMU2019108-1"),
        ("cynomolgus_macaque", "T2T-MFA8v1.1", "GCF_037993035.2", "582-1"),
    ],
)
def test_accepts_every_registered_target_species(
    tmp_path, species_profile, genome_build, assembly_accession, reference_context
):
    report = import_crispritz_targets(
        _metadata(
            species_profile=species_profile,
            genome_build=genome_build,
            assembly_accession=assembly_accession,
            reference_strain_or_isolate=reference_context,
        ),
        _targets(tmp_path),
    )

    assert report.species_profile == species_profile
    assert report.assembly_accession == assembly_accession


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"species_profile": "monkey"}, "unknown species"),
        ({"genome_build": "old-build"}, "must match"),
        ({"assembly_accession": "GCF_bad"}, "assembly_accession"),
        ({"reference_strain_or_isolate": "Wistar"}, "reference_strain_or_isolate"),
        ({"crispritz_commit": "bad"}, "verified only"),
        ({"variant_aware": True}, "variant_snapshot"),
    ],
)
def test_rejects_unverifiable_run_metadata(tmp_path, overrides, message):
    with pytest.raises(ValueError, match=message):
        import_crispritz_targets(_metadata(**overrides), _targets(tmp_path))


def test_rejects_rows_inconsistent_with_declared_thresholds(tmp_path):
    path = _targets(
        tmp_path,
        "DNA\tAAAAAAAAAAAAAAAAAAAANNN\tAAAAAAAAAAAAAAAAAA-ATGG\tchr2\t20\t20\t-\t1\t2\t3\n",
    )

    with pytest.raises(ValueError, match="max_dna_bulge"):
        import_crispritz_targets(_metadata(), path)


def test_rejects_inconsistent_total(tmp_path):
    path = _targets(
        tmp_path,
        "X\tAAAAAAAAAAAAAAAAAAAANNN\tAAAAAAAAAAAAAAAAAAAATGG\tchr1\t10\t10\t+\t1\t0\t2\n",
    )

    with pytest.raises(ValueError, match="must equal"):
        import_crispritz_targets(_metadata(), path)


def test_rejects_truncated_result_row(tmp_path):
    path = _targets(
        tmp_path,
        "X\tAAAAAAAAAAAAAAAAAAAANNN\tAAAAAAAAAAAAAAAAAAAATGG\tchr1\t10\n",
    )

    with pytest.raises(ValueError, match="does not match"):
        import_crispritz_targets(_metadata(), path)


def test_metadata_example_is_json_serializable():
    assert json.loads(json.dumps(_metadata()))["species_profile"] == "rat"
