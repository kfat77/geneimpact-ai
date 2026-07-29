import json

import pytest

from geneimpact.mgi import normalize_phenotypic_alleles, parse_phenotypic_alleles


EDITED = (
    "MGI:1\tGene<em1>\tedit 1\tEndonuclease-mediated\tNull/knockout\t"
    "\tMGI:2\tGene\tNM_1\tENSMUSG1\tMP:1,MP:2\tGene<old>\tgene name\n"
)
TARGETED = (
    "MGI:3\tGene<tm1>\ttargeted 1\tTargeted\tReporter\t123\tMGI:2\t"
    "Gene\tNM_1\tENSMUSG1\tMP:3\t\tgene name\n"
)
IMPC_EDITED = EDITED.replace("Gene<em1>", "Gene<em1(IMPC)Mbp>")


def test_parser_normalizes_multivalue_fields():
    record = next(parse_phenotypic_alleles([EDITED]))

    assert record.is_endonuclease_mediated
    assert record.allele_attributes == ("Null/knockout",)
    assert record.high_level_mp_ids == ("MP:1", "MP:2")
    assert record.origin_program is None


def test_parser_marks_impc_origin():
    assert next(parse_phenotypic_alleles([IMPC_EDITED])).origin_program == "IMPC"


def test_normalizer_filters_to_endonuclease_mediated_records(tmp_path):
    source = tmp_path / "MGI_PhenotypicAllele.rpt"
    output = tmp_path / "normalized.jsonl"
    source.write_text("# comment\n" + EDITED + TARGETED, encoding="utf-8")

    summary = normalize_phenotypic_alleles(source, output)
    records = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]

    assert summary.total_records == 2
    assert summary.genome_edited_records == 1
    assert summary.output_records == 1
    assert summary.output_phenotype_annotated_records == 1
    assert records[0]["allele_id"] == "MGI:1"
    manifest = json.loads(
        (tmp_path / "normalized.jsonl.manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["output_sha256"] == summary.output_sha256


def test_parser_rejects_schema_drift():
    with pytest.raises(ValueError, match="expected 13"):
        list(parse_phenotypic_alleles(["too\tfew\n"]))
