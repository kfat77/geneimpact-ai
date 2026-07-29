import hashlib
import json

import pytest

from geneimpact.snapshots import create_mgi_snapshot


def test_mgi_snapshot_writes_report_and_checksum_manifest(tmp_path):
    content = b"allele\tphenotype\nMGI:1\tMP:1\n"
    manifest = create_mgi_snapshot(
        "all-phenotypes", tmp_path, downloader=lambda _: content
    )

    assert (tmp_path / "ALL_Phenotype.rpt").read_bytes() == content
    saved = json.loads(
        (tmp_path / "ALL_Phenotype.rpt.manifest.json").read_text(encoding="utf-8")
    )
    assert saved["sha256"] == hashlib.sha256(content).hexdigest()
    assert manifest.byte_count == len(content)


def test_empty_snapshot_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="empty"):
        create_mgi_snapshot("all-phenotypes", tmp_path, downloader=lambda _: b"")
