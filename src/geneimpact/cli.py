"""Command-line interface for auditable assessment reports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .datasources import check_ensembl_profile
from .mgi import normalize_phenotypic_alleles
from .snapshots import MGI_REPORTS, create_mgi_snapshot
from .species import PROFILES
from .workflow import DEFAULT_MODEL_VERSION, assess_request


def main() -> None:
    parser = argparse.ArgumentParser(description="Animal genome-edit impact assessment (research use only).")
    subparsers = parser.add_subparsers(dest="command", required=True)
    assess = subparsers.add_parser("assess", help="Create an assessment report from a JSON request.")
    assess.add_argument("request", type=Path, help="Path to an assessment request JSON file.")
    assess.add_argument("--output", "-o", type=Path, help="Write the report to this JSON file.")
    assess.add_argument("--model-version", default=DEFAULT_MODEL_VERSION)
    source_check = subparsers.add_parser(
        "source-check", help="Verify a registered species profile against Ensembl."
    )
    source_check.add_argument("--species", default="mouse", choices=sorted(PROFILES))
    snapshot = subparsers.add_parser(
        "snapshot-mgi", help="Download a versioned MGI report with a checksum manifest."
    )
    snapshot.add_argument("--report", required=True, choices=sorted(MGI_REPORTS))
    snapshot.add_argument("--output-dir", required=True, type=Path)
    normalize = subparsers.add_parser(
        "normalize-mgi", help="Normalize an MGI phenotypic allele report to JSONL."
    )
    normalize.add_argument("--input", required=True, type=Path)
    normalize.add_argument("--output", required=True, type=Path)
    normalize.add_argument(
        "--all-alleles",
        action="store_true",
        help="Include non-endonuclease-mediated alleles.",
    )
    args = parser.parse_args()

    if args.command == "source-check":
        try:
            result = check_ensembl_profile(PROFILES[args.species])
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps({
            "source": result.source,
            "matches": result.matches,
            "checked_release": result.checked_release,
            "errors": list(result.errors),
        }, indent=2))
        return
    if args.command == "snapshot-mgi":
        try:
            manifest = create_mgi_snapshot(args.report, args.output_dir)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps({
            "source": manifest.source,
            "report_key": manifest.report_key,
            "filename": manifest.filename,
            "sha256": manifest.sha256,
            "byte_count": manifest.byte_count,
            "retrieved_at": manifest.retrieved_at,
        }, indent=2))
        return
    if args.command == "normalize-mgi":
        try:
            summary = normalize_phenotypic_alleles(
                args.input,
                args.output,
                genome_edited_only=not args.all_alleles,
            )
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps({
            "input_sha256": summary.input_sha256,
            "output_sha256": summary.output_sha256,
            "total_records": summary.total_records,
            "genome_edited_records": summary.genome_edited_records,
            "phenotype_annotated_records": summary.phenotype_annotated_records,
            "output_phenotype_annotated_records": summary.output_phenotype_annotated_records,
            "output_records": summary.output_records,
            "output": str(args.output),
        }, indent=2))
        return

    try:
        request = json.loads(args.request.read_text(encoding="utf-8"))
        report = assess_request(request, model_version=args.model_version)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        parser.error(str(error))

    rendered = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(f"Assessment report written to {args.output}")
    else:
        print(rendered, end="")
