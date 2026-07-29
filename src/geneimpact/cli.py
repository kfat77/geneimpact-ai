"""Command-line interface for auditable assessment reports."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path

from .behive import normalize_behive_efficiency
from .behive_bystander import normalize_behive_bystander
from .behive_validation import evaluate_behive_validation
from .capabilities import capabilities_for_species
from .crispritz import import_crispritz_targets
from .crisprscan import score_crisprscan
from .crisprscan_validation import PREDICTION_FIELDS, evaluate_crisprscan_transfer
from .datasources import check_ncbi_profile
from .dossier import build_research_dossier, verify_dossier_integrity
from .benchmark import build_mgi_benchmark
from .baseline import evaluate_benchmark
from .impc import ImpcClient
from .impc_validation import build_impc_validation
from .impc_calibration import evaluate_impc_calibration
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
        "source-check", help="Verify a registered species profile against its authoritative assembly record."
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
    impc_gene = subparsers.add_parser(
        "impc-gene", help="Fetch significant IMPC phenotype results for one mouse gene."
    )
    impc_gene.add_argument("--gene", required=True)
    impc_gene.add_argument("--output", type=Path)
    benchmark = subparsers.add_parser(
        "benchmark-mgi", help="Build leakage-aware grouped benchmark splits from normalized MGI JSONL."
    )
    benchmark.add_argument("--input", required=True, type=Path)
    benchmark.add_argument("--output-dir", required=True, type=Path)
    benchmark.add_argument(
        "--include-impc-origin",
        action="store_true",
        help="Include MGI rows marked as IMPC-derived; unsafe for independent IMPC validation.",
    )
    baseline = subparsers.add_parser(
        "evaluate-baseline", help="Fit and evaluate the transparent phenotype-frequency baseline."
    )
    baseline.add_argument("--benchmark-dir", required=True, type=Path)
    baseline.add_argument("--k", type=int, default=5)
    impc_validation = subparsers.add_parser(
        "benchmark-impc", help="Build a bounded assay-level IMPC validation dataset."
    )
    impc_validation.add_argument("--gene", action="append", required=True)
    impc_validation.add_argument("--output", required=True, type=Path)
    impc_calibration = subparsers.add_parser(
        "calibrate-impc", help="Evaluate a probability baseline on gene-disjoint IMPC data."
    )
    impc_calibration.add_argument("--calibration", required=True, type=Path)
    impc_calibration.add_argument("--test", required=True, type=Path)
    impc_calibration.add_argument("--output", required=True, type=Path)
    behive_import = subparsers.add_parser(
        "import-behive-efficiency",
        help="Validate and normalize an externally executed BE-Hive efficiency result.",
    )
    behive_import.add_argument("--input", required=True, type=Path)
    behive_import.add_argument("--output", type=Path)
    behive_validation = subparsers.add_parser(
        "validate-behive-efficiency",
        help="Evaluate BE-Hive on a leakage-audited independent mES dataset.",
    )
    behive_validation.add_argument("--input", required=True, type=Path)
    behive_validation.add_argument("--output", required=True, type=Path)
    behive_bystander = subparsers.add_parser(
        "import-behive-bystander",
        help="Validate and normalize externally executed BE-Hive bystander outcomes.",
    )
    behive_bystander.add_argument("--input", required=True, type=Path)
    behive_bystander.add_argument("--output", type=Path)
    capabilities = subparsers.add_parser(
        "capabilities",
        help="Show available and candidate predictors for a registered species.",
    )
    capabilities.add_argument("--species", required=True, choices=sorted(PROFILES))
    crispritz_import = subparsers.add_parser(
        "import-crispritz",
        help="Validate a version-locked external CRISPRitz targets file.",
    )
    crispritz_import.add_argument("--metadata", required=True, type=Path)
    crispritz_import.add_argument("--targets", required=True, type=Path)
    crispritz_import.add_argument("--output", required=True, type=Path)
    crisprscan_score = subparsers.add_parser(
        "score-crisprscan",
        help="Score guides in the declared zebrafish embryo CRISPRscan domain.",
    )
    crisprscan_score.add_argument("--input", required=True, type=Path)
    crisprscan_score.add_argument("--output", required=True, type=Path)
    crisprscan_validation = subparsers.add_parser(
        "validate-crisprscan-transfer",
        help="Evaluate reported CRISPRscan scores on an independent transfer dataset.",
    )
    crisprscan_validation.add_argument("--input", required=True, type=Path)
    crisprscan_validation.add_argument(
        "--prediction-field",
        choices=PREDICTION_FIELDS,
        default="reported_crisprscan_score",
    )
    crisprscan_validation.add_argument("--output", required=True, type=Path)
    dossier = subparsers.add_parser(
        "dossier",
        help="Build one unified, integrity-checkable research dossier.",
    )
    dossier.add_argument("request", type=Path)
    dossier.add_argument("--output", "-o", required=True, type=Path)
    dossier.add_argument("--model-version", default=DEFAULT_MODEL_VERSION)
    verify_dossier = subparsers.add_parser(
        "verify-dossier",
        help="Verify a dossier content hash (not an authenticity signature).",
    )
    verify_dossier.add_argument("report", type=Path)
    args = parser.parse_args()

    if args.command == "source-check":
        try:
            result = check_ncbi_profile(PROFILES[args.species])
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
    if args.command == "impc-gene":
        try:
            evidence = ImpcClient().significant_gene_phenotypes(args.gene)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(evidence), indent=2, ensure_ascii=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"IMPC evidence written to {args.output}")
        else:
            print(rendered, end="")
        return
    if args.command == "benchmark-mgi":
        try:
            manifest = build_mgi_benchmark(
                args.input,
                args.output_dir,
                include_impc_origin=args.include_impc_origin,
            )
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps(asdict(manifest), indent=2))
        return
    if args.command == "evaluate-baseline":
        try:
            report = evaluate_benchmark(args.benchmark_dir, k=args.k)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps(asdict(report), indent=2))
        return
    if args.command == "benchmark-impc":
        try:
            manifest = build_impc_validation(args.gene, args.output)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps(asdict(manifest), indent=2))
        return
    if args.command == "calibrate-impc":
        try:
            report = evaluate_impc_calibration(
                args.calibration,
                args.test,
                output_path=args.output,
            )
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps(asdict(report), indent=2))
        return
    if args.command == "import-behive-efficiency":
        try:
            document = json.loads(args.input.read_text(encoding="utf-8"))
            prediction = normalize_behive_efficiency(document)
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(prediction), indent=2, ensure_ascii=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"BE-Hive audit record written to {args.output}")
        else:
            print(rendered, end="")
        return
    if args.command == "validate-behive-efficiency":
        try:
            report = evaluate_behive_validation(args.input)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(report), indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"BE-Hive validation report written to {args.output}")
        return
    if args.command == "import-behive-bystander":
        try:
            document = json.loads(args.input.read_text(encoding="utf-8"))
            prediction = normalize_behive_bystander(document)
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(prediction), indent=2, ensure_ascii=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"BE-Hive bystander audit record written to {args.output}")
        else:
            print(rendered, end="")
        return
    if args.command == "capabilities":
        rows = [
            {
                **asdict(item),
                "status": item.status.value,
            }
            for item in capabilities_for_species(args.species)
        ]
        print(json.dumps(rows, indent=2, ensure_ascii=False))
        return
    if args.command == "import-crispritz":
        try:
            metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
            if not isinstance(metadata, dict):
                raise ValueError("CRISPRitz metadata must be a JSON object.")
            report = import_crispritz_targets(metadata, args.targets)
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(report), indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"CRISPRitz audit report written to {args.output}")
        return
    if args.command == "score-crisprscan":
        try:
            request = json.loads(args.input.read_text(encoding="utf-8"))
            if not isinstance(request, dict):
                raise ValueError("CRISPRscan request must be a JSON object.")
            report = score_crisprscan(request)
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(report), indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"CRISPRscan score report written to {args.output}")
        return
    if args.command == "validate-crisprscan-transfer":
        try:
            dataset = json.loads(args.input.read_text(encoding="utf-8"))
            if not isinstance(dataset, dict):
                raise ValueError("CRISPRscan validation dataset must be a JSON object.")
            report = evaluate_crisprscan_transfer(
                dataset,
                prediction_field=args.prediction_field,
            )
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(report), indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"CRISPRscan transfer report written to {args.output}")
        return
    if args.command == "dossier":
        try:
            request_bytes = args.request.read_bytes()
            request = json.loads(request_bytes.decode("utf-8"))
            if not isinstance(request, dict):
                raise ValueError("dossier request must be a JSON object.")
            report = build_research_dossier(
                request,
                attachment_base_dir=args.request.parent,
                source_request_sha256=sha256(request_bytes).hexdigest(),
                model_version=args.model_version,
            )
        except (
            OSError,
            UnicodeDecodeError,
            json.JSONDecodeError,
            TypeError,
            ValueError,
        ) as error:
            parser.error(str(error))
        rendered = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"Research dossier written to {args.output}")
        return
    if args.command == "verify-dossier":
        try:
            report = json.loads(args.report.read_text(encoding="utf-8"))
            if not isinstance(report, dict):
                raise ValueError("dossier report must be a JSON object.")
            verification = verify_dossier_integrity(report)
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
            parser.error(str(error))
        print(json.dumps(asdict(verification), indent=2, ensure_ascii=False))
        if not verification.matches:
            raise SystemExit(1)
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
