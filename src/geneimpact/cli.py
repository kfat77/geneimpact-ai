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
from .housden import normalize_housden
from .fruit_fly_cas12a import (
    PORT_2026_CAS12A_SOURCE,
    audit_fruit_fly_cas12a_evidence,
    lookup_fruit_fly_cas12a_array,
)
from .benchmark import build_mgi_benchmark
from .baseline import evaluate_benchmark
from .impc import ImpcClient
from .impc_validation import build_impc_validation
from .impc_calibration import evaluate_impc_calibration
from .indelphi import normalize_indelphi
from .mgi import normalize_phenotypic_alleles
from .readiness import readiness_for_species, readiness_matrix
from .rat_validation import (
    evaluate_rat_guide_transfer,
    prepare_rat_guide_transfer_template,
)
from .snapshots import MGI_REPORTS, create_mgi_snapshot
from .species import PROFILES
from .workflow import DEFAULT_MODEL_VERSION, assess_request
from .sgrna_design import NucleaseType, design_sgrnas
from .offtarget import find_offtargets
from .efficiency import predict_efficiency, predict_indel_outcomes
from .pipeline import PipelineConfig, run_pipeline
from .provenance import StudyContext
from .visualization import generate_html_report
from .genomics import FastaReader


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
    indelphi_import = subparsers.add_parser(
        "import-indelphi",
        help="Validate a version-locked external inDelphi mESC result.",
    )
    indelphi_import.add_argument("--input", required=True, type=Path)
    indelphi_import.add_argument("--output", type=Path)
    housden_import = subparsers.add_parser(
        "import-housden",
        help="Validate an official FlyRNAi Housden result envelope.",
    )
    housden_import.add_argument("--input", required=True, type=Path)
    housden_import.add_argument(
        "--source-response",
        required=True,
        type=Path,
        help="Retained XLS response downloaded from the official FlyRNAi service.",
    )
    housden_import.add_argument("--output", type=Path)
    capabilities = subparsers.add_parser(
        "capabilities",
        help="Show available and candidate predictors for a registered species.",
    )
    capabilities.add_argument("--species", required=True, choices=sorted(PROFILES))
    readiness = subparsers.add_parser(
        "readiness",
        help="Show qualified evidence maturity without promoting hazard observations.",
    )
    readiness_scope = readiness.add_mutually_exclusive_group(required=True)
    readiness_scope.add_argument("--species", choices=sorted(PROFILES))
    readiness_scope.add_argument("--all", action="store_true")
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
    rat_template = subparsers.add_parser(
        "prepare-rat-guide-transfer",
        help="Prepare a sequence-redacted prediction template from pinned rat sources.",
    )
    rat_template.add_argument("--table1", required=True, type=Path)
    rat_template.add_argument("--table5", required=True, type=Path)
    rat_template.add_argument("--output", required=True, type=Path)
    rat_validation = subparsers.add_parser(
        "validate-rat-guide-transfer",
        help="Evaluate external guide-activity scores on the bounded rat benchmark.",
    )
    rat_validation.add_argument("--table1", required=True, type=Path)
    rat_validation.add_argument("--table5", required=True, type=Path)
    rat_validation.add_argument("--predictions", required=True, type=Path)
    rat_validation.add_argument("--output", required=True, type=Path)
    fruit_fly_cas12a = subparsers.add_parser(
        "audit-fruit-fly-cas12a-evidence",
        help=(
            "Audit pinned Port 2026 in-vivo Cas12a array-level LOH evidence."
        ),
    )
    fruit_fly_cas12a.add_argument("--library", required=True, type=Path)
    fruit_fly_cas12a.add_argument("--genotypes", required=True, type=Path)
    fruit_fly_cas12a.add_argument("--source-data", required=True, type=Path)
    fruit_fly_cas12a.add_argument(
        "--line-id",
        help="Optionally include one indivisible HD12aCFD array lookup.",
    )
    fruit_fly_cas12a.add_argument("--output", required=True, type=Path)
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

    # --- New prediction pipeline commands ---

    # design-sgrna: Design sgRNAs from a sequence file
    design_sgrna_cmd = subparsers.add_parser(
        "design-sgrna",
        help="Design sgRNA candidates from a target DNA sequence.",
    )
    design_sgrna_cmd.add_argument(
        "--input", required=True, type=Path,
        help="Path to a FASTA or plain text file containing the target sequence.",
    )
    design_sgrna_cmd.add_argument("--output", "-o", type=Path, help="Write results to this JSON file.")
    design_sgrna_cmd.add_argument(
        "--nuclease", default="SpCas9",
        choices=[n.value for n in NucleaseType],
        help="CRISPR nuclease type (default: SpCas9).",
    )
    design_sgrna_cmd.add_argument("--guide-length", type=int, default=20)
    design_sgrna_cmd.add_argument("--max-candidates", type=int, default=50)

    # offtarget: Search for off-target sites
    offtarget_cmd = subparsers.add_parser(
        "offtarget",
        help="Search for off-target sites for a guide RNA.",
    )
    offtarget_cmd.add_argument(
        "--guide", required=True,
        help="20-nt guide RNA sequence.",
    )
    offtarget_cmd.add_argument(
        "--reference", type=Path,
        help="FASTA file of reference sequences to search.",
    )
    offtarget_cmd.add_argument(
        "--nuclease", default="SpCas9",
        choices=[n.value for n in NucleaseType],
    )
    offtarget_cmd.add_argument("--max-mismatches", type=int, default=4)
    offtarget_cmd.add_argument("--output", "-o", type=Path)

    # predict: Predict editing efficiency for a guide
    predict_cmd = subparsers.add_parser(
        "predict",
        help="Predict editing efficiency for a guide RNA.",
    )
    predict_cmd.add_argument("--guide", required=True, help="20-nt guide RNA sequence.")
    predict_cmd.add_argument(
        "--species", default="mouse",
        choices=sorted(PROFILES),
        help="Target species (default: mouse).",
    )
    predict_cmd.add_argument(
        "--context-35nt",
        help="35-nt context sequence (for CRISPRscan zebrafish scoring).",
    )
    predict_cmd.add_argument("--output", "-o", type=Path)

    # pipeline: Run the full end-to-end pipeline
    pipeline_cmd = subparsers.add_parser(
        "pipeline",
        help="Run the full prediction pipeline: design -> predict -> off-target -> assess.",
    )
    pipeline_cmd.add_argument(
        "--input", required=True, type=Path,
        help="FASTA file containing the target sequence.",
    )
    pipeline_cmd.add_argument("--chrom", help="Chromosome/sequence ID (default: first sequence).")
    pipeline_cmd.add_argument("--start", type=int, help="1-based start position (default: 1).")
    pipeline_cmd.add_argument("--end", type=int, help="1-based end position (default: full sequence).")
    pipeline_cmd.add_argument(
        "--reference", type=Path,
        help="Reference genome FASTA for off-target search (default: use input file).",
    )
    pipeline_cmd.add_argument(
        "--species", default="mouse", choices=sorted(PROFILES),
    )
    pipeline_cmd.add_argument(
        "--nuclease", default="SpCas9",
        choices=[n.value for n in NucleaseType],
    )
    pipeline_cmd.add_argument("--strain", default="C57BL/6J")
    pipeline_cmd.add_argument("--genome-build", default="GRCm39")
    pipeline_cmd.add_argument("--edit-class", default="knockout")
    pipeline_cmd.add_argument("--gene-essentiality", type=float, default=0.0)
    pipeline_cmd.add_argument("--phenotype-severity", type=float, default=0.0)
    pipeline_cmd.add_argument("--top-k", type=int, default=10)
    pipeline_cmd.add_argument("--max-candidates", type=int, default=50)
    pipeline_cmd.add_argument("--max-mismatches", type=int, default=4)
    pipeline_cmd.add_argument("--min-efficiency", type=float, default=0.3)
    pipeline_cmd.add_argument("--min-specificity", type=float, default=0.5)
    pipeline_cmd.add_argument("--output", "-o", required=True, type=Path)
    pipeline_cmd.add_argument(
        "--html", type=Path,
        help="Also generate an HTML visualization report.",
    )

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
    if args.command == "import-indelphi":
        try:
            input_bytes = args.input.read_bytes()
            document = json.loads(input_bytes.decode("utf-8"))
            if not isinstance(document, dict):
                raise ValueError("inDelphi input must be a JSON object.")
            prediction = normalize_indelphi(
                document,
                source_document_sha256=sha256(input_bytes).hexdigest(),
            )
        except (
            OSError,
            UnicodeDecodeError,
            json.JSONDecodeError,
            TypeError,
            ValueError,
        ) as error:
            parser.error(str(error))
        rendered = json.dumps(
            asdict(prediction),
            indent=2,
            ensure_ascii=False,
        ) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"inDelphi audit record written to {args.output}")
        else:
            print(rendered, end="")
        return
    if args.command == "import-housden":
        try:
            input_bytes = args.input.read_bytes()
            document = json.loads(input_bytes.decode("utf-8"))
            if not isinstance(document, dict):
                raise ValueError("Housden input must be a JSON object.")
            prediction = normalize_housden(
                document,
                source_response=args.source_response.read_bytes(),
                source_document_sha256=sha256(input_bytes).hexdigest(),
            )
        except (
            OSError,
            UnicodeDecodeError,
            json.JSONDecodeError,
            TypeError,
            ValueError,
        ) as error:
            parser.error(str(error))
        rendered = json.dumps(
            asdict(prediction),
            indent=2,
            ensure_ascii=False,
        ) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"Housden audit record written to {args.output}")
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
    if args.command == "readiness":
        if args.all:
            report = {
                key: asdict(value)
                for key, value in readiness_matrix().items()
            }
        else:
            report = asdict(readiness_for_species(args.species))
        print(json.dumps(report, indent=2, ensure_ascii=False))
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
    if args.command == "prepare-rat-guide-transfer":
        try:
            template = prepare_rat_guide_transfer_template(
                args.table1,
                args.table5,
            )
        except (OSError, TypeError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(template, indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"Rat guide-transfer prediction template written to {args.output}")
        return
    if args.command == "validate-rat-guide-transfer":
        try:
            predictions = json.loads(
                args.predictions.read_text(encoding="utf-8")
            )
            if not isinstance(predictions, dict):
                raise ValueError("rat guide-transfer predictions must be a JSON object.")
            report = evaluate_rat_guide_transfer(
                args.table1,
                args.table5,
                predictions,
            )
        except (
            OSError,
            json.JSONDecodeError,
            TypeError,
            ValueError,
        ) as error:
            parser.error(str(error))
        rendered = json.dumps(asdict(report), indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"Rat guide-transfer report written to {args.output}")
        return
    if args.command == "audit-fruit-fly-cas12a-evidence":
        try:
            audit = audit_fruit_fly_cas12a_evidence(
                args.library,
                args.genotypes,
                args.source_data,
                source=PORT_2026_CAS12A_SOURCE,
            )
            report = {"audit": asdict(audit)}
            if args.line_id:
                report["array_evidence"] = asdict(
                    lookup_fruit_fly_cas12a_array(
                        args.library,
                        args.genotypes,
                        args.source_data,
                        args.line_id,
                        source=PORT_2026_CAS12A_SOURCE,
                    )
                )
        except (OSError, TypeError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"Fruit-fly Cas12a evidence audit written to {args.output}")
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

    # --- New prediction pipeline command handlers ---

    if args.command == "design-sgrna":
        try:
            raw = args.input.read_text(encoding="utf-8").strip()
            # Strip FASTA header if present
            if raw.startswith(">"):
                lines = raw.split("\n")
                sequence = "".join(lines[1:]).replace(" ", "").replace("\n", "").upper()
            else:
                sequence = raw.replace(" ", "").replace("\n", "").upper()
            nuclease = NucleaseType(args.nuclease)
            result = design_sgrnas(
                sequence=sequence,
                chrom="target",
                nuclease=nuclease,
                guide_length=args.guide_length,
                max_candidates=args.max_candidates,
            )
            output = {
                "target_id": result.target_id,
                "nuclease": result.nuclease.value,
                "candidate_count": result.count,
                "warnings": result.warnings,
                "candidates": [
                    {
                        "guide_id": c.guide_id,
                        "guide_sequence": c.guide_sequence,
                        "pam": c.pam,
                        "strand": c.strand,
                        "start": c.start,
                        "end": c.end,
                        "gc_content": round(c.gc_content, 4),
                        "features": {k: round(v, 4) for k, v in c.features.items()},
                    }
                    for c in result.top_candidates(20)
                ],
            }
        except (OSError, ValueError) as error:
            parser.error(str(error))
        rendered = json.dumps(output, indent=2, ensure_ascii=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"sgRNA design results written to {args.output}")
        else:
            print(rendered, end="")
        return

    if args.command == "offtarget":
        try:
            guide = args.guide.upper().replace(" ", "")
            if len(guide) != 20:
                raise ValueError(f"guide must be 20 nt, got {len(guide)}")
            ref_seqs = None
            if args.reference:
                reader = FastaReader(args.reference)
                ref_seqs = {
                    sid: reader[sid].sequence
                    for sid in reader.sequence_ids[:10]
                }
            nuclease = NucleaseType(args.nuclease)
            report = find_offtargets(
                guide_sequence=guide,
                reference_sequences=ref_seqs,
                nuclease=nuclease,
                max_mismatches=args.max_mismatches,
            )
            output = {
                "guide_sequence": report.guide_sequence,
                "nuclease": report.nuclease.value,
                "total_sites_scanned": report.total_sites_scanned,
                "high_risk_count": report.high_risk_count,
                "moderate_risk_count": report.moderate_risk_count,
                "low_risk_count": report.low_risk_count,
                "specificity_score": round(report.specificity_score, 4),
                "warnings": report.warnings,
                "off_targets": [
                    {
                        "chrom": ot.chrom,
                        "start": ot.start,
                        "end": ot.end,
                        "strand": ot.strand,
                        "sequence": ot.off_target_sequence,
                        "pam": ot.pam,
                        "mismatches": ot.mismatch_count,
                        "mismatch_positions": list(ot.mismatch_positions),
                        "score": round(ot.score, 4),
                        "risk_level": ot.risk_level,
                    }
                    for ot in report.off_targets[:50]
                ],
            }
        except (OSError, ValueError, KeyError) as error:
            parser.error(str(error))
        rendered = json.dumps(output, indent=2, ensure_ascii=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"Off-target report written to {args.output}")
        else:
            print(rendered, end="")
        return

    if args.command == "predict":
        try:
            guide = args.guide.upper().replace(" ", "")
            if len(guide) != 20:
                raise ValueError(f"guide must be 20 nt, got {len(guide)}")
            from .sgrna_design import SgrnaCandidate, compute_guide_features
            features = compute_guide_features(guide)
            candidate = SgrnaCandidate(
                guide_id="input",
                guide_sequence=guide,
                pam="NGG",
                pam_strand="+",
                chrom="N/A",
                start=0,
                end=0,
                strand="+",
                context_30nt="",
                context_35nt=args.context_35nt.upper() if args.context_35nt else "",
                gc_content=features["gc_content"],
                nuclease=NucleaseType.SPCAS9,
                features=features,
            )
            eff = predict_efficiency(candidate, species_key=args.species)
            indel = predict_indel_outcomes(guide, species_key=args.species)
            output = {
                "guide_sequence": eff.guide_sequence,
                "species": eff.species,
                "nuclease": eff.nuclease,
                "efficiency": {
                    "score": round(eff.efficiency_score, 4),
                    "confidence": round(eff.confidence, 4),
                    "model": eff.model_name,
                    "model_version": eff.model_version,
                    "warnings": list(eff.warnings),
                },
                "indel_outcome": {
                    "insertion_rate": round(indel.insertion_rate, 4),
                    "deletion_rate": round(indel.deletion_rate, 4),
                    "no_edit_rate": round(indel.no_edit_rate, 4),
                    "most_likely_outcome": indel.most_likely_outcome,
                    "predicted_indel_size": indel.predicted_indel_size,
                },
                "features": {k: round(v, 4) for k, v in eff.features.items()},
            }
        except (ValueError,) as error:
            parser.error(str(error))
        rendered = json.dumps(output, indent=2, ensure_ascii=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
            print(f"Efficiency prediction written to {args.output}")
        else:
            print(rendered, end="")
        return

    if args.command == "pipeline":
        try:
            reader = FastaReader(args.input)
            if args.chrom:
                chrom = args.chrom
            else:
                chrom = reader.sequence_ids[0]
            seq_len = reader[chrom].length
            start = args.start or 1
            end = args.end or seq_len
            target_seq = reader.fetch(chrom, start, end)

            # Build reference sequences for off-target search
            ref_seqs = None
            if args.reference:
                ref_reader = FastaReader(args.reference)
                ref_seqs = {
                    sid: ref_reader[sid].sequence
                    for sid in ref_reader.sequence_ids[:10]
                }
            else:
                # Use the input file as reference
                ref_seqs = {
                    sid: reader[sid].sequence
                    for sid in reader.sequence_ids[:10]
                }

            species_profile = PROFILES[args.species]
            study_context = StudyContext(
                species=args.species,
                strain_or_breed=args.strain,
                genome_build=args.genome_build,
                edit_class=args.edit_class,
                evidence_snapshot=f"pipeline_run_{args.species}_{args.nuclease}",
            )

            config = PipelineConfig(
                species=args.species,
                nuclease=NucleaseType(args.nuclease),
                max_candidates=args.max_candidates,
                max_offtargets=args.max_mismatches,
                gene_essentiality=args.gene_essentiality,
                phenotype_severity=args.phenotype_severity,
                min_efficiency=args.min_efficiency,
                min_specificity=args.min_specificity,
                top_k=args.top_k,
            )

            from .pipeline import run_pipeline as _run
            report = _run(
                sequence=target_seq,
                config=config,
                study_context=study_context,
                reference_sequences=ref_seqs,
            )
        except (OSError, ValueError, KeyError) as error:
            parser.error(str(error))

        json_text = report.to_json(args.output)
        print(f"Pipeline report written to {args.output}")

        if args.html:
            html_text = generate_html_report(report)
            args.html.parent.mkdir(parents=True, exist_ok=True)
            args.html.write_text(html_text, encoding="utf-8")
            print(f"HTML visualization written to {args.html}")
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
