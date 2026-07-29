"""Command-line interface for auditable assessment reports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .workflow import assess_request


def main() -> None:
    parser = argparse.ArgumentParser(description="Animal genome-edit impact assessment (research use only).")
    subparsers = parser.add_subparsers(dest="command", required=True)
    assess = subparsers.add_parser("assess", help="Create an assessment report from a JSON request.")
    assess.add_argument("request", type=Path, help="Path to an assessment request JSON file.")
    assess.add_argument("--output", "-o", type=Path, help="Write the report to this JSON file.")
    assess.add_argument("--model-version", default="0.2.0")
    args = parser.parse_args()

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
