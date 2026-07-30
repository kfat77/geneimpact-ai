"""Flask web application for interactive gene editing prediction.

Provides a REST API and interactive frontend for:
- sgRNA design from genomic sequences
- Editing efficiency prediction
- Off-target detection
- Full end-to-end pipeline execution
- Genome sequence download
- Result visualization and export

Run::

    pip install flask
    python -m geneimpact.webapp
    # Open http://localhost:5000

API Endpoints:
    POST /api/design-sgrna      — Design sgRNA candidates
    POST /api/predict            — Predict editing efficiency
    POST /api/offtarget          — Off-target search
    POST /api/pipeline           — Full pipeline run
    POST /api/download-genome    — Download genome sequence
    GET  /api/species            — List supported species
    GET  /api/health             — Health check
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from flask import Flask, jsonify, request, send_from_directory

from .genomics import FastaReader
from .sgrna_design import NucleaseType, design_sgrnas
from .offtarget import find_offtargets
from .fast_offtarget import fast_find_offtargets, build_seed_index, SeedIndex
from .efficiency import predict_efficiency, predict_indel_outcomes, compute_evidence_scores
from .pipeline import PipelineConfig, run_pipeline
from .provenance import StudyContext
from .visualization import generate_html_report
from .species import PROFILES
from .advanced_models import score_ruleset2, MODEL_INFO
from .genome_downloader import download_sequence, list_species as list_download_species

__all__ = ["create_app", "run_app"]

app = Flask(__name__, static_folder=None)


# ---------------------------------------------------------------------------
# Frontend serving
# ---------------------------------------------------------------------------

_FRONTEND_DIR = Path(__file__).parent / "webapp_static"


@app.route("/")
def index():
    """Serve the interactive frontend."""
    index_path = _FRONTEND_DIR / "index.html"
    if index_path.exists():
        return send_from_directory(str(_FRONTEND_DIR), "index.html")
    return jsonify({"error": "Frontend not built. See webapp/README.md"}), 404


@app.route("/<path:filename>")
def static_files(filename: str):
    """Serve static assets."""
    return send_from_directory(str(_FRONTEND_DIR), filename)


# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------

@app.route("/api/health")
def health():
    """Health check endpoint."""
    return jsonify({
        "status": "ok",
        "version": "1.1.0",
        "model": MODEL_INFO["name"],
        "model_version": MODEL_INFO["version"],
    })


@app.route("/api/species")
def species_list():
    """List supported species."""
    species_data = []
    for key, profile in sorted(PROFILES.items()):
        species_data.append({
            "key": key,
            "common_name": profile.common_name,
            "scientific_name": profile.scientific_name,
            "assembly": getattr(profile, "assembly_name", ""),
        })
    downloadable = list_download_species()
    return jsonify({
        "supported_species": species_data,
        "downloadable_species": downloadable,
        "nucleases": [n.value for n in NucleaseType],
    })


@app.route("/api/design-sgrna", methods=["POST"])
def api_design_sgrna():
    """Design sgRNA candidates from a target sequence."""
    data = request.get_json()
    if not data or "sequence" not in data:
        return jsonify({"error": "Missing 'sequence' field"}), 400

    sequence = data["sequence"].upper().replace(" ", "").replace("\n", "")
    nuclease_str = data.get("nuclease", "SpCas9")
    max_candidates = data.get("max_candidates", 50)
    guide_length = data.get("guide_length", 20)

    try:
        nuclease = NucleaseType(nuclease_str)
    except ValueError:
        return jsonify({"error": f"Unknown nuclease: {nuclease_str}"}), 400

    try:
        result = design_sgrnas(
            sequence=sequence,
            chrom=data.get("chrom", "target"),
            nuclease=nuclease,
            guide_length=guide_length,
            max_candidates=max_candidates,
        )
    except (ValueError, TypeError) as e:
        return jsonify({"error": str(e)}), 400

    return jsonify({
        "count": result.count,
        "nuclease": nuclease.value,
        "candidates": [
            {
                "guide_id": c.guide_id,
                "guide_sequence": c.guide_sequence,
                "pam": c.pam,
                "strand": c.strand,
                "chrom": c.chrom,
                "start": c.start,
                "end": c.end,
                "gc_content": round(c.gc_content, 4),
                "features": {k: round(v, 4) if isinstance(v, float) else v
                             for k, v in c.features.items()},
            }
            for c in result.candidates
        ],
        "warnings": result.warnings,
    })


@app.route("/api/predict", methods=["POST"])
def api_predict():
    """Predict editing efficiency for a guide RNA."""
    data = request.get_json()
    if not data or "guide" not in data:
        return jsonify({"error": "Missing 'guide' field"}), 400

    guide = data["guide"].upper().replace(" ", "")
    species = data.get("species", "mouse")
    nuclease_str = data.get("nuclease", "SpCas9")

    try:
        nuclease = NucleaseType(nuclease_str)
    except ValueError:
        return jsonify({"error": f"Unknown nuclease: {nuclease_str}"}), 400

    from .sgrna_design import SgrnaCandidate, compute_guide_features

    features = compute_guide_features(guide)
    candidate = SgrnaCandidate(
        guide_id="api_input",
        guide_sequence=guide,
        pam=data.get("pam", "NGG"),
        pam_strand="+",
        chrom="N/A",
        start=0,
        end=0,
        strand="+",
        context_30nt="",
        context_35nt=data.get("context_35nt", "").upper(),
        gc_content=features.get("gc_content", 0.5),
        nuclease=nuclease,
        features=features,
    )

    try:
        efficiency = predict_efficiency(candidate, species_key=species)
        indel = predict_indel_outcomes(guide, species_key=species)

        # Also compute Rule Set 2 detailed score for non-zebrafish
        rs2_detail = None
        if species != "zebrafish":
            rs2 = score_ruleset2(guide, species)
            rs2_detail = {
                "raw_score": round(rs2.raw_score, 4),
                "calibrated_score": round(rs2.calibrated_score, 4),
                "confidence": round(rs2.confidence, 4),
                "pwm_contribution": round(rs2.pwm_contribution, 4),
                "thermo_contribution": round(rs2.thermo_contribution, 4),
                "composition_contribution": round(rs2.composition_contribution, 4),
                "features": {k: round(v, 4) if isinstance(v, float) else v
                             for k, v in rs2.features.items() if not k.startswith("pos_")},
            }
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    result = {
        "guide": guide,
        "species": species,
        "nuclease": nuclease.value,
        "efficiency": {
            "score": round(efficiency.efficiency_score, 4),
            "confidence": round(efficiency.confidence, 4),
            "model": efficiency.model_name,
            "model_version": efficiency.model_version,
            "warnings": list(efficiency.warnings),
        },
        "indel_outcome": {
            "insertion_rate": round(indel.insertion_rate, 4),
            "deletion_rate": round(indel.deletion_rate, 4),
            "no_edit_rate": round(indel.no_edit_rate, 4),
            "most_likely": indel.most_likely_outcome,
            "predicted_size": indel.predicted_indel_size,
            "confidence": round(indel.confidence, 4),
        },
    }
    if rs2_detail:
        result["ruleset2_detail"] = rs2_detail

    return jsonify(result)


@app.route("/api/offtarget", methods=["POST"])
def api_offtarget():
    """Search for off-target sites."""
    data = request.get_json()
    if not data or "guide" not in data:
        return jsonify({"error": "Missing 'guide' field"}), 400

    guide = data["guide"].upper().replace(" ", "")
    nuclease_str = data.get("nuclease", "SpCas9")
    max_mismatches = data.get("max_mismatches", 4)
    use_fast = data.get("use_fast", True)
    reference = data.get("reference", {})

    try:
        nuclease = NucleaseType(nuclease_str)
    except ValueError:
        return jsonify({"error": f"Unknown nuclease: {nuclease_str}"}), 400

    if not reference:
        return jsonify({"error": "Missing 'reference' field (dict of {seq_id: sequence})"}), 400

    try:
        if use_fast:
            report = fast_find_offtargets(
                guide_sequence=guide,
                reference_sequences=reference,
                nuclease=nuclease,
                max_mismatches=max_mismatches,
            )
        else:
            report = find_offtargets(
                guide_sequence=guide,
                reference_sequences=reference,
                nuclease=nuclease,
                max_mismatches=max_mismatches,
            )
    except (ValueError, TypeError) as e:
        return jsonify({"error": str(e)}), 400

    return jsonify({
        "guide": guide,
        "nuclease": nuclease.value,
        "total_sites_scanned": report.total_sites_scanned,
        "high_risk_count": report.high_risk_count,
        "moderate_risk_count": report.moderate_risk_count,
        "low_risk_count": report.low_risk_count,
        "specificity_score": round(report.specificity_score, 4),
        "off_targets": [
            {
                "chrom": ot.chrom,
                "start": ot.start,
                "end": ot.end,
                "strand": ot.strand,
                "sequence": ot.off_target_sequence,
                "pam": ot.pam,
                "mismatches": ot.mismatch_count,
                "score": round(ot.score, 4),
                "risk": ot.risk_level,
                "mismatch_positions": list(ot.mismatch_positions),
            }
            for ot in report.off_targets[:50]
        ],
        "warnings": report.warnings,
        "algorithm": "seed-and-extend" if use_fast else "brute-force",
    })


@app.route("/api/pipeline", methods=["POST"])
def api_pipeline():
    """Run the full prediction pipeline."""
    data = request.get_json()
    if not data or "sequence" not in data:
        return jsonify({"error": "Missing 'sequence' field"}), 400

    sequence = data["sequence"].upper().replace(" ", "").replace("\n", "")
    species = data.get("species", "mouse")
    nuclease_str = data.get("nuclease", "SpCas9")
    reference = data.get("reference", {})

    try:
        nuclease = NucleaseType(nuclease_str)
    except ValueError:
        return jsonify({"error": f"Unknown nuclease: {nuclease_str}"}), 400

    # Validate species
    if species not in PROFILES:
        return jsonify({"error": f"Unknown species: {species}"}), 400

    profile = PROFILES[species]
    study_context = StudyContext(
        species=species,
        strain=data.get("strain", getattr(profile, "default_strain", "")),
        genome_build=data.get("genome_build", getattr(profile, "assembly_name", "")),
        edit_class=data.get("edit_class", "knockout"),
        nuclease=nuclease.value,
    )

    config = PipelineConfig(
        species=species,
        nuclease=nuclease,
        max_candidates=data.get("max_candidates", 50),
        max_offtargets=data.get("max_mismatches", 4),
        top_k=data.get("top_k", 10),
        min_efficiency=data.get("min_efficiency", 0.3),
        min_specificity=data.get("min_specificity", 0.5),
        gene_essentiality=data.get("gene_essentiality", 0.0),
        phenotype_severity=data.get("phenotype_severity", 0.0),
    )

    try:
        report = run_pipeline(
            sequence=sequence,
            config=config,
            study_context=study_context,
            reference_sequences=reference if reference else None,
        )
    except (ValueError, TypeError, KeyError) as e:
        return jsonify({"error": str(e)}), 400

    return jsonify(report.to_dict())


@app.route("/api/download-genome", methods=["POST"])
def api_download_genome():
    """Download a genome sequence from Ensembl/NCBI."""
    data = request.get_json()
    if not data or "species" not in data:
        return jsonify({"error": "Missing 'species' field"}), 400

    species = data["species"]
    chrom = data.get("chrom", "1")
    start = data.get("start")
    end = data.get("end")
    source = data.get("source", "ensembl")
    cache_dir = data.get("cache_dir", "./genome_cache")

    try:
        result = download_sequence(
            species=species,
            chrom=chrom,
            start=start,
            end=end,
            cache_dir=cache_dir,
            source=source,
        )
    except (ValueError, OSError) as e:
        return jsonify({"error": str(e)}), 400

    # Read the downloaded sequence for response
    seq = ""
    if result.local_path and Path(result.local_path).exists():
        content = Path(result.local_path).read_text(encoding="utf-8")
        lines = content.strip().split("\n")
        seq = "".join(l for l in lines if not l.startswith(">"))

    return jsonify({
        "species": result.species,
        "chrom": result.chrom,
        "source": result.source,
        "assembly": result.assembly,
        "sequence_length": result.sequence_length,
        "sha256": result.sha256,
        "start": result.start,
        "end": result.end,
        "cached": result.cached,
        "warnings": list(result.warnings),
        "sequence_preview": seq[:200] + "..." if len(seq) > 200 else seq,
    })


# ---------------------------------------------------------------------------
# Error handlers
# ---------------------------------------------------------------------------

@app.errorhandler(404)
def not_found(e):
    return jsonify({"error": "Not found"}), 404


@app.errorhandler(500)
def server_error(e):
    return jsonify({"error": "Internal server error", "detail": str(e)}), 500


def create_app() -> Flask:
    """Factory function for creating the Flask app."""
    return app


def run_app(host: str = "0.0.0.0", port: int = 5000, debug: bool = False) -> None:
    """Run the web application server."""
    print(f"\n  GeneImpact AI Web Application")
    print(f"  =============================")
    print(f"  Server: http://{host}:{port}")
    print(f"  API docs: http://{host}:{port}/api/health")
    print(f"  Model: {MODEL_INFO['name']} v{MODEL_INFO['version']}")
    print(f"  Features: {MODEL_INFO['features']} feature model")
    print()
    app.run(host=host, port=port, debug=debug)


if __name__ == "__main__":
    run_app()
