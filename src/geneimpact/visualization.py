"""Result visualization: self-contained HTML/SVG reports and charts.

Generates interactive HTML reports with embedded SVG visualizations
for guide RNA positions, off-target heatmaps, efficiency distributions,
and evidence radar charts. No external dependencies required for
core HTML output; matplotlib is optional for static PNG export.
"""

from __future__ import annotations

import html
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Sequence

from .pipeline import GuideResult, PipelineReport

__all__ = [
    "generate_html_report",
    "generate_sgrna_diagram",
    "generate_offtarget_heatmap",
    "generate_efficiency_chart",
    "generate_evidence_radar",
]


def generate_html_report(report: PipelineReport) -> str:
    """Generate a complete self-contained HTML report.

    The report includes:
    - Study context and species validation
    - Guide ranking table with sortable columns
    - SVG diagrams for each top guide
    - Off-target heatmap
    - Efficiency distribution chart
    - Evidence radar chart
    """
    data = report.to_dict()
    guides = data["guides"]

    # Build guide table rows
    table_rows = []
    for g in guides:
        eff_score = g["efficiency"]["score"]
        specificity = g.get("offtarget", {}).get("specificity_score", 0.5)
        tier = g.get("assessment", {}).get("tier", "N/A")
        concern = g.get("assessment", {}).get("concern_score", 0)

        tier_color = {
            "high_concern_review": "#e74c3c",
            "enhanced_review": "#f39c12",
            "standard_review": "#27ae60",
        }.get(tier, "#95a5a6")

        table_rows.append(f"""
        <tr>
            <td>{g['rank']}</td>
            <td class="mono">{html.escape(g['guide_sequence'])}</td>
            <td class="mono">{html.escape(g['pam'])}</td>
            <td>{g['strand']}</td>
            <td>{g['chrom']}:{g['start']}-{g['end']}</td>
            <td>{g['gc_content']:.1%}</td>
            <td><span class="score-bar"><span class="score-fill" style="width:{eff_score*100:.0f}%;background:#3498db"></span></span>{eff_score:.3f}</td>
            <td><span class="score-bar"><span class="score-fill" style="width:{specificity*100:.0f}%;background:#2ecc71"></span></span>{specificity:.3f}</td>
            <td><span class="tier-badge" style="background:{tier_color}">{tier.replace('_review','').replace('_',' ')}</span></td>
            <td>{g['recommendation']}</td>
        </tr>""")

    # Build off-target data for heatmap
    offtarget_data = []
    for g in guides[:5]:  # Top 5 guides
        ot = g.get("offtarget", {})
        offtarget_data.append({
            "guide_id": g["guide_id"],
            "guide_seq": g["guide_sequence"],
            "high_risk": ot.get("high_risk", 0),
            "moderate_risk": ot.get("moderate_risk", 0),
            "low_risk": ot.get("low_risk", 0),
            "specificity": ot.get("specificity_score", 1.0),
        })

    # Build evidence radar data for top guide
    radar_data = {}
    if guides:
        ev = guides[0].get("evidence_scores", {})
        radar_data = {
            "On-target Uncertainty": ev.get("on_target_uncertainty", 0),
            "Off-target Evidence": ev.get("off_target_evidence", 0),
            "Network Impact": ev.get("network_impact_evidence", 0),
            "Welfare Relevance": ev.get("welfare_relevance", 0),
        }

    # Build efficiency distribution
    eff_values = [g["efficiency"]["score"] for g in guides]

    warnings_html = ""
    if data["warnings"]:
        warnings_html = '<div class="warnings"><h3>Warnings</h3><ul>' + \
            "".join(f"<li>{html.escape(w)}</li>" for w in data["warnings"]) + \
            "</ul></div>"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GeneImpact AI - Prediction Report</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       background: #f5f7fa; color: #2c3e50; line-height: 1.6; padding: 20px; }}
.container {{ max-width: 1200px; margin: 0 auto; }}
.header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
           color: white; padding: 30px; border-radius: 12px; margin-bottom: 24px; }}
.header h1 {{ font-size: 28px; margin-bottom: 8px; }}
.header .meta {{ opacity: 0.9; font-size: 14px; }}
.card {{ background: white; border-radius: 12px; padding: 24px; margin-bottom: 20px;
         box-shadow: 0 2px 8px rgba(0,0,0,0.06); }}
.card h2 {{ font-size: 20px; margin-bottom: 16px; color: #2c3e50;
            border-bottom: 2px solid #ecf0f1; padding-bottom: 8px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th {{ text-align: left; padding: 10px 8px; background: #f8f9fa; font-weight: 600;
      border-bottom: 2px solid #dee2e6; white-space: nowrap; }}
td {{ padding: 8px; border-bottom: 1px solid #ecf0f1; vertical-align: middle; }}
tr:hover {{ background: #f8f9fa; }}
.mono {{ font-family: 'SF Mono', Monaco, Consolas, monospace; font-size: 12px; }}
.score-bar {{ display: inline-block; width: 60px; height: 8px; background: #ecf0f1;
              border-radius: 4px; margin-right: 6px; overflow: hidden; vertical-align: middle; }}
.score-fill {{ display: block; height: 100%; border-radius: 4px; }}
.tier-badge {{ display: inline-block; padding: 3px 10px; border-radius: 12px;
               color: white; font-size: 11px; font-weight: 600; text-transform: uppercase; }}
.warnings {{ background: #fff3cd; border: 1px solid #ffeaa7; border-radius: 8px;
             padding: 16px; margin-bottom: 20px; }}
.warnings h3 {{ color: #856404; margin-bottom: 8px; }}
.warnings ul {{ padding-left: 20px; }}
.warnings li {{ color: #856404; margin-bottom: 4px; }}
.notice {{ background: #e8f4fd; border-left: 4px solid #3498db; padding: 12px 16px;
           border-radius: 4px; margin-bottom: 20px; font-size: 13px; color: #2980b9; }}
.grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
@media (max-width: 768px) {{ .grid-2 {{ grid-template-columns: 1fr; }} }}
.summary-stats {{ display: flex; gap: 16px; flex-wrap: wrap; }}
.stat {{ background: #f8f9fa; padding: 12px 20px; border-radius: 8px; text-align: center; }}
.stat .value {{ font-size: 24px; font-weight: 700; color: #2c3e50; }}
.stat .label {{ font-size: 12px; color: #7f8c8d; text-transform: uppercase; }}
</style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>Gene Editing Prediction Report</h1>
        <div class="meta">
            Species: {html.escape(data['study_context']['species'])} |
            Nuclease: {data['config']['nuclease']} |
            Generated: {data['timestamp'][:19]} |
            Pipeline v{data['pipeline_version']}
        </div>
    </div>

    <div class="notice">{html.escape(data['report_notice'])}</div>

    {warnings_html}

    <div class="card">
        <h2>Summary</h2>
        <div class="summary-stats">
            <div class="stat"><div class="value">{len(guides)}</div><div class="label">Guides Analyzed</div></div>
            <div class="stat"><div class="value">{sum(1 for g in guides if g.get('assessment',{}).get('tier')=='standard_review')}</div><div class="label">Standard Review</div></div>
            <div class="stat"><div class="value">{sum(1 for g in guides if g.get('assessment',{}).get('tier')=='enhanced_review')}</div><div class="label">Enhanced Review</div></div>
            <div class="stat"><div class="value">{sum(1 for g in guides if g.get('assessment',{}).get('tier')=='high_concern_review')}</div><div class="label">High Concern</div></div>
        </div>
    </div>

    <div class="card">
        <h2>Guide RNA Rankings</h2>
        <div style="overflow-x:auto;">
        <table>
            <thead>
                <tr>
                    <th>#</th>
                    <th>Guide Sequence (20nt)</th>
                    <th>PAM</th>
                    <th>Strand</th>
                    <th>Location</th>
                    <th>GC%</th>
                    <th>Efficiency</th>
                    <th>Specificity</th>
                    <th>Review Tier</th>
                    <th>Recommendation</th>
                </tr>
            </thead>
            <tbody>
                {''.join(table_rows)}
            </tbody>
        </table>
        </div>
    </div>

    <div class="grid-2">
        <div class="card">
            <h2>Evidence Profile (Top Guide)</h2>
            {_render_radar_svg(radar_data)}
        </div>
        <div class="card">
            <h2>Efficiency Distribution</h2>
            {_render_efficiency_svg(eff_values)}
        </div>
    </div>

    <div class="card">
        <h2>Off-Target Analysis (Top 5 Guides)</h2>
        {_render_offtarget_svg(offtarget_data)}
    </div>

    <div class="card">
        <h2>Study Context</h2>
        <table>
            <tr><td><strong>Species</strong></td><td>{html.escape(data['study_context']['species'])}</td></tr>
            <tr><td><strong>Strain/Breed</strong></td><td>{html.escape(data['study_context']['strain_or_breed'])}</td></tr>
            <tr><td><strong>Genome Build</strong></td><td>{html.escape(data['study_context']['genome_build'])}</td></tr>
            <tr><td><strong>Edit Class</strong></td><td>{html.escape(data['study_context']['edit_class'])}</td></tr>
            <tr><td><strong>Evidence Snapshot</strong></td><td>{html.escape(data['study_context']['evidence_snapshot'])}</td></tr>
            <tr><td><strong>Species Validated</strong></td><td>{'Yes' if data['species_validation']['supported'] else 'No'}</td></tr>
        </table>
    </div>
</div>
</body>
</html>"""
    return html_content


def _render_radar_svg(data: dict[str, float]) -> str:
    """Render an SVG radar chart for evidence scores."""
    if not data:
        return "<p>No data available</p>"

    labels = list(data.keys())
    values = list(data.values())
    n = len(labels)

    cx, cy, r = 150, 150, 100
    angles = [(-90 + i * 360 / n) * 3.14159 / 180 for i in range(n)]

    # Grid circles
    grid_circles = ""
    for radius_pct in [0.25, 0.5, 0.75, 1.0]:
        r_val = r * radius_pct
        grid_circles += f'<circle cx="{cx}" cy="{cy}" r="{r_val}" fill="none" stroke="#ecf0f1" stroke-width="1"/>'

    # Axes and labels
    axes = ""
    labels_svg = ""
    for i, (label, angle) in enumerate(zip(labels, angles)):
        x = cx + r * __import__("math").cos(angle)
        y = cy + r * __import__("math").sin(angle)
        axes += f'<line x1="{cx}" y1="{cy}" x2="{x:.1f}" y2="{y:.1f}" stroke="#bdc3c7" stroke-width="1"/>'
        lx = cx + (r + 25) * __import__("math").cos(angle)
        ly = cy + (r + 25) * __import__("math").sin(angle)
        labels_svg += f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" alignment-baseline="middle" font-size="10" fill="#7f8c8d">{label}</text>'

    # Data polygon
    points = []
    for i, (value, angle) in enumerate(zip(values, angles)):
        x = cx + r * value * __import__("math").cos(angle)
        y = cy + r * value * __import__("math").sin(angle)
        points.append(f"{x:.1f},{y:.1f}")
    polygon = f'<polygon points="{" ".join(points)}" fill="rgba(102,126,234,0.3)" stroke="#667eea" stroke-width="2"/>'

    # Value labels
    value_labels = ""
    for i, (value, angle) in enumerate(zip(values, angles)):
        x = cx + r * value * 0.7 * __import__("math").cos(angle)
        y = cy + r * value * 0.7 * __import__("math").sin(angle)
        value_labels += f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" font-size="9" fill="#2c3e50" font-weight="bold">{value:.2f}</text>'

    return f"""<svg viewBox="0 0 300 300" style="width:100%;max-width:350px;margin:auto;display:block">
        {grid_circles}
        {axes}
        {polygon}
        {value_labels}
        {labels_svg}
    </svg>"""


def _render_efficiency_svg(values: list[float]) -> str:
    """Render an SVG bar chart for efficiency scores."""
    if not values:
        return "<p>No data available</p>"

    bar_width = 30
    gap = 8
    chart_height = 200
    chart_width = max(300, len(values) * (bar_width + gap) + 40)
    max_val = 1.0

    bars = ""
    for i, val in enumerate(values):
        x = 20 + i * (bar_width + gap)
        h = (val / max_val) * chart_height
        y = chart_height + 20 - h
        color = "#27ae60" if val > 0.6 else ("#f39c12" if val > 0.3 else "#e74c3c")
        bars += f'<rect x="{x}" y="{y:.1f}" width="{bar_width}" height="{h:.1f}" fill="{color}" rx="3"/>'
        bars += f'<text x="{x + bar_width/2}" y="{y - 5:.1f}" text-anchor="middle" font-size="10" fill="#2c3e50">{val:.2f}</text>'

    # Axis line
    axis = f'<line x1="10" y1="{chart_height + 20}" x2="{chart_width - 10}" y2="{chart_height + 20}" stroke="#bdc3c7" stroke-width="1"/>'

    # Grid lines
    grid = ""
    for pct in [0.25, 0.5, 0.75, 1.0]:
        y = chart_height + 20 - pct * chart_height
        grid += f'<line x1="10" y1="{y:.1f}" x2="{chart_width - 10}" y2="{y:.1f}" stroke="#ecf0f1" stroke-width="1" stroke-dasharray="3,3"/>'
        grid += f'<text x="5" y="{y + 3:.1f}" text-anchor="end" font-size="9" fill="#95a5a6">{pct:.0%}</text>'

    return f"""<svg viewBox="0 0 {chart_width} {chart_height + 40}" style="width:100%">
        {grid}
        {bars}
        {axis}
    </svg>"""


def _render_offtarget_svg(data: list[dict]) -> str:
    """Render an SVG grouped bar chart for off-target analysis."""
    if not data:
        return "<p>No off-target data available</p>"

    n_guides = len(data)
    bar_width = 20
    group_width = bar_width * 3 + 6
    chart_height = 180
    chart_width = max(300, n_guides * (group_width + 20) + 40)

    bars = ""
    for i, guide in enumerate(data):
        x_base = 30 + i * (group_width + 20)
        high = guide["high_risk"]
        mod = guide["moderate_risk"]
        low = guide["low_risk"]

        # High risk (red)
        h = min(high / 10, 1.0) * chart_height
        bars += f'<rect x="{x_base}" y="{chart_height + 20 - h:.1f}" width="{bar_width}" height="{h:.1f}" fill="#e74c3c" rx="2"/>'
        if high > 0:
            bars += f'<text x="{x_base + bar_width/2}" y="{chart_height + 20 - h - 4:.1f}" text-anchor="middle" font-size="9" fill="#e74c3c">{high}</text>'

        # Moderate risk (orange)
        h = min(mod / 10, 1.0) * chart_height
        bars += f'<rect x="{x_base + bar_width + 2}" y="{chart_height + 20 - h:.1f}" width="{bar_width}" height="{h:.1f}" fill="#f39c12" rx="2"/>'
        if mod > 0:
            bars += f'<text x="{x_base + bar_width + 2 + bar_width/2}" y="{chart_height + 20 - h - 4:.1f}" text-anchor="middle" font-size="9" fill="#f39c12">{mod}</text>'

        # Low risk (green)
        h = min(low / 10, 1.0) * chart_height
        bars += f'<rect x="{x_base + (bar_width + 2) * 2}" y="{chart_height + 20 - h:.1f}" width="{bar_width}" height="{h:.1f}" fill="#27ae60" rx="2"/>'
        if low > 0:
            bars += f'<text x="{x_base + (bar_width + 2) * 2 + bar_width/2}" y="{chart_height + 20 - h - 4:.1f}" text-anchor="middle" font-size="9" fill="#27ae60">{low}</text>'

        # Guide label
        label = f"#{i+1}"
        bars += f'<text x="{x_base + group_width/2}" y="{chart_height + 35}" text-anchor="middle" font-size="10" fill="#7f8c8d">{label}</text>'

    axis = f'<line x1="20" y1="{chart_height + 20}" x2="{chart_width - 10}" y2="{chart_height + 20}" stroke="#bdc3c7" stroke-width="1"/>'

    legend = """
        <rect x="20" y="0" width="12" height="12" fill="#e74c3c" rx="2"/>
        <text x="36" y="10" font-size="10" fill="#7f8c8d">High Risk</text>
        <rect x="100" y="0" width="12" height="12" fill="#f39c12" rx="2"/>
        <text x="116" y="10" font-size="10" fill="#7f8c8d">Moderate</text>
        <rect x="180" y="0" width="12" height="12" fill="#27ae60" rx="2"/>
        <text x="196" y="10" font-size="10" fill="#7f8c8d">Low Risk</text>
    """

    return f"""<svg viewBox="0 0 {chart_width} {chart_height + 50}" style="width:100%">
        {legend}
        {bars}
        {axis}
    </svg>"""


def generate_sgrna_diagram(
    sequence: str,
    candidates: list[dict],
    start_pos: int = 1,
) -> str:
    """Render an SVG diagram showing sgRNA positions on a target sequence."""
    seq_len = len(sequence)
    width = max(800, seq_len * 3)
    height = 120
    scale = width / seq_len

    # Draw sequence line
    line = f'<line x1="10" y1="60" x2="{width - 10}" y2="60" stroke="#34495e" stroke-width="2"/>'

    # Draw guide annotations
    annotations = ""
    for i, cand in enumerate(candidates[:20]):
        x = 10 + (cand["start"] - start_pos) * scale
        w = (cand["end"] - cand["start"] + 1) * scale
        y_offset = 35 if i % 2 == 0 else 70
        color = "#3498db" if cand["strand"] == "+" else "#e74c3c"
        annotations += f'<rect x="{x:.1f}" y="{y_offset}" width="{w:.1f}" height="8" fill="{color}" rx="2" opacity="0.8"/>'
        if w > 20:
            annotations += f'<text x="{x + w/2:.1f}" y="{y_offset - 3}" text-anchor="middle" font-size="8" fill="#7f8c8d">#{cand["rank"]}</text>'

    # Scale markers
    scale_marks = ""
    for pos in range(0, seq_len + 1, max(1, seq_len // 10)):
        x = 10 + pos * scale
        scale_marks += f'<line x1="{x:.1f}" y1="58" x2="{x:.1f}" y2="62" stroke="#bdc3c7" stroke-width="1"/>'
        scale_marks += f'<text x="{x:.1f}" y="80" text-anchor="middle" font-size="9" fill="#95a5a6">{start_pos + pos}</text>'

    return f"""<svg viewBox="0 0 {width} {height}" style="width:100%">
        {line}
        {annotations}
        {scale_marks}
    </svg>"""


def generate_offtarget_heatmap(
    guide_sequence: str,
    offtarget_sites: list[dict],
) -> str:
    """Render an SVG heatmap of off-target mismatches by position."""
    if not offtarget_sites:
        return "<p>No off-target sites found.</p>"

    guide_len = len(guide_sequence)
    cell_size = 24
    width = guide_len * cell_size + 80
    height = len(offtarget_sites) * cell_size + 40

    # Header (guide sequence)
    header = ""
    for i, base in enumerate(guide_sequence):
        x = 80 + i * cell_size
        header += f'<rect x="{x}" y="5" width="{cell_size}" height="{cell_size}" fill="#ecf0f1" stroke="#bdc3c7" stroke-width="0.5"/>'
        header += f'<text x="{x + cell_size/2}" y="{5 + cell_size/2 + 4}" text-anchor="middle" font-size="11" font-family="monospace" fill="#2c3e50">{base}</text>'

    # Mismatch cells
    rows = ""
    for row_idx, site in enumerate(offtarget_sites[:30]):
        y = 30 + row_idx * cell_size
        # Row label
        rows += f'<text x="75" y="{y + cell_size/2 + 4}" text-anchor="end" font-size="9" fill="#7f8c8d">{row_idx + 1}</text>'

        target = site.get("off_target_sequence", "")
        mismatches = set(site.get("mismatch_positions", []))

        for i in range(guide_len):
            x = 80 + i * cell_size
            if i in mismatches:
                # Mismatch - color by risk
                risk = site.get("risk_level", "low")
                color = {"high": "#e74c3c", "moderate": "#f39c12", "low": "#f1c40f"}.get(risk, "#f1c40f")
                rows += f'<rect x="{x}" y="{y}" width="{cell_size}" height="{cell_size}" fill="{color}" stroke="white" stroke-width="1"/>'
                base = target[i] if i < len(target) else "?"
                rows += f'<text x="{x + cell_size/2}" y="{y + cell_size/2 + 4}" text-anchor="middle" font-size="11" font-family="monospace" fill="white" font-weight="bold">{base}</text>'
            else:
                rows += f'<rect x="{x}" y="{y}" width="{cell_size}" height="{cell_size}" fill="#2ecc71" stroke="white" stroke-width="1"/>'
                base = target[i] if i < len(target) else "?"
                rows += f'<text x="{x + cell_size/2}" y="{y + cell_size/2 + 4}" text-anchor="middle" font-size="11" font-family="monospace" fill="white">{base}</text>'

    return f"""<svg viewBox="0 0 {width} {height}" style="width:100%;overflow:auto">
        {header}
        {rows}
    </svg>"""


def generate_efficiency_chart(values: list[float]) -> str:
    """Public alias for the efficiency bar chart renderer."""
    return _render_efficiency_svg(values)


def generate_evidence_radar(data: dict[str, float]) -> str:
    """Public alias for the evidence radar chart renderer."""
    return _render_radar_svg(data)
