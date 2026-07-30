# GeneImpact AI — Researcher Guide

> **Research use only.** Predictions are computational and must be validated experimentally.
> This tool does not establish safety, authorize an edit, or replace ethics/biosafety/veterinary review.

## Quick Start

### Installation

```bash
git clone https://github.com/kfat77/geneimpact-ai.git
cd geneimpact-ai
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install openpyxl xlrd
pip install -e .
```

### 30-Second Pipeline

```bash
# 1. Prepare your target sequence in FASTA format
# 2. Run the full pipeline
geneimpact pipeline \
  --input target.fa \
  --chrom chr1 \
  --species mouse \
  --nuclease SpCas9 \
  --output report.json \
  --html report.html

# 3. Open report.html in your browser
```

## Core Commands

### 1. `pipeline` — Full End-to-End Analysis

Runs: sgRNA design → efficiency prediction → off-target search → evidence scoring → assessment → report.

```bash
geneimpact pipeline \
  --input target.fa \
  --chrom chr1 \
  --start 1000 \
  --end 5000 \
  --species mouse \
  --nuclease SpCas9 \
  --strain "C57BL/6J" \
  --genome-build GRCm39 \
  --edit-class knockout \
  --gene-essentiality 0.8 \
  --phenotype-severity 0.6 \
  --top-k 10 \
  --max-candidates 50 \
  --max-mismatches 4 \
  --min-efficiency 0.3 \
  --min-specificity 0.5 \
  --output report.json \
  --html report.html
```

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | FASTA file with target sequence |
| `--chrom` | first seq | Chromosome/sequence ID |
| `--start` | 1 | 1-based start position |
| `--end` | full length | 1-based end position |
| `--species` | mouse | Target species |
| `--nuclease` | SpCas9 | CRISPR nuclease |
| `--strain` | C57BL/6J | Animal strain |
| `--genome-build` | GRCm39 | Genome assembly |
| `--edit-class` | knockout | Edit type (knockout, knockin, base_editing) |
| `--gene-essentiality` | 0.0 | Gene essentiality score (0-1) |
| `--phenotype-severity` | 0.0 | Expected phenotype severity (0-1) |
| `--top-k` | 10 | Number of top guides to report |
| `--max-candidates` | 50 | Max sgRNA candidates to evaluate |
| `--max-mismatches` | 4 | Max mismatches for off-target search |
| `--output` | required | JSON report path |
| `--html` | optional | HTML visualization path |

### 2. `design-sgrna` — sgRNA Design Only

```bash
geneimpact design-sgrna \
  --input target.fa \
  --nuclease SpCas9 \
  --guide-length 20 \
  --max-candidates 50 \
  --output guides.json
```

### 3. `predict` — Efficiency Prediction for a Single Guide

```bash
geneimpact predict \
  --guide GAGTCTGCTGACAGAGCTCG \
  --species mouse \
  --context-35nt AAAAGAGAGTCTGCTGACAGAGCTCGAGGAAAAA \
  --output prediction.json
```

### 4. `offtarget` — Off-Target Search

```bash
geneimpact offtarget \
  --guide GAGTCTGCTGACAGAGCTCG \
  --reference genome.fa \
  --nuclease SpCas9 \
  --max-mismatches 4 \
  --output offtargets.json
```

## Supported Species

| Species | Key | Genome Build | Efficiency Model |
|---------|-----|--------------|-----------------|
| Mouse | `mouse` | GRCm39 | RuleSet2-Transfer |
| Rat | `rat` | GRCr8 | Rat-Transfer |
| Zebrafish | `zebrafish` | GRCz12tu | CRISPRscan (calibrated) |
| Fruit Fly | `fruit_fly` | BDGP6.54 | Drosophila-Heuristic |
| Rhesus Macaque | `rhesus_macaque` | T2T-MMU8v2.0 | Macaque-Transfer |
| Cynomolgus Macaque | `cynomolgus_macaque` | T2T-MFA8v1.1 | Macaque-Transfer |

## Supported Nucleases

| Nuclease | PAM | PAM Position |
|----------|-----|-------------|
| SpCas9 | NGG | 3' of guide |
| SaCas9 | NNGRRT | 3' of guide |
| Cas12a | TTTV | 5' of guide |

## Output Format

### JSON Report Structure

```json
{
  "pipeline_version": "1.0.0",
  "timestamp": "2026-07-30T...",
  "config": { ... },
  "study_context": { ... },
  "species_validation": { ... },
  "guides": [
    {
      "rank": 1,
      "guide_sequence": "GAGTCTGCTGACAGAGCTCG",
      "pam": "AGG",
      "strand": "+",
      "start": 100,
      "end": 119,
      "gc_content": 0.55,
      "efficiency": {
        "score": 0.85,
        "confidence": 0.35,
        "model": "RuleSet2-Transfer"
      },
      "indel_outcome": {
        "deletion_rate": 0.50,
        "insertion_rate": 0.11,
        "most_likely": "deletion",
        "predicted_size": -3
      },
      "offtarget": {
        "high_risk": 0,
        "specificity_score": 1.0,
        "top_sites": [...]
      },
      "evidence_scores": {
        "on_target_uncertainty": 0.55,
        "off_target_evidence": 0.0,
        "network_impact_evidence": 0.48,
        "welfare_relevance": 0.42
      },
      "assessment": {
        "concern_score": 0.55,
        "tier": "enhanced_review",
        "rationale": ["on-target uncertainty"]
      },
      "recommendation": "Recommended candidate; requires enhanced review."
    }
  ],
  "report_notice": "Research decision-support only..."
}
```

### Review Tiers

| Tier | Score | Meaning |
|------|-------|---------|
| `standard_review` | < 0.4 | Low concern signals |
| `enhanced_review` | 0.4 - 0.7 | Moderate concern, requires review |
| `high_concern_review` | ≥ 0.7 | High concern signals, careful review needed |

## Python API

```python
from geneimpact import (
    design_sgrnas, predict_efficiency, find_offtargets,
    run_pipeline, PipelineConfig, NucleaseType,
)
from geneimpact.provenance import StudyContext

# Design sgRNAs
result = design_sgrnas("ATCGATCGATCGATCGATCGGG", nuclease=NucleaseType.SPCAS9)
for candidate in result.candidates:
    print(f"{candidate.guide_id}: {candidate.guide_sequence} PAM={candidate.pam}")

# Predict efficiency
eff = predict_efficiency(result.candidates[0], species_key="mouse")
print(f"Efficiency: {eff.efficiency_score:.3f}")

# Search off-targets
report = find_offtargets(
    guide_sequence="GAGTCTGCTGACAGAGCTCG",
    reference_sequences={"chr1": "..."},
    nuclease=NucleaseType.SPCAS9,
    max_mismatches=3,
)
print(f"Specificity: {report.specificity_score:.3f}")

# Full pipeline
config = PipelineConfig(species="mouse", top_k=10)
study = StudyContext(
    species="mouse", strain_or_breed="C57BL/6J",
    genome_build="GRCm39", edit_class="knockout",
    evidence_snapshot="my_experiment_v1",
)
report = run_pipeline(
    sequence="ATCGATCG...",
    config=config,
    study_context=study,
    reference_sequences={"chr1": "..."},
)
report.to_json("output.json")
```

## Understanding the Results

### Efficiency Score (0-1)
Estimated probability of successful editing at the on-target site. Higher = better.

- **> 0.7**: High efficiency expected
- **0.3 - 0.7**: Moderate efficiency
- **< 0.3**: Low efficiency; consider alternative guides

### Specificity Score (0-1)
How unique the guide is in the reference genome. Higher = fewer off-targets.

- **> 0.8**: Highly specific, minimal off-target risk
- **0.5 - 0.8**: Moderate specificity
- **< 0.5**: Low specificity, significant off-target risk

### Evidence Scores (4 dimensions, 0-1 each)
- **on_target_uncertainty**: How uncertain we are about the on-target outcome
- **off_target_evidence**: Concern from off-target sites
- **network_impact_evidence**: Biological network disruption potential
- **welfare_relevance**: Animal welfare relevance

### Indel Outcome
Predicted distribution of insertion/deletion outcomes:
- `deletion_rate`: Fraction of edits expected to be deletions
- `insertion_rate`: Fraction expected to be insertions
- `no_edit_rate`: Fraction expected to have no edit
- `predicted_indel_size`: Most common indel size (negative = deletion)

## Limitations

1. **Efficiency models are heuristic** for non-zebrafish species. Always validate experimentally.
2. **Off-target search** uses exact matching with position-weighted penalties. For comprehensive analysis, use dedicated tools (CRISPRitz, Cas-OFFinder) and import results.
3. **Indel prediction** is a simplified model. For precise outcome prediction, use inDelphi or similar tools.
4. **No genome download**: The system processes user-provided FASTA files. Download reference genomes from Ensembl/NCBI separately.

## Integration with External Tools

GeneImpact AI can import results from external predictors:

```bash
# Import CRISPRscan scores (zebrafish)
geneimpact score-crisprscan --input request.json --output scores.json

# Import inDelphi results
geneimpact import-indelphi --input indelphi_result.json --output normalized.json

# Import CRISPRitz off-target results
geneimpact import-crispritz --metadata meta.json --targets targets.txt --output audit.json

# Import BE-Hive efficiency results
geneimpact import-behive-efficiency --input behive.json --output normalized.json
```

## FAQ

**Q: Can I use this for human gene editing?**
A: This tool is designed for animal genome editing research. Human use is not supported.

**Q: How accurate are the efficiency predictions?**
A: Zebrafish predictions use the published CRISPRscan model (validated). Other species use transfer heuristics with lower confidence (~0.35). Always validate experimentally.

**Q: Do I need to download a reference genome?**
A: Yes. Download FASTA files from Ensembl or NCBI. The tool processes user-provided sequences.

**Q: Can I run this offline?**
A: Yes. The prediction pipeline runs entirely offline. Only MGI/IMPC data fetching requires internet.

**Q: How do I cite this tool?**
A: See the project README for citation information.
