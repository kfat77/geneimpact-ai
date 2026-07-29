# Mouse data contract

Mouse is the first registered benchmark species. The initial profile is intentionally narrow:

- species: `mouse` / `mus_musculus`
- taxonomy ID: `10090`
- reference assembly: `GRCm39`
- INSDC assembly accession: `GCA_000001635.9`
- reference strain boundary: `C57BL/6J` (`C57BL/6` accepted as an input alias)

## Authoritative sources

The platform uses [Ensembl REST species metadata](https://rest.ensembl.org/documentation/info/species) to detect assembly or accession drift. Functional and positional lookups will use versioned Ensembl endpoints.

Mouse phenotype and allele evidence will be sourced from [Mouse Genome Informatics data reports](https://www.informatics.jax.org/downloads/reports/index.html) and its documented programmatic access. MGI exposes phenotype, allele, genotype, strain, ontology, and orthology reports. Each imported snapshot must retain its source URL, retrieval timestamp, release identifier when available, file checksum, and licensing/provenance note.

IMPC results are planned as an independent high-throughput phenotype source. They are not treated as interchangeable with MGI annotations; dataset and procedure identifiers must remain attached to every observation.

## Applicability rule

An assessment is rejected if its genome build does not match the registered profile. A non-reference strain is allowed only with an explicit warning that strain-specific validation is still required. No coordinate or annotation is silently lifted between assemblies.

## Live metadata check

```bash
python -m geneimpact source-check --species mouse
```

The command verifies the locally registered taxonomy ID, assembly name, and assembly accession against the current Ensembl REST response.

## Reproducible MGI snapshot

```bash
python -m geneimpact snapshot-mgi \
  --report all-phenotypes \
  --output-dir data/mgi
```

The download is accompanied by a JSON manifest containing the exact source URL, UTC retrieval time, byte count, and SHA-256 checksum. The `data/` directory is ignored by Git because public source files may be large and restricted study data must never be committed.

## Normalize genome-editing evidence

```bash
python -m geneimpact normalize-mgi \
  --input data/mgi/MGI_PhenotypicAllele.rpt \
  --output data/mgi/endonuclease-alleles.jsonl
```

By default, normalization retains only MGI rows classified as `Endonuclease-mediated`. Each JSONL record preserves the allele, marker, Ensembl gene, allele attributes, publication reference, and high-level Mammalian Phenotype identifiers. A companion manifest records input/output SHA-256 values and record counts for audit.

These records are prior phenotype evidence, not direct measurements of a proposed edit and not causal labels. Missing phenotype annotations remain missing; they are not interpreted as “no phenotype.”

## Query IMPC knockout phenotypes

```bash
python -m geneimpact impc-gene \
  --gene Prkdc \
  --output data/impc/Prkdc-significant.json
```

The command uses the official [IMPC statistical-result Solr API](https://www.ebi.ac.uk/training/online/courses/impc-solr-api/introduction-to-the-solr-api-accessing-impc-data-programmatically/using-simple-solr-syntax-in-your-browser/) and requests only significant results for one marker. The result retains MP terms, effect sizes, P values, sex, zygosity, procedure, parameter, source query, and retrieval timestamp.

IMPC evidence is generated from standardized knockout mouse phenotyping. It can inform plausible phenotypic consequences, but it is not a prediction of every edit type, molecular context, strain, or off-target event.

## Build a bounded IMPC validation set

```bash
python -m geneimpact benchmark-impc \
  --gene Prkdc \
  --gene Kit \
  --output data/benchmarks/impc-validation.jsonl
```

This dataset retains both significant and non-significant tested procedure/parameter results. A non-significant result means that a specific IMPC comparison did not meet its statistical significance criteria; it does not mean the edited gene has no phenotype. Runs are limited to 50 genes to keep public API usage bounded.

## Build leakage-aware benchmark splits

```bash
python -m geneimpact benchmark-mgi \
  --input data/mgi/endonuclease-alleles.jsonl \
  --output-dir data/benchmarks/mgi-v1
```

All records for the same gene are assigned to one deterministic split, preventing the same gene from appearing in training and evaluation. The default builder also excludes MGI alleles marked as IMPC-derived, because using those records for training and IMPC results for validation would create source leakage.

The current benchmark contains positive high-level MP associations only. Missing annotations are not negative labels, so this dataset supports coverage and ranking evaluation but is not yet a valid binary safety classifier dataset.

## Establish the baseline hurdle

```bash
python -m geneimpact evaluate-baseline \
  --benchmark-dir data/benchmarks/mgi-v1 \
  --k 5
```

The baseline ranks MP terms only by their frequency in the training split, then reports macro Recall@K and gene hit rate without adapting to validation or test data. Its report is bound to the benchmark manifest checksum.

This is intentionally a weak but honest baseline. A future predictor must beat it on untouched genes. Calibration is marked not applicable because the current benchmark contains positive associations rather than positive and tested-negative labels.
