# Fruit-fly Cas12a array evidence audit

GeneImpact AI v0.15.0 can verify and summarize the Port et al. 2026
*Drosophila melanogaster* Cas12a source files. This is an in-vivo,
array-level evidence workflow. It is not a guide-design model or a safety
prediction.

## Qualified biological scope

- organism: *Drosophila melanogaster*;
- reagent: HD12aCFD arrays containing three or four 23-nt spacers;
- nuclease: LbCas12a-D156R;
- evidence: wing-disc loss-of-heterozygosity observations in the published
  reporter contexts;
- source: Port et al. 2026, DOI
  `10.1038/s41467-026-68434-z`.

The verified library contains 845 arrays and 3,373 component spacers. The
source workbook contains 8,197 `Fig5i` rows and 600 `Fig5j` rows, of which 490
have scored category counts covering 8,478 discs.

## Obtain the source files

Download these files directly from the publisher:

1. Supplementary Data 1, the HD12aCFD library CSV;
2. Supplementary Data 3, the genotype workbook;
3. the Source Data workbook.

The command rejects any file whose bytes do not match the pinned publisher
checksum. The repository does not redistribute these files.

## Run the audit

```bash
python -m geneimpact audit-fruit-fly-cas12a-evidence \
  --library 41467_2026_68434_MOESM3_ESM.csv \
  --genotypes 41467_2026_68434_MOESM5_ESM.xlsx \
  --source-data 41467_2026_68434_MOESM9_ESM.xlsx \
  --output fruit-fly-cas12a-audit.json
```

Add `--line-id HD12aCFD0001` to include one array lookup. The lookup reports
the FlyBase target, number of component spacers, SHA-256 hashes of those
spacers, screen membership, and array-level LOH summaries. It does not emit
raw spacer sequences.

## Interpretation boundary

The paper reports 168 active arrays among 169 arrays whose targets lie within
the monitored interval. That extreme imbalance does not support useful
discrimination or probability calibration. Furthermore, the experimental
unit is a multiplex array: redundancy or cooperation between its component
guides cannot be disentangled.

Consequently:

- never label one component guide from an array-level observation;
- never interpret LOH as a repair-spectrum, phenotype, welfare, or safety
  probability;
- never treat outside-interval observations as genome-wide off-target recall;
- never transfer the result to another nuclease, stock, tissue, stage, or
  delivery system without new evidence;
- retain the publication's 525/490 screen-line and 168/169 on-target
  aggregates as publisher claims unless an interval manifest allows exact
  reconstruction.

The machine-readable status is `usable_bounded_benchmark`, while
`predictive_adapter_available` remains false for this evidence record.
