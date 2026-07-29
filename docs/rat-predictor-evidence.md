# Rat predictor evidence: public-data qualification

**Review date:** 2026-07-29

**Species:** *Rattus norvegicus*

**Intended platform use:** calibrated prediction or bounded benchmarking of CRISPR guide activity, off-target activity, and repair outcomes in rat research.

## Decision

No public rat dataset located in this review is sufficient to train or claim a broadly calibrated rat gene-editing predictor.

The strongest direct evidence is Anderson et al. (2018). Its public supplement supports a useful but very small, context-specific benchmark:

- 25 rat SpCas9 guide records across 12 named targets in the `rn5` assembly;
- quantitative mean on-target measurements for 16 of those rat guide records;
- candidate-site off-target testing, with exact sequences and per-animal allele counts for confirmed events.

This is enough to implement a **version-pinned rat embryo benchmark/importer**, provided that its limitations are visible in every result. It is not enough to expose an unqualified “rat prediction” capability. The exact in-vivo rat strain is not disclosed per project, the on-target sample is small, candidate off-target negatives were preselected by an algorithm, and no explicit dataset reuse licence was located.

The large Lewis-rat T-cell screen by Kendirli et al. is valuable biological evidence but is **not an editing-efficiency dataset**. Its labels are guide abundance and migration/fitness phenotypes, which are jointly affected by target-gene biology, cell state, viral representation, proliferation, and editing. It must not be relabelled as cleavage efficiency.

Recommended current readiness:

| Capability | Rat status |
|---|---|
| Context and assembly validation | Ready |
| General SpCas9 guide-activity predictor | Insufficient public data |
| Bounded in-vivo embryo SpCas9 benchmark | Candidate: Anderson 2018 |
| Genome-wide off-target predictor | Insufficient public data |
| Candidate-site off-target benchmark | Candidate with selection-bias warning: Anderson 2018 |
| Indel/repair-spectrum predictor | Insufficient public data |
| Large-rearrangement risk evidence | Hazard evidence only: CRISMERE |
| HDR/delivery-condition evidence | Method/transfer evidence only |
| Gene-level functional consequence evidence | Phenotype evidence only; not editing efficiency |

## Qualification criteria

A dataset qualifies for a calibrated predictor only when it provides:

1. exact species and, for in-vivo work, strain or stock;
2. reference assembly and target sequences or unambiguous coordinates;
3. nuclease, edit class, delivery method, tissue/cell type, and developmental stage;
4. public per-guide or per-outcome quantitative labels with denominators and replicates;
5. enough target diversity to separate sequence effects from locus, delivery, and gene-function effects;
6. raw or sufficiently granular processed measurements;
7. an explicit reuse licence compatible with the intended distribution;
8. an independent test partition with a documented audit against training-data overlap.

Method papers and hazard observations remain useful evidence, but they do not become prediction datasets merely because they report a percentage.

## Candidate comparison

| Candidate | Rat context | Public labels and scale | Licence/access | Qualified use | Decision |
|---|---|---|---|---|---|
| [Anderson et al. 2018](https://doi.org/10.1038/s41592-018-0011-5) | One-cell embryos; SpCas9 mRNA plus IVT sgRNA; `rn5`; project strain not disclosed | Rat subset: 25 guide records/12 targets; 16 guides with quantitative mean on-target values; candidate off-target loci and selected per-animal alleles | Supplements public; article is not NCBI open access; no explicit dataset licence located | Small external activity/off-target benchmark | **Best available direct candidate, not a calibrated predictor** |
| [Kendirli et al. 2023](https://doi.org/10.1038/s41593-023-01432-2) | Lewis-rat MBP-reactive CD4+ T cells; ex-vivo lentiviral Cas9/sgRNA; adoptive transfer | 87,690 guides; 21,410 genes, 396 miRNAs, 800 non-targeting controls; guide counts across tissues and replicates | Article CC BY; counts in GEO GSE232340 | Gene-function/migration phenotype screen | **Phenotype evidence only** |
| [Birling et al. 2017, CRISMERE](https://doi.org/10.1038/srep43331) | Sprague Dawley zygotes; Cas9 mRNA plus four guides per project | Five rat project configurations; founder deletion/inversion/duplication outcomes over 37 kb to 24.5 Mb intervals | CC BY; article and supplement public | Large structural-variant hazard/method evidence | **Not per-guide prediction data** |
| [Remy et al. 2017](https://doi.org/10.1038/s41598-017-16328-y) | SD/Crl zygotes; SpCas9 RNP electroporation; optional ssODN | Three loci (`Rosa26`, `Ephx2`, `FlnA`); condition-level embryo/pup NHEJ and knock-in counts | CC BY; article and supplement public | Delivery/HDR transfer evidence | **Too few loci; labels confounded by conditions** |
| [Challa et al. 2021](https://doi.org/10.1016/j.mex.2021.101419) | Sprague Dawley embryos; SpCas9 RNP electroporation and embryo culture | Two guides at one `Chrna7` locus; HMA/RFLP validation | CC BY; article public | Guide-validation workflow | **Method evidence only** |
| [Jin et al. 2023](https://doi.org/10.1016/j.stemcr.2022.11.012) | Primarily DA2B rat ESCs; RNP-assisted large-vector HDR; SD blastocysts for chimeras | Project-level targeted-colony rates; seven paired vector-only versus vector+CRISPR projects plus additional CRISPR-assisted projects | CC BY-NC-ND; supplement public; authors state no large dataset was generated | Large-HDR method evidence | **Too few, highly confounded project labels** |
| [Yoshimi et al. 2014](https://doi.org/10.1038/ncomms5240) | F344, DA, and Wistar rats; zygote microinjection; allele-specific CRISPR | Small number of guides at `Tyr`, `Asip`, and `Kit`; founder/embryo editing and correction outcomes | CC BY; article and supplement public | Allele-specific feasibility and locus evidence | **Insufficient target diversity** |
| [Savell et al. 2019](https://doi.org/10.1523/ENEURO.0495-18.2019) | Embryonic rat hippocampal neurons; dCas9-VPR activation, not nuclease editing | Target-gene expression and a 9-sample Bdnf RNA-seq series | CC BY; GEO GSE117961 and SRA SRP155892 | CRISPRa expression evidence | **Out of scope for cleavage or repair prediction** |

## Detailed evidence

### 1. Anderson et al. 2018: strongest bounded benchmark

The study reports 81 mouse and rat projects, 119 sgRNAs, 1,423 algorithm-predicted off-target sites, and 32 confirmed off-target sites. Twenty-one of the 119 guides showed at least one confirmed off-target event. The authors used IVT sgRNA and Cas9 mRNA in one-cell embryos and screened candidate sites selected by an MIT off-target algorithm. The primary article and author manuscript describe the design and limitations ([Nature Methods article](https://www.nature.com/articles/s41592-018-0011-5), [PMC author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6558654/)).

Exact public artifacts:

- [Supplementary Table 1: full guide and predicted off-target list](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-018-0011-5/MediaObjects/41592_2018_11_MOESM3_ESM.xlsx)
- [Supplementary Table 2: allele breakdown for G0 animals](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-018-0011-5/MediaObjects/41592_2018_11_MOESM4_ESM.xlsx)
- [Supplementary Table 5: raw mean on-target values](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-018-0011-5/MediaObjects/41592_2018_11_MOESM7_ESM.xlsx)
- [SRA study SRP124981](https://www.ncbi.nlm.nih.gov/sra/?term=SRP124981)

Inspection of the spreadsheets gives the following rat-specific scope:

- Supplementary Table 1 contains 25 `rn5` guide records across 12 named rat targets: `Il13`, `Map3k14`, `Trpa1`, `Usp30`, `Esr1`, `Il33`, `Jag1`, `Map4k1`, `Mapt`, `Ripk1`, `Rorc`, and `Tdo2`.
- Supplementary Table 5 contains quantitative mean on-target editing values for 16 rat guide records. The remaining rat entries do not form a uniform quantitative table.
- Supplementary Table 2 gives granular read/allele fractions for selected confirmed off-target-positive rat projects, rather than a complete per-animal record for all 25 rat guides.

Why it is useful:

- It is direct *R. norvegicus* embryo evidence rather than a human or mouse transfer assumption.
- It preserves guide and candidate-site sequences, genome build, molecular context, and several quantitative labels.
- It can test whether a generic ranking method transfers at all to rat embryos.

Why it is not a general predictor dataset:

- The paper does not identify the exact in-vivo rat strain for every project.
- Sixteen quantitative on-target guide records cannot support stable calibration across sequence motifs, chromatin states, genes, edit classes, or delivery conditions.
- Candidate off-target sites were selected by an existing algorithm; unobserved sites and algorithm-missed sites are not valid negatives. A model evaluated only on these sites would measure re-ranking of a selected candidate set, not genome-wide detection.
- Project objectives include different edit designs and founder-selection practices. The published mean is not a common, controlled per-cell cleavage assay.
- The SRA accession does not provide a complete, uniformly linked raw rat dataset for all 25 rat guide records.
- The article is not listed as open access by the NCBI OA service. The spreadsheets are downloadable, but no explicit dataset reuse licence was located. Redistribution inside a software release therefore needs rights review or author permission.

Required implementation safeguards:

- Pin the three source files by URL, retrieval date, byte count, and SHA-256.
- Import only rows explicitly marked `rn5`; never pool the mouse rows.
- Separate `on_target_activity` and `candidate_site_off_target` tasks.
- Retain project, guide, target, tested-site rank, sequencing method, and denominator.
- Never label untested genomic sites as negatives.
- Group all rows from the same guide and gene into one split.
- Run a training-overlap audit before calling the benchmark independent.
- Surface `strain_unknown`, `legacy_assembly_rn5`, `small_sample`, and `candidate_selection_bias` on every report.

**Qualification:** `transfer_evidence_only`. GeneImpact implements a
non-redistributing, checksum-pinned external-transfer evaluator; explicit
licence resolution is still required before publishing a transformed dataset.
This does not justify `usable_adapter` for general rat prediction.

### 2. Kendirli et al. 2023: large public screen, wrong label for activity

The study used Lewis-rat myelin-basic-protein-reactive CD4+ T cells, introduced Cas9 and a lentiviral genome-wide sgRNA library ex vivo, transferred the cells into naïve Lewis rats, and measured guide representation in blood, spleen, meninges, and CNS parenchyma. The genome-wide library contains 87,690 guides targeting 21,410 genes and 396 miRNAs plus 800 non-targeting controls, with approximately four guides per target. A 12,000-guide validation library was also used. The primary article is [open under CC BY](https://pmc.ncbi.nlm.nih.gov/articles/PMC10545543/).

Exact public artifacts:

- [GEO series GSE232340](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232340)
- [Genome-wide raw guide counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232340/suppl/GSE232340_GWscreen_RawCounts.txt.gz)
- [Genome-wide normalized guide counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232340/suppl/GSE232340_GWscreen_NormalizedCounts.txt.gz)
- [Validation-screen raw guide counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232340/suppl/GSE232340_ValidationScreen_RawCounts.txt.gz)
- [Validation-screen normalized guide counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232340/suppl/GSE232340_ValidationScreen_NormalizedCounts.txt.gz)

The abundance changes are not direct cleavage measurements. They combine editing probability with gene essentiality, activation, proliferation, survival, trafficking, bottlenecks, and sampling. Even plasmid-to-culture depletion cannot be interpreted as guide activity without independent edit-rate measurements for each guide.

**Qualification:** `phenotype_evidence_only`. Suitable for modelling T-cell migration biology or for testing gene-level functional priors; unsuitable as on-target efficiency, off-target, or repair-outcome training labels.

### 3. Birling et al. 2017: large-rearrangement hazard evidence

CRISMERE uses two guide pairs flanking a region to generate deletions, inversions, and duplications. In Sprague Dawley rat zygotes, the paper reports five project configurations spanning approximately 37.2 kb, 121.7 kb, 3.513 Mb, and 24.499 Mb regions. Founder-level structural outcomes are reported for `Cbs`, `Dyrk1a`, `Umodl1-Prmt2`, and `Lipi-Zfp295`. The study is [CC BY and publicly available](https://pmc.ncbi.nlm.nih.gov/articles/PMC5339700/); the [Europe PMC supplementary-file bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5339700/supplementaryFiles) contains the guide and junction material.

The outcome belongs to a four-guide/project configuration. It cannot be assigned to any one guide, and five configurations cannot identify sequence determinants separately from interval length, locus, embryo handling, and founder sampling.

**Qualification:** `hazard_evidence_only`. It supports an explicit warning that paired or multiplexed cuts can yield large deletions, inversions, duplications, and mosaic founders. It cannot calibrate their probability for a new design.

### 4. Remy et al. 2017: delivery and HDR transfer evidence

This study used Sprague-Dawley SD/Crl zygotes and SpCas9 RNP electroporation, with or without donor DNA. It reports condition-level embryo or pup outcomes for three loci (`Rosa26`, `Ephx2`, and `FlnA`) across voltage, pulse duration, RNP concentration, and donor configurations. The exact target sequences are disclosed in the paper. The [article is CC BY](https://pmc.ncbi.nlm.nih.gov/articles/PMC5707420/), and its [supplementary-file bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5707420/supplementaryFiles) is public.

The labels are valuable for estimating operational ranges under a named protocol, but three loci cannot separate guide-sequence effects from delivery and donor effects. The study does not expose a broad per-guide raw read-count dataset.

**Qualification:** `transfer_evidence_only` for SD/Crl zygote RNP electroporation; not a sequence predictor dataset.

### 5. Challa et al. 2021: guide-validation workflow

The MethodsX paper demonstrates an accessible embryo-culture assay in Sprague Dawley rats. It tests two guides at one `Chrna7` exon-11 locus and uses HMA and RFLP to detect non-homologous end joining and donor-mediated editing. NCBI identifies the article as [CC BY](https://pmc.ncbi.nlm.nih.gov/articles/PMC8374522/).

The study answers “can these particular guides cut in this embryo workflow?” It does not publish a diverse quantitative guide panel or a repair spectrum.

**Qualification:** `method_evidence_only`.

### 6. Jin et al. 2023: large HDR in rat ESCs

The study uses primarily DA2B rat embryonic stem cells, RNP-assisted homology-directed repair with large targeting vectors, drug selection, and SD blastocysts for chimera production. Project-level targeted-colony rates are reported, including seven vector-only versus vector+CRISPR comparisons and additional CRISPR-assisted targeting projects. The authors explicitly state that no large dataset was generated. The [article and associated data](https://pmc.ncbi.nlm.nih.gov/articles/PMC9860120/) are public, and the [supplementary-file bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9860120/supplementaryFiles) includes Tables S1-S7.

NCBI records the licence as **CC BY-NC-ND**, not CC BY. Commercial reuse, derivative redistribution, and inclusion of transformed supplementary material require legal review.

The targeting-rate label is a whole-workflow outcome affected by vector size, homology arms, modification class, locus, guide set, cell line, selection, and colony calling. The small number of projects cannot identify a guide-activity function.

**Qualification:** `method_evidence_only` for large-vector HDR in rat ESCs.

### 7. Yoshimi et al. 2014: allele-specific feasibility

The study demonstrates allele-specific editing and phenotype correction in F344, DA, and Wistar rats at a small number of loci, including `Tyr`, `Asip`, and `Kit`. It reports embryo/founder outcomes and exact allele-specific designs. The [article is CC BY](https://pmc.ncbi.nlm.nih.gov/articles/PMC4083438/), and the [supplementary-file bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4083438/supplementaryFiles) is public.

Its scientific value is high for allele-specific design, but the low locus count and heterogeneous edit objectives do not support a general guide or repair-outcome model.

**Qualification:** `transfer_evidence_only` for allele-specific feasibility; insufficient for calibration.

### 8. Savell et al. 2019: CRISPRa, not nuclease editing

This study develops dCas9-VPR transcriptional activation in embryonic rat hippocampal neurons. Public expression data include [GEO GSE117961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117961) and [SRA SRP155892](https://www.ncbi.nlm.nih.gov/sra/?term=SRP155892). The [article is CC BY](https://pmc.ncbi.nlm.nih.gov/articles/PMC6412672/).

Expression activation by a catalytically inactive Cas9 effector is a different task from nuclease cleavage, indel formation, HDR, or off-target editing. Its data should remain in a distinct CRISPRa capability domain.

**Qualification:** `out_of_scope` for the current predictor.

## Training-overlap and independence rules

“Public” is not equivalent to “independent.” Before any imported model or score is benchmarked:

1. inventory the model's named training datasets, supplementary files, and publication dates;
2. hash guide and target sequences after canonical orientation;
3. detect exact and near-duplicate protospacers, shared genes, and shared projects;
4. exclude every overlapping guide/project from the test set;
5. split remaining data by guide and gene, not by individual sequencing row;
6. report the final rat-only test count after exclusions;
7. refuse performance claims when the remaining sample is too small for a meaningful confidence interval.

Anderson et al. can be independent evidence only for a model demonstrably not trained or tuned on its projects. Kendirli et al. cannot become an editing-activity benchmark through a split strategy because its endpoint is a different biological label.

## Minimum acquisition target for a real rat predictor

The following is a proposed engineering acceptance target, not a universal biological law:

- at least 200-500 rat guides across at least 50 genes for an initial narrow activity model;
- one declared strain, assembly, nuclease, delivery method, tissue/cell context, and developmental stage per calibrated domain;
- amplicon-NGS edit fractions with raw edited/total read counts, biological replicates, and negative controls;
- exact protospacer/PAM and assay amplicon sequences;
- predefined quality filters and complete reporting of failed/low-coverage guides;
- a gene-grouped, laboratory-independent holdout;
- a separate, much larger multi-target dataset for repair spectra, retaining individual indel sequences and counts;
- explicit machine-readable reuse terms that permit the intended hosted and/or commercial use.

For off-target prediction, candidate sites must come from an assay that can discover sites independently of the scoring algorithm being evaluated, followed by targeted validation and documented detection limits. Algorithm-selected PCR panels alone cannot establish genome-wide recall.

## Blocking gaps and next actions

1. **Exact strain metadata:** request project-level rat strain and colony information for the Anderson dataset from the authors.
2. **Complete raw labels:** request a guide-to-sample manifest and raw rat amplicon data for all 25 rat guide records, including the nine without quantitative means in Supplementary Table 5.
3. **Reuse rights:** obtain explicit permission or a data licence for the Anderson spreadsheets before redistributing normalized derivatives.
4. **Scale:** seek a collaboration with rat transgenic cores to release a homogeneous, multi-gene embryo-editing dataset with failed guides retained.
5. **Repair outcomes:** acquire or generate rat-specific per-allele amplicon data across hundreds of targets; no public candidate located here approaches the scale needed for a calibrated repair-spectrum model.
6. **External validation:** reserve a laboratory- and gene-disjoint rat cohort before model fitting.

Until those gaps are closed, the scientifically defensible product claim is:

> Rat sequence, assembly, strain declarations, and external results can be validated and audited. Public rat embryo evidence can be benchmarked in a narrow SpCas9 context. General rat editing efficiency, off-target risk, repair outcome, and phenotypic safety are not yet calibrated predictions.
