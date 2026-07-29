# Fruit-fly in-vivo predictor evidence

**Review date:** 2026-07-29

**Species:** *Drosophila melanogaster*

**Scope:** embryo, germline, somatic-tissue, and whole-animal CRISPR nuclease activity, repair outcomes, off-targets, multiplex interactions, and phenotype hazards.

## Decision

No public dataset located in this review supports a generally calibrated fruit-fly guide-activity, repair-outcome, off-target, or phenotypic-safety predictor.

One newly published resource does justify a narrow new capability. [Port et al. 2026](https://doi.org/10.1038/s41467-026-68434-z) provides exact multiplex Cas12a guide arrays, exact transgenic contexts, raw loss-of-heterozygosity (LOH) observations, quantitative wing and toxicity measurements, and a CC BY 4.0 licence. It supports a reproducible **array-level in-vivo Cas12a benchmark**. It does not support per-guide activity prediction because each tested reagent contains multiple guides, only one of 169 on-target arrays was inactive, and activity can arise from redundancy or synergy between guides.

For SpCas9, the best public in-vivo evidence is split across:

- a 66-single-guide, seven-gene transgenic study with mainly phenotype and germline-transmission labels;
- a 104-guide/mismatch-variant study whose quantitative sequence-efficiency panel covers only four visible-marker genes;
- a 1,738-line paired-guide library whose large screen labels viability rather than direct cutting.

These are useful bounded benchmarks and transfer evidence. None provides a large, homogeneous, per-guide amplicon-sequencing dataset across diverse genes.

Recommended readiness:

| Capability | Scientifically qualified status |
|---|---|
| Housden SpCas9 score in S2R+ cells | Existing narrow cell-culture ranking only |
| Single-guide SpCas9 activity in embryos, germline, or whole animals | `validation_candidate`; no calibrated predictor |
| Dual-guide SpCas9 whole-animal knockout | Array/resource audit plus phenotype evidence; no edit-probability calibration |
| Quadruple-guide Cas12a+ somatic activity | `usable_bounded_benchmark` at the array level |
| Per-guide Cas12a+ activity | Unsupported by the multiplex labels |
| Cas9 or Cas12a repair spectrum | Insufficient public multi-target data |
| Genome-wide off-target probability | Insufficient positive labels |
| Multiplex interaction or phenotype prediction | Phenotype/hazard evidence only |
| “Safe edit” prediction | Unsupported |

## Qualification standard

A predictor dataset must preserve exact stock/genotype, assembly, nuclease, guide-expression cassette, Cas-expression cassette, tissue and developmental context, guide sequence, quantitative edit labels and denominators, replicates, raw measurements, and reusable licence. It must include enough active and inactive reagents across enough genes to estimate calibration and must reserve a gene-, reagent-, and laboratory-disjoint test set.

Phenotype is not a direct editing label. Viability, tissue size, pigmentation, or migration may reflect guide activity, repair outcome, target-gene function, maternal RNA/protein perdurance, mosaicism, and selection. Likewise, multiplex-array labels cannot be assigned to their component guides.

## Candidate comparison

| Primary candidate | Exact context and scale | Public label | Licence/access | Qualified use | Decision |
|---|---|---|---|---|---|
| [Port et al. 2026](https://doi.org/10.1038/s41467-026-68434-z) | LbCas12a-D156R (`Cas12a+`); predominantly four-guide HD12aCFD arrays; transgenic larval tissues/whole flies; >800 genes | 525 arrays/2,197 guides screened on 2L and 490 arrays/1,957 guides on 2R; 168/169 in-window arrays reproducibly active; raw LOH, wing, apoptosis, and other source data | CC BY 4.0; source files public | Array-level activity, specificity-observation, toxicity and phenotype benchmark | **Strongest bounded in-vivo candidate; not per-guide calibration** |
| [Port et al. 2020](https://doi.org/10.7554/eLife.53865) | SpCas9; conditional paired-guide HD_CFD lines at `attP40`; BDGP6 design; somatic tissues and germline | 1,738 line/guide-pair mapping; whole-animal viability for 639 lines; direct ICE/Sanger editing for a small subset | CC BY; five XLSX supplements | Paired-guide resource and phenotype benchmark | **Not a direct activity predictor dataset** |
| [Port et al. 2015](https://doi.org/10.1534/g3.115.019083) | SpCas9; 66 transgenic single guides at `attP40`; act-Cas9/nos-Cas9; seven genes | Somatic phenotype for 66 guides; quantitative germline transmission mainly for 5 `yellow` and 18 `ebony` guides | CC BY 4.0; supplement public | Narrow single-guide in-vivo benchmark | **Useful but only seven genes; phenotype confounding** |
| [Ren et al. 2014](https://doi.org/10.1016/j.celrep.2014.09.044) | SpCas9 transgenic embryos; injected U6b guide plasmids; visible-marker loci | 104 tested guide/mismatch configurations; quantitative on-target panel includes 27 `white`, 4 `vermilion`, 4 `ebony`, and 4 `yellow` guides; germline mutation rates | Public article/supplements; NCBI does not list it as open access; no explicit dataset licence located | Narrow embryo/germline sequence benchmark and mismatch study | **Best early Cas9 design evidence; target diversity inadequate** |
| [Housden et al. 2015](https://doi.org/10.1126/scisignal.aab3729) | SpCas9 in Drosophila S2R+ cell culture; luciferase reporters | 75 reporter-targeting guides, three biological replicates; relative reporter mutation activity | Public author manuscript and official FlyRNAi service; no clear reusable coefficient/data licence located | S2R+ ranking only | **Explicitly out of scope for in-vivo prediction** |
| [Champer et al. 2017](https://doi.org/10.1371/journal.pgen.1006796) | SpCas9 homing drives at `yellow`; nanos- and vasa-Cas9; `w1118`, Canton-S, and five Global Diversity Lines | Individual-cross progeny counts, drive conversion, germline/embryonic resistance formation, selected resistance-allele sequences | CC BY; XLSX and PDF supplements | Repair-timing, maternal carryover, stock-dependence and gene-drive hazard evidence | **Strong repair/hazard evidence at one locus; not general repair prediction** |
| [Bassett et al. 2013](https://doi.org/10.1016/j.celrep.2013.06.020) | Cas9 mRNA plus in-vitro-transcribed sgRNA injected into syncytial embryos; `yellow` and `white` | Founder mutation and germline-transmission rates; targeted sequence/off-target checks | CC BY-NC-ND; article/supplement public | Historical embryo-delivery transfer evidence | **Two genes; method evidence only** |
| [Yu et al. 2013](https://doi.org/10.1534/genetics.113.153825) | Embryo Cas9/gRNA targeting seven loci across euchromatin and heterochromatin | Locus-level founder/germline outcomes, up to 100% at selected loci | Public article/supplements; NCBI does not list it as open access | Historical germline feasibility | **Seven loci; insufficient for calibration** |
| [Yang et al. 2024](https://doi.org/10.1016/j.heliyon.2024.e29061) | nos-Cas9 transgenic background; two `white` guides; pooled >50 edited and >50 comparator males for WGS | Pooled WGS variants; SRA PRJNA1083262 | CC BY-NC-ND; raw WGS public | Hypothesis-generating genome-instability observation | **Cannot attribute variants to off-target cutting; hazard evidence only** |

## 1. Port et al. 2026: valid array-level Cas12a benchmark

The study uses a D156R variant of *Lachnospiraceae bacterium* Cas12a and HD12aCFD transgenes designed to encode four guides against independent sites in a common exon. The published library contains 838 rows with four nonblank spacers and seven rows with three. The guide array is inserted at a defined landing site and Cas12a is supplied from controlled transgenes. Exact genotypes are published in Supplementary Data 3; the library table gives line ID, up to four spacer sequences, FlyBase gene ID, and gene symbol.

Exact artifacts:

- [Supplementary Data 1: 845 HD12aCFD array records and up to four component spacers](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-68434-z/MediaObjects/41467_2026_68434_MOESM3_ESM.csv)
- [Supplementary Data 2: plasmid and primer sequences](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-68434-z/MediaObjects/41467_2026_68434_MOESM4_ESM.xlsx)
- [Supplementary Data 3: exact Drosophila genotypes](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-68434-z/MediaObjects/41467_2026_68434_MOESM5_ESM.xlsx)
- [Source Data workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-68434-z/MediaObjects/41467_2026_68434_MOESM9_ESM.xlsx)
- [Europe PMC complete supplementary bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12827956/supplementaryFiles)
- [VDRC stock collection](https://stockcenter.vdrc.at/control/main)

The source workbook is unusually valuable. It contains:

- 8,197 raw disc observations for the 2L LOH screen (`Fig5i`);
- line-level LOH category counts for the 2R screen (`Fig5j`);
- 3,699 and 2,725 wing observations in the main Cas12a/Cas9 comparisons (`Fig6g`, `Fig6h`);
- per-disc apoptosis, LOH gain/loss, proliferation, protein-loss, vein, and wing measurements.

The LOH assay screened 525 lines encoding 2,197 guides on 2L and 490 lines encoding 1,957 guides on 2R. For the on-target intervals between reporter and centromere, 168 of 169 arrays caused reproducible LOH. Lines targeting outside the interval serve as a search for unintended cutting within the monitored interval; 16 showed sporadic observations, but none reproduced.

What this supports:

- reproducible import and audit of exact three- or four-guide arrays;
- continuous or ordinal **array-level LOH observation** benchmarking in the reported larval reporter context, with per-array interval relationships unresolved until an interval manifest is obtained;
- comparison of Cas12a quadruple arrays with the named Cas9 single/dual-guide lines;
- empirical hazard modules for apoptosis after simultaneous breaks, nearby-element effects, LOH, and multiplex-associated deletions.

What it does not support:

- a label for any individual component guide;
- a calibrated probability of activity, because only 1/169 in-window arrays was inactive;
- genome-wide off-target recall, because only 33% of the genome was monitored and there were no reproducible positive off-target arrays;
- generalization to SpCas9, embryo injection, germline editing, another Cas12a variant, another expression cassette, or a wild stock;
- repair-spectrum prediction: sequence-resolved repair was limited to selected reporter/endogenous loci rather than hundreds of independent targets.

**Qualification:** implement as `usable_bounded_benchmark` with task name such as `dmel_cas12a_array_loh_evidence_v1`. Every result must say “multiplex-array LOH observation,” not “guide efficiency,” and a zero score must not be called inactive while the target-to-interval relationship is unresolved.

## 2. Port et al. 2020: large paired-guide resource, phenotype labels

The HD_CFD collection uses pCFD6 to express two SpCas9 guides conditionally. The guide library was designed against BDGP6, integrated at `attP40`, and combined with Gal4/Cas9 lines for tissue-specific, germline, or ubiquitous mutagenesis.

Exact artifacts are contained in the [Europe PMC supplementary bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7062466/supplementaryFiles):

- `elife-53865-supp1.xlsx`: 1,738 line records with both guide sequences and target FlyBase IDs;
- `elife-53865-supp2.xlsx`: viability results for 639 tested lines, including curated known-lethal and known-viable subsets;
- `elife-53865-supp3.xlsx`: uncharacterized essential-gene results;
- `elife-53865-supp4.xlsx`: oligonucleotides;
- `elife-53865-supp5.xlsx`: exact strains.

In the ubiquitous screen, 290/639 guide-pair lines produced no viable offspring, 53 were semi-lethal, and 269 produced viable adults. Among 210 lines targeting genes curated as essential, 23 gave viable offspring; the paper shows that some apparent false negatives still edited efficiently and attributes part of the discrepancy to maternal transcript or protein perdurance. Among 54 lines targeting genes curated as viable, five were lethal and one semi-lethal.

These labels measure the product of editing, repair, gene essentiality, maternal contribution, mosaicism, and developmental timing. Direct ICE/Sanger editing measurements cover only a small selected subset and are not supplied as a uniform 1,738-line read-count matrix.

**Qualification:** public reagent/resource adapter and whole-animal phenotype benchmark, not an edit-efficiency predictor.

## 3. Port et al. 2015: bounded single-guide in-vivo evidence

This CC BY study generated 66 transgenic U6:3 single-guide lines targeting seven genes, with every guide inserted at `attP40` and crossed to a strong act-Cas9 source. Sixty-five produced a somatic phenotype; the inactive line carried a target-site polymorphism. Quantitative germline transmission was reported mainly for five `yellow` and eighteen `ebony` guides. Six of the eighteen `ebony` guides and all five `yellow` guides transmitted nonfunctional alleles to more than 93% of offspring; other `ebony` guides ranged from 13% to 90%.

The exact guide/phenotype material is in the [Europe PMC supplementary bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4502383/supplementaryFiles).

This is valuable because it uses individual guides in whole animals and controls transgene position. It remains too narrow for a general sequence model: numerical germline labels are concentrated in two pigmentation genes, while most of the 66-guide evidence is categorical phenotype. In-frame indels can leave gene function intact, so “no phenotype” is not equivalent to “no cutting.”

**Qualification:** a narrow single-guide external benchmark after digitizing Table S2; not a calibrated probability model.

## 4. Ren et al. 2014: early in-vivo Cas9 sequence evidence

Ren et al. systematically altered guide length, mismatch positions, concentration, and target sequence in Cas9-transgenic embryos. The paper reports 104 tested guide/mismatch configurations. Its main quantitative efficiency analysis comprises 39 on-target guides across only four visible-marker genes: 27 `white`, four `vermilion`, four `ebony`, and four `yellow`. Heritable mutation rate is mutant F1 divided by all observed F1. The study also tested four-guide multiplex injection and reported concentration/fertility trade-offs.

This evidence is direct and useful, but guide concentration and mismatch variants are mixed with distinct on-target sequences. The four-gene concentration creates severe locus/phenotype dependence. Candidate off-target statements are based on selected mismatch constructs or predicted similar sites, not an unbiased genome-wide assay.

The NCBI OA service does not classify the article as open access. Public access to the paper and supplements is not an explicit machine-readable data licence; redistribution requires rights review.

**Qualification:** `validation_candidate` for an embryo/germline SpCas9 ranking benchmark; not enough genes for general calibration.

## 5. Housden 2015 must remain cell-line-only

The Housden activity matrix was derived from 75 guides targeting luciferase reporters in S2R+ cells, with three biological replicates and a relative reporter activity endpoint. It was compared with guide efficiencies from earlier Drosophila publications and is exposed through the [official FlyRNAi scoring service](https://www.flyrnai.org/evaluateCrispr/).

This is not embryo, germline, somatic-tissue, or whole-animal evidence. Reporter plasmid accessibility, cell-line karyotype, transfection, repair, and expression differ from transgenic or injected flies. The score is a ranking, not an edit probability.

Training-overlap risk is material: the Housden paper used prior Drosophila studies for external comparison. Those same guides must not be presented as an untouched independent test without guide-sequence overlap analysis.

**Qualification:** retain the existing S2R+ adapter boundary and reject it for declared in-vivo dossiers.

## 6. Champer et al. 2017: repair timing and stock dependence

The study publishes individual-cross progeny phenotypes and calculated drive parameters for nanos- and vasa-Cas9 homing drives targeting `yellow`. It compares `w1118`, Canton-S, and five Global Diversity Lines. Germline conversion, functional/nonfunctional resistance, embryo mosaicism, and selected indel sequences are reported.

Exact artifacts:

- [S1 Dataset and appendices bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5518997/supplementaryFiles)
- within the bundle, `pgen.1006796.s004.xlsx` contains individual-cross phenotypes and calculations;
- `pgen.1006796.s003.pdf` contains guide/primer and resistance-allele sequences.

The results show that maternal Cas9 carryover can create post-fertilization resistance alleles and that embryo resistance formation varies strongly across genetic backgrounds. These are important warnings against treating an editing percentage as a fixed property of a guide.

All observations use one `yellow` locus and a gene-drive construct. Selection, homing, end joining, maternal deposition, promoter, and locus biology are inseparable.

**Qualification:** `hazard_evidence_only` for repair timing, mosaicism, stock dependence, and gene-drive containment.

## 7. Early germline studies: feasibility, not prediction

[Bassett et al. 2013](https://doi.org/10.1016/j.celrep.2013.06.020) injected Cas9 mRNA and in-vitro-transcribed guides into syncytial embryos and reported founder/germline outcomes at `yellow` and `white`. The article is CC BY-NC-ND. [Yu et al. 2013](https://doi.org/10.1534/genetics.113.153825) targeted seven loci and reported germline efficiencies up to 100%.

Both establish embryo and germline feasibility and preserve useful protocol denominators. Their locus counts are too small, conditions are heterogeneous, and labels are not large uniform amplicon-sequencing matrices.

**Qualification:** `method_evidence_only`.

## 8. Yang et al. 2024: do not treat pooled WGS differences as off-target labels

The study used two `white` guides in a nos-Cas9 transgenic genotype and pooled more than 50 edited and more than 50 comparator males for WGS. Raw data are at [NCBI BioProject PRJNA1083262](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1083262). The article is CC BY-NC-ND.

The reported between-pool variants cannot distinguish Cas9-induced mutations from inherited stock variation, segregating variants, drift, recombination, bottleneck effects, or pooling. There are only two guides at one locus and no per-founder matched parental design. Absence of overlap with predicted Cas-OFFinder sites also does not establish genome-wide specificity.

**Qualification:** hypothesis-generating genome-instability observation only; never a positive or negative off-target training set.

## Training-overlap and split rules

1. Canonicalize and hash every spacer, target, and array before model fitting.
2. Keep all guides from one array, gene, transgenic line, and publication in one split.
3. Do not split the four Cas12a guides in one HD12aCFD array across train and test.
4. Treat Port 2026 comparisons with 154 BDSC and 103 VDRC Cas9 lines as overlapping with their source libraries.
5. Audit Housden against Ren, Port, Bassett, and Yu sequences because earlier fly datasets were used in score evaluation.
6. Remove target-site polymorphism cases from sequence-only evaluation or represent the actual stock haplotype.
7. Report performance separately for nuclease, promoter, delivery, tissue, stage, and stock.
8. Never derive “no cutting” from viability, normal morphology, or the absence of a visible phenotype.

## Implementable recommendation

Implement two narrowly named evidence imports:

1. **`dmel_cas12a_array_loh_evidence_v1`**
   - ingest Supplementary Data 1, Supplementary Data 3, and the Source Data workbook;
   - retain all three or four component guides as one indivisible reagent;
   - expose raw and mean LOH score, reporter, disc count, genotype, source row,
     and `interval_relationship=unresolved` until a qualified interval manifest
     can supply the chromosome relationship;
   - classify 168/169 in-window arrays as observed reproducible activity;
   - use outside-interval arrays only for monitored-region specificity observations, not genome-wide negatives;
   - report discrimination metrics as not applicable unless future inactive arrays supply a meaningful negative class.

2. **`dmel_spcas9_in_vivo_transfer_benchmarks_v1`**
   - digitize Port 2015 single-guide germline outcomes and Ren 2014 heritable mutation rates;
   - keep the studies separate because delivery and Cas-expression contexts differ;
   - restrict claims to ranking/transfer evidence;
   - do not merge phenotype-only HD_CFD viability labels into edit-efficiency training.

Continue to reject Housden results for embryo, germline, or whole-animal contexts.

## Acquisition blockers

- No public homogeneous SpCas9 dataset with hundreds of individual guides, many genes, exact stock haplotypes, and per-guide amplicon read counts was located.
- The Cas12a resource has excellent scale but multiplex labels and almost no inactive arrays, preventing per-guide attribution and probability calibration.
- No fruit-fly repair dataset approaches the target diversity required for sequence-resolved repair-spectrum prediction.
- Available off-target studies either have no reproducible positives, monitor only part of the genome, test predicted candidates, or lack matched-founder controls.
- Public phenotype screens do not provide counterfactual safety labels and cannot predict animal welfare or absence of unintended effects.
- Several historically important supplements lack explicit redistribution licences compatible with packaging normalized data in a software distribution.

Until these gaps are closed, the defensible product statement is:

> GeneImpact AI can audit fruit-fly reference and experimental context, import the Housden S2R+ ranking within its cell-line domain, and benchmark a declared multiplex Cas12a array against public in-vivo LOH evidence. General in-vivo guide efficiency, repair outcome, off-target probability, multiplex phenotype, and safety remain uncalibrated research outputs.
