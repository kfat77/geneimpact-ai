# Research aims and success criteria for animal genome-edit impact prediction

## Overall objective

Build a reproducible framework that ranks the plausible unintended consequences of a proposed animal genome edit, records the evidence behind each concern, and identifies edits that require deeper review before experimental progression.

This project does **not** promise a consequence-free edit. The expected result is a calibrated, species- and edit-class-specific prioritization of validation questions in a defined experimental context.

## Aim 1 — On-target consequence and uncertainty baseline

For each registered proposal, consolidate reference-genome version, target-locus annotation, expected functional change, and evidence gaps. Report the evidence provenance, uncertainty level, and applicability limits by species, strain/breed, and edit class.

**Success criterion:** the system's uncertainty is calibrated on held-out historical outcomes and does not label poorly evidenced edits as low concern.

## Aim 2 — Evidence-gated unintended-consequence prioritization

For a proposed edit, combine independently versioned evidence for candidate off-target effects, functional-network context, and welfare-relevant endpoints. The system must retain the source and strength of each signal rather than returning a black-box score.

**Success criterion:** high-priority concerns have better recall than the baseline in a held-out dataset and remain directionally consistent in an independent laboratory or dataset.

## Aim 3 — Generalizability, welfare, and failure analysis

Evaluate calibration and uncertainty in every adequately represented species, strain/breed, genome build, and edit class. Document missing evidence, selection bias, transfer failures, animal-welfare outcomes, and negative results.

**Success criterion:** subgroup results and limitations are released alongside aggregate results; unresolved material performance gaps or welfare risks block promotion beyond exploratory use.

## Decision rule

| Evidence level | Required evidence | Permitted wording |
| --- | --- | --- |
| Exploratory | Mechanistic or computational rationale | “Candidate concern” |
| Replicated | Held-out result plus independent dataset consistency | “Replicated concern” |
| Strong support | Replication plus corroborating outcome evidence | “Strongly supported concern” |

No output may use “safe”, “consequence-free”, or an authorization interpretation. High uncertainty is itself a reason for further review.
