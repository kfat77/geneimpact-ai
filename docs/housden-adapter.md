# Housden fruit-fly guide-ranking adapter

GeneImpact AI can normalize a Housden score obtained from the official
DRSC/FlyRNAi service. The adapter is deliberately limited to the evidence
domain used to establish the score: *Drosophila* S2R+ cell culture with
SpCas9 and a 20-nt guide.

It is not an embryo, germline, whole-animal, phenotype, welfare, or safety
predictor. An S2R+ result is rejected when a dossier declares a fruit-fly
embryo or in-vivo context.

## Why external import

The official service is live but does not expose a version. Available upstream
implementations and coefficient distributions have licensing constraints or
ambiguities that are not compatible with silently redistributing them under
this repository's MIT license. GeneImpact AI therefore does not include the
Housden coefficient table or upstream scoring code.

Researchers run the official service, retain its original response, calculate
the response file's SHA-256, and enter the reported score and checksum in the
input envelope. The normalized audit record keeps only the protospacer
SHA-256, not the sequence.

```bash
python -m geneimpact import-housden \
  --input examples/housden-fruit-fly-result.json \
  --output housden-audit.json

python -m geneimpact dossier \
  examples/dossier-fruit-fly-housden-request.json \
  --output fruit-fly-cell-dossier.json
```

The example checksum is a placeholder and must be replaced with the SHA-256 of
the retained official response.

## Interpretation

The Housden value is a ranking score, not an editing probability. The original
service help describes scores above 7.5 as high efficiency, while current
FlyRNAi result guidance recommends scores above 5.0. The audit record preserves
both statements instead of selecting one silently.

The score does not assess off-target effects. Pair it with a version-locked
reference search such as the GeneImpact CRISPRitz adapter, and retain empirical
validation for the exact delivery, cell line, locus, and laboratory.

Primary references:

- Housden et al. method: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4642709/>
- official scoring service: <https://www.flyrnai.org/evaluateCrispr/>
- official help: <https://www.flyrnai.org/evaluateCrispr/help.jsp>
