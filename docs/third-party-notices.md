# Third-party notices

## crisprScore coefficient data

GeneImpact AI includes the CRISPRscan numerical coefficient data distributed by
the `crisprScore` project. The GeneImpact scoring implementation is independent
and is tested against the version-locked upstream result.

Copyright (c) 2022, Genentech, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Uribe-Salazar et al. zebrafish benchmark

The transformed benchmark in
`data/benchmarks/crisprscan-nhgri1-2022.json` is derived from Supplementary
Table 2 of:

Uribe-Salazar JM, Kaya G, Sekar A, Weyenberg K, Ingamells C, Dennis MY.
"Evaluation of CRISPR gene-editing tools in zebrafish." BMC Genomics 23, 12
(2022). https://doi.org/10.1186/s12864-021-08238-1

The article and supplementary material are distributed under the Creative
Commons Attribution 4.0 International License:
https://creativecommons.org/licenses/by/4.0/

Changes made by GeneImpact AI: percentages were divided by 100, guide
sequences were replaced with SHA-256 digests, selected source metadata was
normalized, and no raw sequencing data was copied.

The zero-overlap audit also reads canonical training oligos from Supplementary
Table 1 of Moreno-Mateos et al. (2015),
https://doi.org/10.1038/nmeth.3543. The training workbook is not redistributed.

## inDelphi external model and mouse validation evidence

GeneImpact AI does not redistribute the inDelphi source code or model files.
The external-result adapter records the official repository commit
`9ab67ca53ebb91e49aeb4530ec1e999ee9827ca1`. The upstream repository declares a
limited copyright license whose terms name eligible academic users and US
government research institutions for non-commercial research, and requires a
separate agreement for commercial or industrially sponsored use:
<https://github.com/maxwshen/inDelphi-model/blob/9ab67ca53ebb91e49aeb4530ec1e999ee9827ca1/LICENSE.txt>.

Aggregate mouse-embryo transfer evidence is attributed to:

Lkhagvadorj K, Okamura E, Taki T, et al. "Optimizing CRISPR precision in mouse
embryos via microhomology-mediated end joining-dominant targeting."
Communications Biology 9, 371 (2026).
<https://doi.org/10.1038/s42003-026-09771-z>

The article and supplementary workbook are marked CC BY-NC-ND 4.0. GeneImpact
AI does not copy or transform the workbook rows. It stores the source URL,
SHA-256 checksum, and aggregate factual metrics independently recomputed from
the unmodified source.

## Housden score and FlyRNAi service

GeneImpact AI does not redistribute the Housden coefficient table, CRISPOR
implementation, or FlyRNAi service code. The adapter accepts a researcher-
declared result from the official DRSC/FlyRNAi service and retains a checksum
of the original response.

The method is attributed to:

Housden BE, Valvezan AJ, Kelley C, et al. "Identification of potential drug
targets for Tuberous Sclerosis Complex by synthetic screens combining CRISPR-
based knockouts with RNAi." Science Signaling 8, rs9 (2015).
<https://doi.org/10.1126/scisignal.aab3729>

Official service: <https://www.flyrnai.org/evaluateCrispr/>.

The GeneImpact adapter contains independently written validation and
normalization logic. It does not include the upstream position coefficients.
