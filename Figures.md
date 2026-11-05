---
title: "SELEX source data for manuscript figures"
output:
  html_document:
    df_print: paged
    toc: true
    toc_depth: 3
    self_contained: true
    keep_md: true
    mathjax: null
    highlight: null
---

This report presents the SELEX-derived motif logos, DNA-binding specificity
trees, projections and affinity comparisons for the manuscript figures.
Figures 2 and 3 contain the main SELEX panels; Figures S2-S19 contain the
corresponding supplementary panels. Numerical source data and software
versions are provided at the end.

<style>
body { color: #222; font-size: 15px; line-height: 1.6; }
.main-container { max-width: 1240px; padding: 24px 36px 64px; }
h1, h2, h3 { color: #111; line-height: 1.25; }
h1 { margin-top: 48px; font-size: 28px; }
h2 { margin-top: 32px; font-size: 22px; }
h3 { font-size: 18px; }
.title { margin-top: 10px; font-size: 32px; }
img { max-width: 100%; height: auto; }
table { font-size: 13px; margin-bottom: 24px; }
th, td { padding: 7px 10px !important; }
#TOC { background: #f7f8fa; padding: 20px 24px; margin: 24px 0 36px; }
#TOC > ul { columns: 2; column-gap: 48px; padding-left: 20px; }
#TOC li { break-inside: avoid; }
@media (max-width: 768px) {
  .main-container { padding: 16px; }
  #TOC > ul { columns: 1; }
}
@media print { .main-container { max-width: none; padding: 0; } }
</style>



## Binding models

Position-specific affinity matrices describe sequence-dependent binding
energies. For each sample, candidate models are evaluated against the HOX-PBX
seed, and the model, orientation and alignment with the lowest motif-matching
score define its specificity signature. The 12-mer spans positions 1-12;
the HOX-PBX 8-mer spans positions 2-9, and the HOX 4-mer spans positions 6-9.











# Figure 2

## Figure 2a: Nematostella HOX-PBX logos

Energy logos for the six Nematostella HOX-PBX complexes.

<img src="outputs/publication_panels/fig2a-dimer-logos-1.png" width="100%" />

## Figure 2b: position 7 T-A-G projection

Base preferences at position 7 projected onto the T-A-G plane. T is at the
top, G at the left and A at the right. Bars show the shared full-motif
replicate discrepancy projected onto the three-HOX reference plane, as
defined in [Source data](#replicate-discrepancy).

<img src="outputs/publication_panels/fig2b-position7-triangle-1.png" width="2160" />

## Figure 2c: 12-mer specificity tree

Unrooted hierarchical clustering of specificity signatures for positions 1-12.

<img src="outputs/publication_panels/fig2c-tree-12mer-1.png" width="2700" />

## Figure 2d: 8-mer specificity tree

Unrooted hierarchical clustering of specificity signatures for positions 2-9.

<img src="outputs/publication_panels/fig2d-tree-8mer-1.png" width="2700" />

## Figure 2e: 4-mer specificity tree

Unrooted hierarchical clustering of specificity signatures for positions 6-9.

<img src="outputs/publication_panels/fig2e-tree-4mer-1.png" width="2700" />

### Shared tree legend

<img src="outputs/publication_panels/fig2-tree-legend-1.png" width="1950" />

# Figure 3

Relative DNA-binding affinities are compared for 2,000 unique 12-mer sequences
generated with random seed 42 from seven 8-mer cores and random flanking bases.
For each sequence, affinities are summed over all binding windows in both the
forward and reverse-complement orientations. Scoring uses four N bases at
each end; an N contributes the least favorable coefficient at that position.
Each window is scored relative to the optimal sequence for its model, and
window affinities are summed without further rescaling. The same sequences
and scoring procedure are used in Figures S7, S12 and S13. Dimer models span
12 positions and monomer models span eight positions. Colors identify the
seven core sequences; lines show ordinary least-squares fits with an intercept.

## Figure 3a: AntHOX1a versus AntHOX6a

<img src="outputs/publication_panels/fig3a-anthox1a-vs-anthox6a-1.png" width="1950" />

## Figure 3b: AntHOX1a versus AntHOX8

<img src="outputs/publication_panels/fig3b-anthox1a-vs-anthox8-1.png" width="1950" />

## Figure 3c: AntHOX6a versus AntHOX8

<img src="outputs/publication_panels/fig3c-anthox6a-vs-anthox8-1.png" width="1950" />

# Supplementary Figure 2

DNA-binding specificity signatures projected onto the plane defined by human
HOXA1, HOXA5 and HOXA9. Panels a-c use positions 1-12, 2-9 and 6-9,
respectively. Each selected position contributes its X, Y and Z tetrahedron
coordinates. For the bars, each of four replicate differences retains all
36 coordinates and is projected onto the plane defined by the full
12-position HOXA1, HOXA5 and HOXA9 signatures. The mean of the four projected
Euclidean lengths is the common bar half-width in panels a-c; total bar
length is twice this mean. The same full-signature error calculation is
used for the 12-mer, 8-mer and 4-mer panels.

## Figure S2a: 12-mer projection

<img src="outputs/publication_panels/figS2a-projection-12mer-1.png" width="2160" />

## Figure S2b: 8-mer projection

<img src="outputs/publication_panels/figS2b-projection-8mer-1.png" width="2160" />

## Figure S2c: 4-mer projection

<img src="outputs/publication_panels/figS2c-projection-4mer-1.png" width="2160" />

# Supplementary Figure 3

Cross-species 12-mer R-squared heatmap.

**Figure S3b**

<img src="outputs/publication_panels/supp3-cross-species-12mer-r2-1.png" width="2700" />

# Supplementary Figure 4

Cross-species HOX-PBX 8-mer R-squared heatmap.

**Figure S4b**

<img src="outputs/publication_panels/supp4-cross-species-8mer-r2-1.png" width="2700" />

# Supplementary Figure 5

Cross-species HOX 4-mer R-squared heatmap.

**Figure S5b**

<img src="outputs/publication_panels/supp5-cross-species-4mer-r2-1.png" width="2700" />

# Supplementary Figure 6

Human and Nematostella HOX-PBX 12-mer clustering.

**Figure S6b**

<img src="outputs/publication_panels/supp6-human-cnidarian-tree-1.png" width="2550" />

**Figure S6c**

<img src="outputs/publication_panels/supp6-human-cnidarian-r2-1.png" width="2550" />

# Supplementary Figure 7

Pairwise synthetic-affinity comparisons among the three human anchors and the
three directive-axis Nematostella HOX-PBX motifs.

**Figure S7a-l**

<img src="outputs/publication_panels/supp7-human-cnidarian-scatters-1.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-2.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-3.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-4.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-5.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-6.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-7.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-8.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-9.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-10.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-11.png" width="32%" /><img src="outputs/publication_panels/supp7-human-cnidarian-scatters-12.png" width="32%" />

# Supplementary Figure 8

Representative Nematostella HOX-PBX 12-mer clustering.

**Figure S8b**

<img src="outputs/publication_panels/supp8-cnidarian-tree-1.png" width="2400" />

**Figure S8c**

<img src="outputs/publication_panels/supp8-cnidarian-r2-1.png" width="2100" />

# Supplementary Figure 9

Drosophila and human HOX-cofactor logos.

**Figure S9 (logo components)**

<img src="outputs/publication_panels/supp9-reference-logos-1.png" width="100%" />

# Supplementary Figure 10

Drosophila HOX-Exd 12-mer clustering.

**Figure S10b**

<img src="outputs/publication_panels/supp10-drosophila-tree-1.png" width="2400" />

**Figure S10c**

<img src="outputs/publication_panels/supp10-drosophila-r2-1.png" width="2100" />

# Supplementary Figure 11

Nematostella HOX monomer motifs and clustering.

**Figure S11b**

<img src="outputs/publication_panels/supp11-monomer-logos-1.png" width="100%" />

**Figure S11c**

<img src="outputs/publication_panels/supp11-monomer-r2-1.png" width="2100" />

# Supplementary Figure 12

Synthetic-affinity comparisons used to contrast monomer and HOX-PBX binding
relationships.

**Figure S12a-k**

<img src="outputs/publication_panels/supp12-monomer-dimer-scatters-1.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-2.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-3.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-4.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-5.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-6.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-7.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-8.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-9.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-10.png" width="32%" /><img src="outputs/publication_panels/supp12-monomer-dimer-scatters-11.png" width="32%" />

# Supplementary Figure 13

Oral-aboral class comparisons in Nematostella and Drosophila.

**Figure S13a-b**

<img src="outputs/publication_panels/supp13-oral-aboral-scatters-1.png" width="49%" /><img src="outputs/publication_panels/supp13-oral-aboral-scatters-2.png" width="49%" />

# Supplementary Figure 14

Replicate logos included in the manuscript reproducibility panel.

**Figure S14a-e**

<img src="outputs/publication_panels/supp14-replicate-logos-1.png" width="70%" />

# Supplementary Figure 15

Clustering of quality-controlled 12-mer replicate signatures with Drosophila
and human reference models.

**Figure S15a**

<img src="outputs/publication_panels/supp15-retained-12mer-tree-1.png" width="3000" />

**Figure S15b**

<img src="outputs/publication_panels/supp15-retained-12mer-r2-1.png" width="3300" />

# Supplementary Figure 16

Clustering of quality-controlled HOX-PBX 8-mer replicate signatures.

**Figure S16a**

<img src="outputs/publication_panels/supp16-retained-8mer-tree-1.png" width="3000" />

**Figure S16b**

<img src="outputs/publication_panels/supp16-retained-8mer-r2-1.png" width="3300" />

# Supplementary Figure 17

Eight-mer clustering of the intein AntHOX1a replicate pair, the AntHOX1,
AntHOX6 and AntHOX8 replicate pairs, and one dataset each for AntHOX6a and
AntHOX7, together with human and Drosophila references.

**Figure S17a**

<img src="outputs/publication_panels/supp17-replicate-subset-8mer-tree-1.png" width="3000" />

**Figure S17b**

<img src="outputs/publication_panels/supp17-replicate-subset-8mer-r2-1.png" width="3300" />

# Supplementary Figure 18

Clustering of quality-controlled HOX 4-mer replicate signatures.

**Figure S18a**

<img src="outputs/publication_panels/supp18-retained-4mer-tree-1.png" width="3000" />

**Figure S18b**

<img src="outputs/publication_panels/supp18-retained-4mer-r2-1.png" width="3300" />

# Supplementary Figure 19

Clustering of Nematostella replicate signatures across the complete 12-mer.

**Figure S19a**

<img src="outputs/publication_panels/supp19-retained-cnidarian-tree-1.png" width="2700" />

**Figure S19b**

<img src="outputs/publication_panels/supp19-retained-cnidarian-r2-1.png" width="2700" />

# Source data

## Replicate discrepancy

Each full 12-position motif has a 36-coordinate tetrahedron signature.
Replicate differences use the complete signatures of four pairs:
AntHOX1-PBX, AntHOX1a-PBX (intein), AntHOX6-PBX and AntHOX8-PBX.
The AntHOX1a GST pair is not part of this estimate.

The plane defined by the full 12-position human HOXA1, HOXA5 and HOXA9
signatures has orthonormal basis vectors u<sub>1</sub> and u<sub>2</sub>, each
with 36 coordinates. For each pair, the full difference vector Δ<sub>i</sub>
is projected onto this plane, giving the length
d<sub>i</sub> = √[(Δ<sub>i</sub> · u<sub>1</sub>)<sup>2</sup> +
(Δ<sub>i</sub> · u<sub>2</sub>)<sup>2</sup>].
All 12 positions are retained for this calculation in every panel.
The common bar half-width is the mean of d<sub>1</sub>, d<sub>2</sub>,
d<sub>3</sub> and d<sub>4</sub>, and the total horizontal or vertical bar
length is twice this mean.

The position-7 triangle and the 12-mer, 8-mer and 4-mer projections display
this shared three-HOX-plane replicate discrepancy for all Nematostella
points. Point coordinates use the positions and reference plane specified
for each panel. The source tables include the 36-coordinate error-plane
basis, full replicate differences and their projected lengths.




|Replicate 1                |Replicate 2                | Projected Euclidean distance|
|:--------------------------|:--------------------------|----------------------------:|
|AntHOX1-PBX-rep1           |AntHOX1-PBX-rep2           |                    0.0098459|
|AntHOX1a-PBX-rep1 (intein) |AntHOX1a-PBX-rep2 (intein) |                    0.0113708|
|AntHOX6-PBX-rep1           |AntHOX6-PBX-rep2           |                    0.0237321|
|AntHOX8-PBX-rep1           |AntHOX8-PBX-rep2           |                    0.1411488|


|Panel      | Bar half-width| Total bar length|
|:----------|--------------:|----------------:|
|Figure 2b  |      0.0465244|        0.0930488|
|Figure S2a |      0.0465244|        0.0930488|
|Figure S2b |      0.0465244|        0.0930488|
|Figure S2c |      0.0465244|        0.0930488|

## Projection coordinates

Coordinates are expressed in orthonormal bases of the planes used for the
corresponding panels. CSV tables contain the coordinates, selected features,
replicate identities and error-bar values. Newick trees and pairwise
R-squared matrices accompany the clustering panels. Tree-label tables link
Newick identifiers to the labels displayed in each figure.


### Figure 2b

|HOX model |        x|        y|
|:---------|--------:|--------:|
|AbdA      | 1.249391| 1.324834|
|AbdB      | 1.072931| 1.138388|
|AntHOX1   | 1.709283| 0.917260|
|AntHOX1a  | 1.695257| 1.467669|
|AntHOX6   | 0.801840| 0.603319|
|AntHOX6a  | 0.965217| 0.399023|
|AntHOX7   | 1.310314| 0.615581|
|AntHOX8   | 1.606131| 0.575602|
|Antp      | 1.227586| 1.170278|
|Dfd       | 1.906094| 0.267796|
|HoxA1     | 0.569909| 0.281616|
|HoxA5     | 1.535094| 0.498138|
|HoxA9     | 1.381922| 1.813929|
|Lab       | 0.511107| 0.262001|
|Pb        | 0.459544| 0.255626|
|Scr       | 1.642509| 0.639324|
|Ubx       | 1.457059| 1.734654|




### Figure S2a

|HOX model |         x|         y|
|:---------|---------:|---------:|
|AbdA      |  1.042315|  0.883264|
|AbdB      |  0.966708|  0.945368|
|AntHOX1   |  1.535876|  0.397532|
|AntHOX1a  |  1.476052|  0.769259|
|AntHOX6   |  0.435962|  0.162584|
|AntHOX6a  |  0.267176|  0.045920|
|AntHOX7   |  0.890526|  0.082284|
|AntHOX8   |  1.253561|  0.054220|
|Antp      |  0.987235|  0.789482|
|Dfd       |  1.556829| -0.169484|
|HoxA1     |  0.000000|  0.000000|
|HoxA5     |  1.158964|  0.000000|
|HoxA9     |  1.506722|  1.690896|
|Lab       |  0.077529|  0.112987|
|Pb        | -0.010120|  0.177303|
|Scr       |  1.139632| -0.040592|
|Ubx       |  1.411660|  1.172249|




### Figure S2b

|HOX model |         x|         y|
|:---------|---------:|---------:|
|AbdA      |  0.942049|  0.922914|
|AbdB      |  0.899512|  0.963678|
|AntHOX1   |  1.455467|  0.359097|
|AntHOX1a  |  1.406399|  0.733865|
|AntHOX6   |  0.425385|  0.190344|
|AntHOX6a  |  0.205797|  0.055126|
|AntHOX7   |  0.846965|  0.046324|
|AntHOX8   |  1.184614|  0.015518|
|Antp      |  0.923388|  0.875005|
|Dfd       |  1.504055| -0.279639|
|HoxA1     |  0.000000|  0.000000|
|HoxA5     |  1.112027|  0.000000|
|HoxA9     |  1.441434|  1.599613|
|Lab       |  0.056828|  0.139825|
|Pb        | -0.029062|  0.223606|
|Scr       |  1.070919| -0.029648|
|Ubx       |  1.348523|  1.223992|




### Figure S2c

|HOX model |         x|         y|
|:---------|---------:|---------:|
|AbdA      |  0.935432|  0.959039|
|AbdB      |  0.923627|  1.082350|
|AntHOX1   |  1.412784|  0.347826|
|AntHOX1a  |  1.412859|  0.807574|
|AntHOX6   |  0.382158|  0.136814|
|AntHOX6a  |  0.194611|  0.073704|
|AntHOX7   |  0.875178|  0.147342|
|AntHOX8   |  1.187124|  0.076011|
|Antp      |  0.923308|  0.929120|
|Dfd       |  1.525072| -0.219174|
|HoxA1     |  0.000000|  0.000000|
|HoxA5     |  1.101099|  0.000000|
|HoxA9     |  1.372767|  1.536264|
|Lab       |  0.032521|  0.083397|
|Pb        | -0.044442|  0.179914|
|Scr       |  1.142835|  0.167307|
|Ubx       |  1.350004|  1.337221|

## Affinity comparisons

CSV tables provide the seeded DNA sequences, separate forward and
reverse-complement affinities and their sums for every model used in the
scatter plots. Coordinates, sequence identities and fitted-line coefficients for all
28 scatter panels accompany the model-level affinity and coefficient tables.



## Software environment

The analysis uses R and the packages listed below.


|Software   |Version |
|:----------|:-------|
|R          |4.3.1   |
|FamilyCode |1.0.1   |
|pheatmap   |1.0.13  |
|ggplot2    |3.5.2   |
|dplyr      |1.1.4   |
|ape        |5.8.1   |
|phangorn   |2.12.1  |
|ggtree     |3.8.2   |
|Rtsne      |0.17    |
|ggrepel    |0.9.6   |
|EnvStats   |3.1.0   |
|stringr    |1.5.1   |
|digest     |0.6.33  |
|rmarkdown  |2.29    |
|knitr      |1.50    |
