# Cnidarian HOX-DNA binding preferences reveal the shared logic of eumetazoan axial patterning

Carlos Sanchez-Higueras, Shaoxun Liu, Melanie Weilert, Julia Zeitlinger,
Harmen J. Bussemaker, James C-G Hombría, Richard S. Mann, and Matthew C. Gibson.
*Nature Communications*, in press.

This repository, which was developed by Shaoxun Liu under the supervision of Harmen Bussemaker at Columbia University, 
provides figure-generation code and numerical source data for all SELEX-related analyses in the paper. 

The official figure source is
[Figures.Rmd](Figures.Rmd). The [HTML report](Figures.html) is self-contained;
the [Markdown report](Figures.md) uses the accompanying panel images.

## Manuscript figures

| Figure | Computational panels |
| --- | --- |
| Figure 2 | a: HOX-PBX energy logos; b: position-7 triangle with T at the top, G at the left and A at the right; c-e: 12-mer, 8-mer and 4-mer trees |
| Figure 3 | a-c: pairwise relative-affinity scatter plots |
| Figure S2 | a-c: 12-mer, 8-mer and 4-mer human-HOX-plane projections |
| Figures S3-S5 | b: cross-species R-squared heatmaps |
| Figures S6, S8 and S10 | b: specificity tree; c: R-squared heatmap |
| Figure S7 | a-l: human and cnidarian affinity comparisons |
| Figure S9 | a-b: Drosophila and human energy logos |
| Figure S11 | b: monomer energy logos; c: monomer R-squared heatmap |
| Figure S12 | a-k: monomer and dimer affinity comparisons |
| Figure S13 | a-b: oral-aboral affinity comparisons |
| Figure S14 | Independent SELEX replicate energy logos |
| Figures S15-S19 | a: replicate specificity tree; b: replicate R-squared heatmap |

[The panel manifest](outputs/manuscript_source_data/panel_manifest.csv) maps
manuscript labels to all 60 PNG/vector-PDF exports, including the Figure 2
legend and the combined Figure S9 logo panel. Experimental images, embryo
schematics and assembled manuscript figures are provided separately.

## Software

The verified figure build uses R 4.3.1. Exact analysis-package versions are
recorded in [software_versions.csv](outputs/manuscript_source_data/software_versions.csv).
Rendering requires Pandoc, available with RStudio or as a standalone installation.
Validation and packaging use Python 3.9 or later and its standard library.

R dependencies are `FamilyCode`, `pheatmap`, `ggplot2`, `dplyr`, `ape`,
`phangorn`, `ggtree`, `Rtsne`, `ggrepel`, `EnvStats`, `stringr`, `digest`,
`rmarkdown` and `knitr`. The tested versions of `FamilyCode` and `ggtree`
are 1.0.1 and 3.8.2, respectively.

```r
install.packages(c("pheatmap", "ggplot2", "dplyr", "ape", "phangorn",
  "Rtsne", "ggrepel", "EnvStats", "stringr", "digest", "rmarkdown", "knitr",
  "BiocManager", "remotes"))
BiocManager::install("ggtree")
remotes::install_github("BussemakerLab/FCpackage")
```

## Reproduce the figure report

From the repository root:

```sh
Rscript --vanilla scripts/render_figures.R
python scripts/validate_publication.py
python scripts/package_publication.py
```

The supplied PSAMs are sufficient for this workflow; raw-read processing and
model fitting are separate upstream analyses. The renderer uses saved source
files, with no dependency on an RStudio project or saved `.RData` workspace.
If Pandoc is not on `PATH`, set `RSTUDIO_PANDOC` to the directory containing
its executable before rendering.

Outputs are:

- `Figures.html` and `Figures.md`: manuscript figure reports.
- `outputs/publication_panels/`: 60 PNG images and 60 vector PDFs.
- `outputs/manuscript_source_data/`: numerical source tables, Newick trees,
  model selections, input checksums, panel mappings and software versions.
- `dist/SELEX_manuscript_source_data.zip`: the figure workflow, required PSAMs,
  reports, panels and numerical source data, with SHA-256 checksums.
- `dist/SELEX_manuscript_panels.zip`: panel images labeled by manuscript figure.

Validation independently checks forward and reverse-complement window scores,
all scatter coordinates and fitted lines, the 36-coordinate projection-bar
calculation, source-input checksums, and consistency of the HTML/Markdown
reports with the exported images. Rendering also validates motif-position
selection, planar projections and triangle orientation.

## Calculations and source data

The 12-mer, 8-mer and 4-mer analyses use positions 1-12, 2-9 and 6-9,
respectively. The corresponding specificity signatures contain 36, 24 and
12 tetrahedron coordinates. Heatmap entries are squared Pearson correlations
between position-centered log-affinity coefficient vectors over the selected
positions. Functional trees use the tetrahedron signatures.

All 28 scatter panels use the same 2,000 unique 12-mer sequences generated
with random seed 42. Relative affinities sum contributions from sliding binding
windows in both the forward and reverse-complement orientations. Four N bases
are added at each end; N is assigned the least favorable coefficient at each
model position. Each window is scored relative to the optimal sequence for
the model, and the window affinities are summed without further rescaling.
Twelve-position dimer models have nine windows per orientation; eight-position
monomer models have thirteen. Source tables provide sequence identities,
separate strand scores, model coefficients, panel coordinates and all 196
ordinary least-squares fits with an intercept.

Projection bars use four complete pairs: AntHOX1-PBX, AntHOX1a-PBX (intein),
AntHOX6-PBX and AntHOX8-PBX. Each full 36-coordinate replicate difference is
projected onto the plane defined by the full 12-position HOXA1, HOXA5 and
HOXA9 signatures. The mean projected Euclidean length is the common bar
half-width; the full bar length is twice this mean. The same value is used
for the position-7 triangle and the 12-mer, 8-mer and 4-mer projections.
Replicate trees label the intein and GST AntHOX1a datasets separately.

## Upstream data and supporting analyses

Raw SELEX sequencing data are associated with
[BioProject PRJNA1439443](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1439443).
`CnDmelHox/`, `CnHoxMonomer/` and `JudithHumanHox/` contain the supplied PSAMs,
count tables and model-training notebooks. Upstream processing uses
[ProBoundTools](https://github.com/BussemakerLab/ProBoundTools) and
[PyProBound](https://github.com/BussemakerLab/PyProBound).

`Alignments/` and `HoxPhylogeneticTree.Rmd` support protein-sequence analyses.
`MotifAnalysis.Rmd` contains exploratory motif analyses. The official manuscript
SELEX panel definitions, model selections and scoring are in `Figures.Rmd` and
its `Rfunctions/` helpers.

## Protein quality control reports

[QC_Protein_Reports/](QC_Protein_Reports/) contains the identification reports
and quality certificates for recombinant proteins used in the SELEX experiments.
The [protein QC README](QC_Protein_Reports/README.pdf) describes the seven protein
folders and nine reports, including the LC-MS and gel-band identification reports,
sequence-coverage reports and certificates of analysis.
