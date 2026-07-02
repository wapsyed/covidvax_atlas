# COVID-19 Vaccination Atlas

[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](LICENSE)
[![Preprint](https://img.shields.io/badge/medRxiv-2024.05.22.24307755v3-blue)](https://www.medrxiv.org/content/10.1101/2024.05.22.24307755v3)

A systems vaccinology analysis of the transcriptomic immune response to COVID-19 vaccination and infection.

---

## About

The COVID-19 Vaccination Atlas provides an integrative characterization of the immune system's transcriptional response across five vaccine platforms and natural SARS-CoV-2 infection. The study comprises **681 bulk RNA-seq samples from 343 participants**, covering healthy vaccinated individuals, infected patients with and without prior vaccination, and recipients of vaccines targeting other pathogens (influenza, smallpox, yellow fever).

**Five COVID-19 vaccines analyzed:**

| Vaccine | Brand Name | Platform |
|---------|-----------|----------|
| BBIBP-CorV | Covilo® | Inactivated virus |
| ZF2001 | Zifivax® | Protein subunit |
| ChAdOx1 | Vaxzebria® / Covishield® | Viral vector |
| mRNA-1273 | Spikevax® | mRNA |
| BNT162b2 | Comirnaty® | mRNA |

**Key markers** distinguishing vaccine types and infection: `LHFPL2`, `FYN`, `TYMP`, `JUP`, `SERINC5`, `KIT`, `GRAMD1C`.

---

## Interactive Tools

| Tool | Description | Link |
|------|-------------|------|
| **COVID-19 Gene Explorer** | Explore gene expression across COVID-19 vaccination and infection conditions | [Launch App](https://wapsyed.shinyapps.io/Interactive_CovidOnly_Heatmap/) |
| **COVID-19 + 13 Vax Atlas** | Compare COVID-19 with 13 other vaccines | [Launch App](https://wapsyed.shinyapps.io/Interactive_Covidand13Vax_Heatmap/) |

---

## Analysis Pipeline

The complete analysis is organized into **13 steps**, implemented in the main R Markdown notebook:

| # | Step | Description |
|---|------|-------------|
| 1 | **Preparing** | Package loading, color scheme standardization, figure themes |
| 2 | **Gene Sets Construction** | Curated immune gene sets: Blood Transcription Modules (BTM), CellMarker, ImmuneGO, VAXSigDB, TCR/BCR repertoires |
| 3 | **Study Selection & Data Retrieval** | Download and process raw count matrices from GEO (GSE189039, GSE199750, GSE201530, GSE201533, GSE206023) |
| 4 | **General Data Analysis** | Metadata annotation, sample unification, TPM/CPM normalization, biotype classification |
| 5 | **Differential Gene Expression** | DESeq2 Wald test per condition, FDR correction, contrast generation |
| 6 | **DEG Analysis** | Up/down regulation classification, immune vs non-immune gene assignment, Fisher-Jaccard overlaps |
| 7 | **Gene Set Enrichment** | GSEA and ORA against ImmuneGO, BTM, CellMarker, and OtherGO gene collections |
| 8 | **Principal Component Analysis** | PCA on immune gene L2FC matrix, variable contribution, KNN clustering |
| 9 | **Machine Learning** | Random Forest classifiers (tidymodels/ranger) to discriminate vaccine types, infection, and conditions; VIP gene identification |
| 10 | **Immune Gene Heatmaps** | ComplexHeatmaps of DEGs across all conditions, faceted by immune process, day, and vaccine type |
| 11 | **13-Vaccine Atlas Comparison** | Integration with Hagan et al. transcriptomic vaccine atlas — shared and unique immune gene responses across 13 vaccines |
| 12 | **MSigDB VAX Collection** | Enrichment analysis against the curated VAXSigDB gene set collection |
| 13 | **Gene Networks** | Network construction (Cytoscape, Gephi), chord diagrams, shared gene overlap visualization |

---

## Repository Structure

```
covidvax_atlas/
├── Article_Covid19VaxAtlas.Rmd        # Main analysis notebook (40k+ lines)
├── ML_Covid19VaxAtlas.Rmd             # Machine learning-focused notebook
├── 13VaxAtlas/                        # Hagan et al. vaccine atlas integration data
├── Annotations and metadata/          # Sample annotations, vaccine conditions, demographics, color palettes
├── Cell populations/                  # xCell deconvolution results
├── DGE analysis/                      # DEG tables, DESeq2 results, count matrices (raw/normalized)
├── Figures/                           # All publication-ready figures (~483 files)
├── Gene set enrichment analysis/      # GSEA results (ImmuneGO, BTM, CellMarker)
├── Interactive_CovidOnly_Heatmap/     # Shiny app for COVID-only gene explorer
├── Interactive_Covidand13Vax_Heatmap/ # Shiny app for COVID + 13 vaccines comparison
├── Networks/                          # Cytoscape/Gephi networks, edge/node tables, chord diagram files
├── Notebooks/                         # Supplementary notebooks
│   ├── VaxGO_Gene_sets_construction.Rmd
│   ├── Gene network.Rmd
│   ├── MachineLearning_Descriptors.Rmd
│   └── RandomEffects.Rmd
├── PCA/                               # PCA results, loadings, contribution plots
├── Tables/                            # Summary statistics tables (demographics, data curation)
├── VaxGO gene sets/                   # Curated immune gene set definitions (ImmuneGO, BTM, VAXSigDB, CellMarker)
└── raw_data/                          # Downloaded GEO datasets (GSE189039, GSE199750, GSE201530, GSE201533, GSE206023)
```

---

## How to Reproduce

### Prerequisites

- R >= 4.2
- Bioconductor packages: `DESeq2`, `clusterProfiler`, `ComplexHeatmap`, `GEOquery`, `biomaRt`, `org.Hs.eg.db`
- CRAN packages: `tidyverse`, `tidymodels`, `ranger`, `FactoMineR`, `factoextra`, `ggsci`, `janitor`, `here`, `edgeR`, etc. (see setup chunk in `Article_Covid19VaxAtlas.Rmd`)

### Steps

1. Clone the repository
2. Open `covidvax_atlas.Rproj` in RStudio
3. Run `Article_Covid19VaxAtlas.Rmd` (chunks are set to `eval=FALSE` by default; set `eval=TRUE` to execute)
4. For machine learning analyses, run `ML_Covid19VaxAtlas.Rmd`
5. For gene set construction, run `Notebooks/VaxGO_Gene_sets_construction.Rmd`
6. For network analyses, run `Notebooks/Gene network.Rmd`

> **Note:** Raw data must be downloaded from GEO before running. See Section 3 of the main notebook for download instructions.

---

## Citation

Prates-Syed WA, Fonseca DLM, Zaki Pour S, et al. **COVID-19 Vaccination Atlas: An Integrative Systems Vaccinology Approach.** *medRxiv* 2024.05.22.24307755v3.

```
https://doi.org/10.1101/2024.05.22.24307755v3
```

---

## License

This project is licensed under **CC0 1.0 Universal** — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Wasim Aluísio Prates-Syed** — wasim.syed@usp.br  
**Gustavo Cabral-Miranda** — gcabral.miranda@usp.br

Department of Immunology, Institute of Biomedical Sciences (ICB), University of São Paulo (USP), Brazil.
