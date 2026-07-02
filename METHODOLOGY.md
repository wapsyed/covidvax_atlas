# Methodology

## COVID-19 Vaccination Atlas: A Systems Vaccinology Analysis

---

## 1. Data Acquisition

Five publicly available bulk RNA-seq datasets were retrieved from the Gene Expression Omnibus (GEO):

| GEO Accession | Vaccine(s) | Platform | Samples |
|---------------|-----------|----------|---------|
| GSE189039 | BNT162b2, mRNA-1273, Ad26.COV2.S | Illumina NovaSeq 6000 | ~200 |
| GSE199750 | BBIBP-CorV, ZF2001 | Illumina NovaSeq 6000 | ~120 |
| GSE201530 | ChAdOx1, BNT162b2 | Illumina NovaSeq 6000 | ~90 |
| GSE201533 | BNT162b2 | Illumina NovaSeq 6000 | ~90 |
| GSE206023 | BBIBP-CorV, ZF2001 | Illumina NovaSeq 6000 | ~180 |

**Total:** 681 samples from 343 participants, including healthy vaccinated individuals, breakthrough infections, and unvaccinated infected controls.

---

## 2. Gene Set Construction (VaxGO)

A comprehensive immune gene set collection (named **VaxGO**) was built from multiple sources:

### 2.1 Blood Transcription Modules (BTM)
- 346 modules from Li et al. (Nature Immunology, 2014)
- Curated and filtered for immune-related processes

### 2.2 CellMarker Database
- Cell-type-specific markers from the CellMarker database (Zhang et al., Nucleic Acids Res., 2019)
- Filtered for blood/immune cell types
- Grouped into: B cells, T cells (CD4, CD8), NK cells, monocytes, macrophages, dendritic cells, neutrophils, eosinophils, mast cells

### 2.3 ImmuneGO
- Manually curated Gene Ontology-based immune process gene sets
- Structured hierarchically into:
  - **General processes:** Innate Immune System, Adaptive Immune System, Antimicrobial, Antiviral & IFN, Inflammation, Signaling, Complement, Immunoglobulin-mediated Response
  - **Specific cell types:** B cell, T cell, NK cell, Monocyte, Macrophage, Dendritic Cell, Neutrophil, Eosinophil, Mast Cell
- Annotated as "Manual" (curated) or by GO term enrichment

### 2.4 VAXSigDB
- Curated vaccine-specific gene signatures from published studies
- Filtered for PBMC samples, "VAC ONLY" study subtype

### 2.5 Non-Immune Gene Sets
- Metabolic processes, cellular processes (apoptosis, cell cycle, etc.)
- Used as negative controls and for non-immune characterization

### 2.6 TCR and BCR Repertoire Genes
- Separate gene lists for T cell receptor and B cell receptor repertoire components

---

## 3. Data Processing and Normalization

### 3.1 Count Matrix Preparation
- Raw count matrices downloaded from GEO supplementary files
- Gene symbols standardized using `biomaRt` (GRCh38/hg38)
- Sample metadata harmonized across studies (day, dose, vaccine type, infection status, severity)

### 3.2 Normalization
- **TPM (Transcripts Per Million):** Gene length-normalized for within-sample comparisons
- **CPM (Counts Per Million):** Library-size-normalized for between-sample comparisons
- **VST (Variance Stabilizing Transformation):** DESeq2 vst for PCA and heatmap visualizations
- `edgeR` TMM normalization for differential expression

### 3.3 Data Integration
- All five datasets merged into unified matrices:
  - `all_matrices_counts.rds` — raw counts
  - `all_matrices_counts_tpm_normalized.rds` — TPM-normalized
  - `all_matrices_cpm_normalized.rds` — CPM-normalized

---

## 4. Differential Gene Expression Analysis

### 4.1 DESeq2 Workflow
- Independent `DESeqDataSet` per study
- Design formula: `~ condition` (where condition encodes vaccine, dose, day, infection status)
- Wald test for pairwise contrasts
- Multiple testing correction via Benjamini-Hochberg (FDR < 0.05)

### 4.2 Weighted FDR Adjustment
To account for small sample sizes, an adaptive FDR threshold was applied:
- If `n >= 5`: `fdr <= 0.05`
- If `n < 5`: `fdr <= 0.05 * 6/n`

### 4.3 Direction Classification
- **UP:** log2 fold change >= 0.58 (1.5-fold)
- **DOWN:** log2 fold change <= -0.58
- **NEUTRAL:** between thresholds

### 4.4 Immune vs Non-Immune Classification
Each DEG was annotated as:
- **Immune:** if present in any ImmuneGO gene set
- **Non-immune:** if present in OtherGOs (metabolism, cellular processes)
- Further subclassified by BTM membership

---

## 5. Gene Set Enrichment Analysis

### 5.1 Over-Representation Analysis (ORA)
- ImmuneGO and BTM gene sets as reference
- Fisher's exact test with FDR correction
- Up- and down-regulated genes tested separately

### 5.2 Gene Set Enrichment Analysis (GSEA)
- Pre-ranked by log2 fold change
- Weighted enrichment statistic
- 10,000 permutations
- Performed against:
  - ImmuneGO general processes
  - ImmuneGO specific (cell-type) processes
  - BTM modules
  - CellMarker cell types

### 5.3 ssGSEA
- Single-sample GSEA for individual sample-level enrichment scores
- Used for correlation analysis with xCell cell-type abundance estimates

---

## 6. Principal Component Analysis

### 6.1 Input Data
- L2FC matrix: genes × conditions (mean log2 fold change per condition)
- Filtered for immune-related DEGs

### 6.2 PCA Workflow (FactoMineR)
- PCA on centered/scaled gene expression matrix
- Variable contribution analysis to identify key discriminating genes
- KNN clustering on PC coordinates to identify condition groupings

### 6.3 Batch/Study Effect Assessment
- PVCA (Principal Variance Component Analysis) to quantify study, vaccine type, day, and sex contributions
- Sex prediction from gene expression (XIST, RPS4Y1) for validation

---

## 7. Machine Learning

### 7.1 Classification Models
- **Random Forest:** `ranger` implementation in `tidymodels` framework
- **Additional algorithms tested:** Decision trees (rpart), logistic regression (glmnet), KNN (kknn), neural networks (nnet)

### 7.2 Workflow
1. **Data splitting:** 70% training, 30% testing (stratified by outcome)
2. **Preprocessing:** `step_nzv`, `step_corr`, `step_normalize`, `step_smote` (themis) for class imbalance
3. **Hyperparameter tuning:** grid search with 10-fold cross-validation
4. **Model evaluation:** accuracy, ROC-AUC, precision-recall, F1-score
5. **Feature importance:** VIP (Variable Importance Plot) scores

### 7.3 Outcomes Predicted
- Vaccine type (BBIBP vs ZF2001 vs BNT vs MO vs ChAd)
- Infection status (infected vs vaccinated only)
- Condition (specific vaccine + dose combinations)

---

## 8. Heatmap Visualization

### 8.1 ComplexHeatmap
- Hierarchically clustered heatmaps of DEG L2FC across all conditions
- Row annotations: immune process, cell type, BTM module
- Column annotations: vaccine type, day, dose, infection status

### 8.2 Interactive Heatmaps (Shiny)
Two interactive Shiny applications built with `InteractiveComplexHeatmap`:
- **COVID-19 Only:** Focused on COVID-19 vaccine and infection conditions
- **COVID-19 + 13 Vaccines:** Expanded comparison with 13 other vaccines from Hagan et al.

Features: gene search, regular expression filtering, sub-heatmap selection, annotation table export.

---

## 9. 13-Vaccine Atlas Comparison

### 9.1 Reference Data
- Integrated with Hagan et al. (2023) vaccine transcriptomic atlas
- 13 additional vaccines (influenza, yellow fever, smallpox, HBV, HPV, etc.)
- Shared condition schema (vaccine type, time points)

### 9.2 Comparative Analysis
- Shared and unique immune genes between COVID-19 and other vaccines
- Jaccard similarity and Fisher's exact test for gene set overlap
- PCA on combined COVID-19 + 13-vaccine dataset
- Separate analysis for vaccine-only and infection-derived signatures

---

## 10. Gene Networks

### 10.1 Network Construction
- Bipartite networks: genes connected to conditions (vaccine or infection)
- Immune and non-immune gene networks constructed separately
- Edge weight: log2 fold change direction and magnitude

### 10.2 Network Analysis Tools
- **Cytoscape:** Force-directed layout, cluster identification
- **Gephi:** Hierarchical clustering, modularity analysis
- **Visualization:** Chord diagrams (`chorddiag` R package) showing condition-gene-process relationships

### 10.3 Overlap Analysis
- Venn diagrams and UpSet plots for shared DEGs across conditions
- Jaccard distance matrices for hierarchical clustering of immune responses

---

## 11. Cell Population Deconvolution

### 11.1 xCell
- Cell-type enrichment scores estimated using xCell (Aran et al., 2017)
- 64 immune and stromal cell types
- ssGSEA-based scoring normalized to a reference

### 11.2 Downstream Analysis
- Correlation of cell-type scores with immune gene expression
- Line charts of cell population dynamics over time post-vaccination
- Comparison of cell-type profiles across vaccine platforms

---

## 12. Random Effects Modeling

To account for repeated measures and individual-level variation:
- Mixed-effects models using `lme4` and `glmmTMB`
- Random intercept per participant
- Fixed effects: time (day), vaccine type, dose
- Applied to GSE206023 (BBIBP/ZF2001 longitudinal data)

---

## Main Files Reference

### Core Analysis Scripts

| File | Description | Lines |
|------|-------------|-------|
| `Article_Covid19VaxAtlas.Rmd` | Main analysis notebook: complete pipeline (sections 1–13), fully commented in English | ~40,500 |
| `ML_Covid19VaxAtlas.Rmd` | Machine learning-focused analysis (RF, hyperparameter tuning, VIP genes) | ~5,600 |

### Code Documentation

All code in `Article_Covid19VaxAtlas.Rmd` is fully documented with:
- **Markdown headers** (`##`, `###`, `####`) describing each section and subsection
- **Explanatory text** before each code chunk describing the analysis purpose
- **Inline comments** (`#`) for key operations within code chunks
- **Plot annotations:** `# Save plot`, `# Print plot` for all graphic outputs
- **Heatmap annotations:** `# Save heatmap to PNG`, `# Build heatmap object`, `# Print heatmap`
- **Data save annotations:** `# Save data` before file writes

### Supplementary Notebooks

| File | Description |
|------|-------------|
| `Notebooks/VaxGO_Gene_sets_construction.Rmd` | BTM, CellMarker, ImmuneGO gene set curation |
| `Notebooks/Gene network.Rmd` | Cytoscape/Gephi network construction and chord diagrams |
| `Notebooks/MachineLearning_Descriptors.Rmd` | Additional ML descriptor analysis |
| `Notebooks/RandomEffects.Rmd` | Mixed-effects modeling (lme4, glmmTMB) |
| `Annotations and metadata/Vaccine landscape/Vaccine landscape.Rmd` | COVID-19 vaccine landscape overview |
| `Annotations and metadata/Colors_paper.R` | Publication color palette definitions |

### Key Data Files

| File | Contents |
|------|----------|
| `DGE analysis/all_studies_degs_weighted.rds` | Master DEG table with weighted FDR correction |
| `DGE analysis/all_matrices_counts.rds` | Unified raw count matrix (all 5 studies) |
| `DGE analysis/all_matrices_counts_tpm_normalized.rds` | TPM-normalized expression matrix |
| `DGE analysis/all_matrices_cpm_normalized.rds` | CPM-normalized expression matrix |
| `VaxGO gene sets/ImmuneGO_Annotated_Genes_2024-03-26.rds` | Curated immune gene-to-process annotations |
| `VaxGO gene sets/ImmuneGO_Annotated_GO_2024_03_26.csv` | ImmuneGO ontology table with process groupings |
| `VaxGO gene sets/btm_annotation_immune.rds` | BTM module immune annotations |
| `VaxGO gene sets/CellMarker_ImmuneCells.csv` | Immune cell-type markers from CellMarker |
| `Annotations and metadata/ann_vaccines_conditions.csv` | Condition-level annotations (vaccine, dose, day) |
| `Annotations and metadata/ann_vaccines_samples.csv` | Sample-level metadata (age, sex, condition) |
| `13VaxAtlas/IS2_FDR05_genes.rds` | Hagan et al. 13-vaccine atlas DEGs |
| `Gene set enrichment analysis/all_degs_p_05_vac_infected_GSEA_ALL_8-1-24.csv` | Complete GSEA results |

### Interactive Applications

| File | Description |
|------|-------------|
| `Interactive_CovidOnly_Heatmap/app.R` | Shiny app: COVID-only immune gene expression explorer |
| `Interactive_Covidand13Vax_Heatmap/app.R` | Shiny app: COVID-19 + 13 vaccine atlas comparison |
| `Interactive_CovidOnly_Heatmap/htShiny.RData` | Precomputed heatmap objects for the COVID-only app |
| `Interactive_Covidand13Vax_Heatmap/htShiny.RData` | Precomputed heatmap objects for the multi-vaccine app |
