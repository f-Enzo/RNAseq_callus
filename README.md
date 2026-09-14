# RNA-seq Analysis of Grapevine Anthocyanin Biosynthesis

RNA-seq analysis scripts for investigating anthocyanin biosynthesis pathways in grapevine (*Vitis vinifera*) cell cultures across multiple cultivars.

## About This Project

Although anthocyanin biosynthesis is well characterized in grapevine, key questions remain regarding the **spatio-temporal regulation of the anthocyanin pathway** and the molecular mechanisms underlying the teinturier phenotype (red-fleshed cultivars with pigmented flesh tissue).

This project analyzes RNA-seq data from five grapevine cultivars to uncover the genomic and transcriptional mechanisms controlling anthocyanin production:

- **Gamay Fréaux** — black-skinned, red-fleshed teinturier (pigment accumulation in tissues)
- **Merlot**, **Mourvèdre**, **Syrah** — black-skinned cultivars (no red flesh)
- **Muscat Blanc** — white-skinned cultivar (control)

The analysis employs **whole-genome and RNA sequencing** to identify differentially expressed genes, enriched metabolic pathways, and genetic differences in anthocyanin regulation.

## Scripts

### `RNAseq_analysis.R`

**Purpose:** Multi-sample global analysis of anthocyanin pathway gene expression

**Inputs:**
- `callus_metadata.csv` — sample metadata (cultivar, phenotype annotations)
- `featureCount_unstranded.tsv` — raw gene counts from featureCounts
- `ch_pn_genes_extended.csv` — gene annotation file (Vitis vinifera gene IDs mapped to NCBI/KEGG identifiers)

**Key Analyses:**

| Analysis | Method | Output |
|----------|--------|--------|
| **Normalization** | TMM (Trimmed Mean of M-values) via edgeR | log-CPM normalized counts |
| **Multivariate Analysis** | PCA + PLS-DA via mixOmics | Sample separation, gene-component correlations |
| **Differential Expression** | DESeq2 (Wald test + LRT) | Significantly expressed genes, MA plots, heatmaps |
| **Pathway Enrichment** | KEGG + GO (clusterProfiler) | Anthocyanin pathway (vvi00942), Flavonoid pathway (vvi00941) |
| **Gene Set Enrichment** | FGSEA (fgsea) | Ranked pathway analysis with leading edge genes |

**Main Outputs:**
- PCA/PLS-DA plots of sample clustering
- Heatmaps of top differentially expressed genes
- KEGG pathway enrichment dot plots
- Gene ontology enrichment visualizations
- Identified genes highly correlated with PCA/PLS-DA components

---

### `pairwise_rnaseq_analysis.R`

**Purpose:** Pairwise differential expression analysis between Gamay Fréaux and other cultivars

**Inputs:**
- Same as above (metadata, count table, gene annotations)

**Key Analyses:**

| Analysis | Method | Output |
|----------|--------|--------|
| **Pairwise DE** | DESeq2 (Gamay Fréaux vs. each cultivar) | log2 fold-changes and p-values for each comparison |
| **Gene Filtering** | Significance thresholds | Genes with padj < 0.01 AND \|log2FC\| > 2 |
| **Pathway Visualization** | pathview | KEGG anthocyanin (vvi00941) and flavonoid (vvi00942) pathway maps annotated with expression data |

**Main Outputs:**
- Aggregated results table with log2FC and padj for each comparison
- Heatmap of significant genes across comparisons
- Annotated KEGG pathway visualizations (one per cultivar pair)
- Ranked gene lists for each comparison (ordered by log2FC)

---

## Key Features

✓ **Multi-method normalization** — TMM normalization for robust count-based analysis  
✓ **Dimensionality reduction** — PCA and supervised PLS-DA for sample separation  
✓ **Differential expression** — Wald test for condition effects, LRT for cultivar effects  
✓ **Functional annotation** — KEGG pathway and GO enrichment analysis  
✓ **Metabolic pathway mapping** — Anthocyanin and flavonoid pathway visualization  
✓ **Pairwise comparisons** — Systematic cultivar-to-cultivar expression differences  

## Required R Packages

- **Data manipulation:** tidyverse, readr, data.table
- **Differential expression:** DESeq2, edgeR
- **Multivariate analysis:** mixOmics, FactoMineR
- **Visualization:** pheatmap, corrplot, ggplot2
- **Functional annotation:** clusterProfiler, fgsea, enrichplot, pathview
- **Genome databases:** org.Vvinifera.eg.db

## Getting Started

1. **Prepare input files** in your working directory:
   - `callus_metadata.csv` (sample metadata)
   - `featureCount_unstranded.tsv` (raw gene counts)
   - `ch_pn_genes_extended.csv` (gene annotation mapping)

2. **Update file paths** in the scripts:
   - Modify `setwd()` to your project directory
   - Adjust input file paths as needed

3. **Run the analyses:**
   ```r
   # Global multi-sample analysis
   source("RNAseq_analysis.R")
   
   # Pairwise cultivar comparisons
   source("pairwise_rnaseq_analysis.R")
   ```

4. **Outputs** will be generated in your working directory, including:
   - PNG/PDF plots (PCA, heatmaps, pathway maps)
   - Console-printed result tables

## Methods Summary

- **RNA-seq processing:** featureCounts for gene quantification
- **Normalization:** TMM normalization (edgeR::calcNormFactors)
- **Statistical testing:** DESeq2 for differential expression
- **FDR control:** Benjamini-Hochberg adjustment (padj < 0.05)
- **Pathway databases:** KEGG (Kyoto Encyclopedia of Genes and Genomes)
- **Organism:** *Vitis vinifera* (grapevine)

## Related Publications

- Kőrösi et al. (2022) — Genomic basis of teinturier phenotype
- Röckel et al. (2020) — Anthocyanin biosynthesis characterization

---

