# Pairwise RNA-seq data analysis of 5 grapewine callus 

# Load packages
library(tidyverse)
library(readr)
library(data.table)
library(DESeq2)
library(mixOmics)
library(FactoMineR)
library(edgeR)
library(corrplot)
library(pheatmap)
library(clusterProfiler)
library(fgsea)
library(org.Vvinifera.eg.db)
library(enrichplot)
library(pathview)

setwd("/home/enzo/Documents/rnaseq_callus/")

# Load data
metadata <- read_csv("callus_metadata.csv")
count_table <- read_delim("featureCount_unstranded.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
ch_pn_genes <- read_csv("/home/enzo/Documents/GeneFamilies/ch_pn_genes_extended.csv")

# Differential expression
# Load data in a deseq2 object
X <- count_table %>%
  column_to_rownames("Geneid") %>%  # Set Geneid as row names
  dplyr::select(-Length) %>%  # Remove Length column
  filter_if(is.numeric, any_vars(. != 0)) %>%  # Remove genes with only zero counts
  t() %>%  # Transpose
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%  # Store rownames in a column
  mutate(Sample = sub("_.*", "", Sample)) %>%
  arrange(parse_number(Sample)) %>%  # Sort numerically by extracted number
  remove_rownames %>%
  column_to_rownames("Sample") # Drop sample columns

dds <- DESeqDataSetFromMatrix(countData = t(X),
                              colData = metadata,
                              design = ~ cultivar)

dds

# Filter very low expressed genes
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Pairwise DE comparison
dds$cultivar <- as.factor(dds$cultivar)

# Get all cultivars except "Gamay Fréau"
cultivars <- levels(dds$cultivar)
cultivars <- cultivars[cultivars != "Gamay Fréaux"]

# Create an empty dataframe to store results
pairwise_res <- data.frame(row.names = rownames(dds))
# Compute average expression for Gamay Fréau samples
pairwise_res$avgExpr_Gamay <- rowMeans(counts(dds, normalized = FALSE)[, colData(dds)$cultivar == "Gamay Fréaux", drop = FALSE])


# Loop through each cultivar and perform differential analysis
for (cultivar in cultivars) {
  # Create a subset of the data for the current contrast
  dds_subset <- dds[, dds$cultivar %in% c(cultivar, "Gamay Fréaux")]
  dds_subset$cultivar <- droplevels(dds_subset$cultivar)  # Drop unused levels
  
  # Re-run DESeq2 for this specific subset
  dds_subset <- DESeq(dds_subset, minReplicatesForReplace=Inf)
  res <- results(dds_subset, contrast = c("cultivar", "Gamay Fréaux", cultivar), cooksCutoff=FALSE, independentFiltering=FALSE) # add cooksCutoff=FALSE to not get NAs when DESeq2 detects an outlier
  
  # Extract log2FoldChange and padj with a suffix for the comparison
  pairwise_res[[paste0("log2FC_", cultivar)]] <- res$log2FoldChange
  pairwise_res[[paste0("padj_", cultivar)]] <- res$padj
}

# Compute max absolute log2 fold change across all comparisons
lfc_columns <- grep("log2FC_", colnames(pairwise_res), value = TRUE)
pairwise_res$maxLFC <- apply(pairwise_res[, lfc_columns], 1, function(x) {
  idx <- which.max(abs(x))  # Get index of max absolute LFC
  return(x[idx])            # Return the real (signed) value
})

# Compute minimum adjusted p-value across all comparisons
pval_columns <- grep("padj_", colnames(pairwise_res), value = TRUE)
pairwise_res$minPvalue <- apply(pairwise_res[, pval_columns], 1, function(x) min(x, na.rm = TRUE))

# View final aggregated results
head(pairwise_res)

pairwise_res_ordered <- pairwise_res[order(pairwise_res$maxLFC, decreasing = TRUE),]

sig_genes <- rownames(pairwise_res[which(pairwise_res$minPvalue < 0.01 & pairwise_res$maxLFC > 2),])

# Visualise log2FC for significant DE genes
lfc_mat <- as.matrix(pairwise_res[, lfc_columns])

gene_annotation <- ch_pn_genes %>% 
  dplyr::filter(chasselas %in% sig_genes) %>%
  arrange(match(chasselas, sig_genes)) %>%
  dplyr::select(chasselas, kegg_ko) %>%
  column_to_rownames("chasselas")

x11()
pheatmap(lfc_mat[sig_genes, ], show_rownames=FALSE, annotation_row = gene_annotation) 

# Visualize anthocyanin pathway for each comparison 

# Translate Vv to NCBI ids
translate_vv_to_id <- function(de_res, translation_df, id_column, lfc_column) {
  
  # Make a data frame with the Vv genes ordered by log2FC and there corresponding id in the other annotation
  ranked_genes_df <- data.frame(chasselas = rownames(de_res), lfc = de_res[[lfc_column]], stringsAsFactors = FALSE) %>%
    left_join(ch_pn_genes %>% 
                mutate(ncbi_id = sub("^LOC", "", ncbi_id)),
              by = "chasselas") %>%
    dplyr::filter(!is.na(!!sym(id_column))) %>%
    arrange(desc(lfc))
  
  if (id_column == 'ncbi_id') {
    ranked_genes_df <- ranked_genes_df%>%
      arrange(desc(abs(lfc))) %>%
      distinct(ncbi_id, .keep_all = TRUE) %>% # remove duplicated ncbi locus, keeping the ones with highest lfc (to avoid error in gseKEGG)
      arrange(desc(lfc))
  }
  
  ranked_gene_list <- pull(ranked_genes_df, !!sym(id_column)) 
  ranked_gene_nv <- setNames(ranked_genes_df$lfc, ranked_gene_list)
  
  return(list(names=ranked_gene_list, nv=ranked_gene_nv))
}

for (lfc_col in lfc_columns) {
  
  # List of genes ranked by lfc
  ranked_genes_vv <- rownames(pairwise_res[order(pairwise_res[[lfc_col]], decreasing = T),])
  
  ranked_genes_vv_nv <- setNames(as.numeric(pairwise_res[order(pairwise_res[[lfc_col]], decreasing = T),][[lfc_col]]), ranked_genes_vv)
  
  # NCBI ids
  ranked_genes_ncbi <- translate_vv_to_id(pairwise_res, ch_pn_genes, 'ncbi_id', lfc_col)
  
  pathview(gene.data = ranked_genes_ncbi$nv, 
           pathway.id = "vvi00941", 
           species = "vvi", 
           out.suffix = lfc_col,
           kegg.dir = "/home/enzo/Documents/rnaseq_callus",
           limit = list(gene=10, cpd=2))
}


ranked_genes_df <- data.frame(chasselas = rownames(pairwise_res), lfc = pairwise_res$maxLFC, stringsAsFactors = FALSE) %>%
  left_join(ch_pn_genes %>% 
              mutate(ncbi_id = sub("^LOC", "", ncbi_id)),
            by = "chasselas") %>%
  dplyr::filter(!is.na(ncbi_id)) %>%
  arrange(desc(abs(lfc))) %>%
  distinct(ncbi_id, .keep_all = TRUE) %>% # remove duplicated ncbi locus, keeping the ones with highest lfc (to avoid error in gseKEGG)
  arrange(desc(lfc))


