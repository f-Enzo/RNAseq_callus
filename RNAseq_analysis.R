# RNA-seq data analysis of 5 grapewine callus 

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

# Load data
metadata <- read_csv("Documents/rnaseq_callus/callus_metadata.csv")
count_table <- read_delim("Documents/rnaseq_callus/featureCount_unstranded.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
ch_pn_genes <- read_csv("Documents/GeneFamilies/ch_pn_genes_extended.csv")

# Multivariate analysis 
X <- count_table %>%
  column_to_rownames("Geneid") %>%  # Set Geneid as row names
  dplyr::select(-Length) %>%  # Remove Length column
  filter_if(is.numeric, any_vars(. != 0)) %>%  # Remove genes with only zero counts
  t() %>%  # Transpose
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%  # Store rownames in a column
  arrange(parse_number(Sample)) %>%  # Sort numerically by extracted number
  dplyr::select(-Sample) #%>% # Drop sample columns
  
# TMM Normalization using edgeR
dge <- DGEList(counts = as.matrix(X))  # Convert to DGEList
dge <- calcNormFactors(dge, method = "TMM")  # Apply TMM normalization
X_norm <- cpm(dge, log = TRUE)  # Obtain log-transformed normalized counts

# PCA With mixOmics
pca_res = pca(X_norm, ncomp = 5)

plot(pca_res)
plotIndiv(pca_res, group = metadata$cultivar, legend = TRUE, comp = c(1,2))
plotVar(pca_res, cutoff = 0.99)
biplot(pca_res, group = metadata$cultivar, cutoff = 0.99, legend.title = "Cultivar")
plotLoadings(pca_res, comp = 1, ndisplay = 20)
plotLoadings(pca_res, comp = 2, ndisplay = 20)


# Get the list of most correlated genes
# Compute the correlations between each gene (variable) and the PCA components.
var_cor <- cor(X_norm, pca_res$variates$X)

# Set the cutoff value
cutoff <- 0.99

# For each component, extract the genes with an absolute correlation greater than the cutoff.
genes_PC1 <- rownames(var_cor)[abs(var_cor[, 1]) > cutoff]
genes_PC2 <- rownames(var_cor)[abs(var_cor[, 2]) > cutoff]
genes_PC3 <- rownames(var_cor)[abs(var_cor[, 3]) > cutoff]
genes_PC4 <- rownames(var_cor)[abs(var_cor[, 4]) > cutoff]
genes_PC5 <- rownames(var_cor)[abs(var_cor[, 5]) > cutoff]


# Optionally, store the results in a list
correlated_genes <- list(PC1 = genes_PC1, PC2 = genes_PC2)
top_correlated <- c(genes_PC1, genes_PC2, genes_PC3, genes_PC4, genes_PC5)

# Display the genes for each component
print(correlated_genes)
corrplot(var_cor[top_correlated,])

# Differential expression
# Load data in a deseq2 object
dds <- DESeqDataSetFromMatrix(countData = t(X),
                              colData = metadata,
                              design = ~ Remarque)

dds

# Filter very low expressed genes
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Perform the DE analsis
dds <- DESeq(dds)
res <- results(dds)
res

# LFC shrinkage (see vignette) (remove the noise associated with log2 fold changes from low count genes)
resLFC <- lfcShrink(dds, coef="Remarque_red.mutation_vs_none", type="apeglm")
resLFC

summary(resLFC)
plotMA(resLFC, ylim=c(-6,6))

# Gene counts of the gene with lowest padj
plotCounts(dds, gene=which.min(res$padj), intgroup="Remarque")

# Order by lowest adjusted p value 
resOrdered <- res[order(res$padj),]

# Heat map of the count matrix
ntd <- normTransform(dds)

# Select mostly expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

# Select most significantly DE genes
select <- rownames(resOrdered)[1:20]

# Select top variance genes
topVarGenes <- head(order(rowVars(assay(ntd)), decreasing = TRUE), 20)
  
df <- as.data.frame(colData(dds)[,c("Remarque","cultivar")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(ntd)[topVarGenes,], annotation_col = df)

# Produce list of DE genes 
resSig <- subset(resOrdered, padj < 0.05)
rownames(resSig)

#Gene set enrichment analysis

# Plot ranked fold changes
barplot(sort(res$log2FoldChange, decreasing = T))
barplot(sort(resLFC$log2FoldChange, decreasing = T)) # Shrinked ones 

# List of genes ranked by lfc
ranked_genes <- rownames(res[order(res$log2FoldChange, decreasing = T),])

ranked_genes_nv <- setNames(as.numeric(res[order(res$log2FoldChange, decreasing = T),]$log2FoldChange), ranked_genes)

# Translate Vv to NCBI ids
ranked_genes_df <- data.frame(chasselas = rownames(res), lfc = res$log2FoldChange, stringsAsFactors = FALSE) %>%
  left_join(ch_pn_genes %>% 
              mutate(ncbi_id = sub("^LOC", "", ncbi_id)),
            by = "chasselas") %>%
  dplyr::filter(!is.na(ncbi_id)) %>%
  arrange(desc(abs(lfc))) %>%
  distinct(ncbi_id, .keep_all = TRUE) %>% # remove duplicated ncbi locus, keeping the ones with highest lfc (to avoid error in gseKEGG)
  arrange(desc(lfc))

ranked_genes_ncbi <- pull(ranked_genes_df, ncbi_id) 
ranked_genes_ncbi_nv <- setNames(ranked_genes_df$lfc, ranked_genes_ncbi)

# KEGG pathway enrichment analysis
top_genes <- names(ranked_genes_ncbi_nv)[abs(ranked_genes_ncbi_nv) > 2]

kk <- enrichKEGG(gene         = top_genes,
                 organism     = 'vvi',
                 pvalueCutoff = 0.05)

kk2 <- gseKEGG(geneList     = ranked_genes_ncbi_nv,
               organism     = 'vvi',
               minGSSize    = 10,
               pvalueCutoff = 0.1,
               verbose      = FALSE)
head(kk2)
# Visualise flavonoid pathway
browseKEGG(kk, 'vvi00941')

pathview(gene.data = ranked_genes_ncbi_nv, 
         pathway.id = "vvi00941", 
         species = "vvi", 
         limit = list(gene=5, cpd=1))

# Visualise anthocyanin pathway 
browseKEGG(kk, 'vvi00942')

# Visualise top enriched pathways
dotplot(kk, showCategory=20)

# Convert kk2 to a data.table and rename columns to match fgsea conventions:
fgseaRes <- as.data.table(as.data.frame(kk2))
setnames(fgseaRes,
         old = c("ID", "Description", "pvalue", "p.adjust", "enrichmentScore", "NES", "setSize", "core_enrichment"),
         new = c("ID", "pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge"))

# Build the pathways list:
# Here, for each pathway, split the 'leadingEdge' (core_enrichment) string by "/"
pathways <- lapply(fgseaRes$leadingEdge, function(x) unlist(strsplit(x, "/")))
names(pathways) <- fgseaRes$pathway

# Plot the GSEA table:
plotGseaTable(pathways, ranked_genes_ncbi_nv, fgseaRes, gseaParam = 1)


# GO enrichment analysis
ego <- enrichGO(gene         = ranked_genes,
                OrgDb        = org.Vvinifera.eg.db,
                keyType      = "GID",      # key type in OrgDb
                ont          = "BP",       # or "CC" / "MF" 
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

goplot(ego)
upsetplot(ego)
dotplot(ego)

