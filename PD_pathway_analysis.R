# ============================================================
# Parkinson's Disease – Differential Expression & Pathway Analysis
# Dataset: GSE68719 (GEO) – Substantia nigra, PD vs Control
# Tools: DESeq2, clusterProfiler, KEGG, GO, GSEA
# Author: Hammad Ali Hassan
# Date: 2025
# ============================================================

# ── 1. INSTALL & LOAD PACKAGES ───────────────────────────────
# Run this block once to install all required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "AnnotationDbi"
))

install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer", "ggrepel"))

# Load libraries
library(GEOquery)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)


# ── 2. DOWNLOAD DATA FROM GEO ────────────────────────────────
# GSE68719: RNA-seq of substantia nigra from PD patients vs healthy controls
gse <- getGEO("GSE68719", GSEMatrix = TRUE, getGPL = FALSE)
gse <- gse[[1]]

# View sample metadata
pData(gse)[, c("title", "characteristics_ch1", "characteristics_ch1.1")]


# ── 3. PREPARE COUNT MATRIX & METADATA ───────────────────────
# Note: GSE68719 provides processed expression values.
# We use the expression matrix as a proxy count matrix for DESeq2.
# For a fully rigorous analysis, raw FASTQ files would be used.

count_matrix <- exprs(gse)

# Round to integers (DESeq2 requires integer counts)
count_matrix <- round(count_matrix)
count_matrix[count_matrix < 0] <- 0

# Build sample metadata
sample_info <- pData(gse)
sample_info$condition <- ifelse(
  grepl("Parkinson", sample_info$characteristics_ch1, ignore.case = TRUE),
  "PD", "Control"
)
sample_info$condition <- factor(sample_info$condition, levels = c("Control", "PD"))

# Confirm group sizes
table(sample_info$condition)


# ── 4. RUN DESeq2 ─────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ condition
)

# Filter low-count genes (keep genes with >= 10 counts in at least 3 samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

# Run DESeq2 differential expression
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "PD", "Control"))
res <- res[order(res$padj), ]

# Summary
summary(res)

# Save full results
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
write.csv(res_df, "DESeq2_PD_results.csv", row.names = FALSE)


# ── 5. VOLCANO PLOT ───────────────────────────────────────────
res_df$significance <- "Not Significant"
res_df$significance[res_df$log2FoldChange > 1  & res_df$padj < 0.05] <- "Upregulated"
res_df$significance[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Downregulated"

# Label top 15 significant genes
top_genes <- res_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(15)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated"     = "#E64B35",
    "Downregulated"   = "#4DBBD5",
    "Not Significant" = "grey70"
  )) +
  geom_text_repel(
    data    = top_genes,
    aes(label = gene),
    size    = 3,
    max.overlaps = 20
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
  labs(
    title    = "Differential Gene Expression in Parkinson's Disease",
    subtitle = "Substantia nigra: PD vs Control (GSE68719)",
    x        = "Log2 Fold Change",
    y        = "-Log10 Adjusted P-value",
    color    = ""
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("volcano_plot_PD.png", volcano_plot, width = 8, height = 6, dpi = 300)
print(volcano_plot)


# ── 6. GENE ID CONVERSION (Symbol → Entrez) ───────────────────
# Required for clusterProfiler KEGG and GO analyses
sig_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

entrez_ids <- bitr(
  sig_genes$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

sig_genes <- merge(sig_genes, entrez_ids, by.x = "gene", by.y = "SYMBOL")


# ── 7. KEGG PATHWAY ENRICHMENT ────────────────────────────────
kegg_res <- enrichKEGG(
  gene         = sig_genes$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_plot <- dotplot(kegg_res, showCategory = 15) +
  labs(
    title    = "KEGG Pathway Enrichment – Parkinson's Disease",
    subtitle = "Significantly dysregulated pathways (padj < 0.05)"
  ) +
  theme_classic(base_size = 12)

ggsave("KEGG_dotplot_PD.png", kegg_plot, width = 9, height = 7, dpi = 300)
print(kegg_plot)


# ── 8. GO ENRICHMENT (Biological Process) ─────────────────────
go_res <- enrichGO(
  gene          = sig_genes$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",           # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

go_plot <- dotplot(go_res, showCategory = 15) +
  labs(
    title    = "GO Biological Process Enrichment – Parkinson's Disease",
    subtitle = "Top enriched biological processes"
  ) +
  theme_classic(base_size = 12)

ggsave("GO_dotplot_PD.png", go_plot, width = 10, height = 7, dpi = 300)
print(go_plot)


# ── 9. GSEA (Gene Set Enrichment Analysis) ────────────────────
# GSEA uses ALL genes ranked by log2FoldChange — not just significant ones
all_entrez <- bitr(
  res_df$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

res_ranked <- merge(res_df, all_entrez, by.x = "gene", by.y = "SYMBOL")
res_ranked <- res_ranked[!is.na(res_ranked$log2FoldChange), ]

# Create named ranked vector
gene_list <- res_ranked$log2FoldChange
names(gene_list) <- res_ranked$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

# Run GSEA on KEGG pathways
gsea_res <- gseKEGG(
  geneList     = gene_list,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# GSEA enrichment plot for top 3 pathways
gsea_plot <- gseaplot2(
  gsea_res,
  geneSetID = 1:3,
  title     = "GSEA – Top Enriched KEGG Pathways in Parkinson's Disease"
)

ggsave("GSEA_plot_PD.png", gsea_plot, width = 10, height = 7, dpi = 300)
print(gsea_plot)


# ── 10. HEATMAP OF TOP 30 DEGs ────────────────────────────────
top30 <- res_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene)

# Variance-stabilised counts for visualisation
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top30, ]
mat <- mat - rowMeans(mat)  # center

annotation_col <- data.frame(
  Condition = sample_info$condition,
  row.names = colnames(mat)
)

ann_colors <- list(
  Condition = c(Control = "#4DBBD5", PD = "#E64B35")
)

png("heatmap_top30_PD.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(
  mat,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 9,
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main              = "Top 30 Differentially Expressed Genes – PD vs Control"
)
dev.off()


# ── 11. SESSION INFO ──────────────────────────────────────────
# Always include for reproducibility
sessionInfo()

# ============================================================
# OUTPUT FILES GENERATED:
#   DESeq2_PD_results.csv     – Full differential expression table
#   volcano_plot_PD.png       – Volcano plot
#   KEGG_dotplot_PD.png       – KEGG pathway enrichment dotplot
#   GO_dotplot_PD.png         – GO biological process dotplot
#   GSEA_plot_PD.png          – GSEA enrichment plot
#   heatmap_top30_PD.png      – Heatmap of top 30 DEGs
# ============================================================
