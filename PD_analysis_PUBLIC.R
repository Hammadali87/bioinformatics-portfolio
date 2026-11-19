# ============================================================
# Parkinson's Disease transcriptomics – differential expression & pathway analysis
# A demonstration of an RNA-seq analysis workflow on PUBLIC data
# Dataset : GSE68719 (NCBI GEO) – post-mortem BA9 prefrontal cortex,
#           Parkinson's disease (n = 29) vs control (n = 44)
# Input   : GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz (DESeq2-normalized counts,
#           protein-coding genes; GEO supplementary file, downloaded automatically)
# Methods : limma on log2(normalized count + 1) (the supplied values are normalized,
#           not raw counts, so DESeq2 is not appropriate here); PCA; KEGG and GO
#           over-representation and KEGG GSEA (clusterProfiler)
# Observed: strong immune/inflammatory enrichment – cytokine signalling, lymphocyte
#           activation, complement – consistent with neuroinflammation in PD. No
#           metabolic pathway stood out. Fold changes correlate with the original
#           authors' results (see results_summary.txt for the exact value).
# Limits  : no adjustment for age, post-mortem interval or RNA quality (the original
#           authors adjusted for these); bulk tissue, so cell-type composition can
#           contribute to immune signals.
# Note    : this is a workflow demonstration. Client analyses are confidential and
#           are not shown here.
# Author  : Hammad Ali Hassan (ORCID 0000-0002-9174-7673)
#
# Setup (first run only):
#   install.packages(c("BiocManager", "RSQLite", "ggplot2", "dplyr", "tidyr",
#                      "pheatmap", "RColorBrewer", "ggrepel"))
#   BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db",
#                          "enrichplot", "AnnotationDbi"))
# Usage: set the working directory to this script's folder, then source() the script.
#   Outputs are written to ./results, downloaded data to ./data. KEGG needs internet.
# ============================================================

DATA_DIR <- "data"
OUT_DIR  <- "results"

# ── 1. LOAD LIBRARIES ────────────────────────────────────────
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
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
})

LFC_CUT <- 0.5    # |log2 fold change| threshold for plots and enrichment
FDR_CUT <- 0.05   # adjusted p-value threshold

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(OUT_DIR, f)


# ── 2. DOWNLOAD & READ THE EXPRESSION FILES ──────────────────
# The GEO series-matrix file is empty for this dataset; the expression values
# are in supplementary files, which are downloaded once into DATA_DIR.
norm_file <- file.path(DATA_DIR, "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")
auth_file <- file.path(DATA_DIR, "GSE68719_mlpd_DESeq2_diffexp_pmi_age_rin_default.txt.gz")
if (!file.exists(norm_file)) {
  getGEOSuppFiles("GSE68719", baseDir = DATA_DIR, makeDirectory = FALSE)
}
norm <- read.delim(norm_file, check.names = FALSE)
cat("Genes x columns:", dim(norm), "\n")


# ── 3. BUILD EXPRESSION MATRIX & SAMPLE LABELS ───────────────
sample_cols <- setdiff(colnames(norm), c("EnsemblID", "symbol"))

# Several Ensembl IDs can share one gene symbol: keep the highest-expressed one
norm$mean_expr <- rowMeans(norm[, sample_cols])
norm <- norm[order(-norm$mean_expr), ]
norm <- norm[!duplicated(norm$symbol) & norm$symbol != "" & !is.na(norm$symbol), ]
expr <- as.matrix(norm[, sample_cols])
rownames(expr) <- norm$symbol

# Labels come from the sample IDs: C_ = control, P_ = Parkinson's disease
sample_info <- data.frame(
  sample    = sample_cols,
  condition = factor(ifelse(grepl("^P_", sample_cols), "PD", "Control"),
                     levels = c("Control", "PD")),
  row.names = sample_cols
)
print(table(sample_info$condition))
stopifnot(sum(sample_info$condition == "PD") == 29,
          sum(sample_info$condition == "Control") == 44)

# Data check: values are normalized (non-integer), so limma on log2 is appropriate
cat("All integers?", all(expr == round(expr)), "\n")
cat("Value range:", range(expr), "\n")
log_expr <- log2(expr + 1)


# ── 4. DIFFERENTIAL EXPRESSION WITH limma ────────────────────
keep <- rowSums(log_expr > 1) >= 10          # expressed in at least 10 samples
log_expr <- log_expr[keep, ]
cat("Genes kept after filtering:", nrow(log_expr), "\n")

design <- model.matrix(~ condition, data = sample_info)
fit <- lmFit(log_expr, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
tt  <- topTable(fit, coef = "conditionPD", number = Inf, sort.by = "P")

res_df <- data.frame(
  gene           = rownames(tt),
  log2FoldChange = tt$logFC,
  pvalue         = tt$P.Value,
  padj           = tt$adj.P.Val
)
rownames(res_df) <- res_df$gene
print(summary(decideTests(fit)))
write.csv(res_df, out("limma_PD_results.csv"), row.names = FALSE)

# Cross-check against the original authors' DESeq2 results
# (they adjusted for post-mortem interval, age and RNA quality; this analysis does not)
auth_cor <- NA
if (file.exists(auth_file)) {
  auth <- read.delim(auth_file)
  auth <- auth[!duplicated(auth$symbol), ]
  common <- intersect(res_df$gene, auth$symbol)
  auth_cor <- cor(res_df[common, "log2FoldChange"],
                  auth$log2FoldChange[match(common, auth$symbol)],
                  use = "complete.obs")
  cat("Genes in common with authors' results:", length(common), "\n")
  cat("Correlation of log2 fold changes (ours vs authors'):", round(auth_cor, 3), "\n")
}


# ── 5. QUALITY CONTROL: PCA ──────────────────────────────────
pca <- prcomp(t(log_expr), scale. = TRUE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Condition = sample_info$condition)
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(Control = "#4DBBD5", PD = "#E64B35")) +
  labs(title = "PCA of log2 expression", x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)")) +
  theme_classic(base_size = 13)
ggsave(out("PCA_PD.png"), pca_plot, width = 7, height = 5.5, dpi = 300)


# ── 6. VOLCANO PLOT ───────────────────────────────────────────
res_df$significance <- "Not Significant"
res_df$significance[res_df$log2FoldChange >  LFC_CUT & res_df$padj < FDR_CUT] <- "Upregulated"
res_df$significance[res_df$log2FoldChange < -LFC_CUT & res_df$padj < FDR_CUT] <- "Downregulated"

top_genes <- res_df %>% filter(padj < FDR_CUT) %>% arrange(padj) %>% head(15)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#E64B35",
                                "Downregulated" = "#4DBBD5",
                                "Not Significant" = "grey70")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", linewidth = 0.4) +
  labs(title = "Differential Gene Expression in Parkinson's Disease",
       subtitle = "PD vs Control (GSE68719)",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", color = "") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")
ggsave(out("volcano_plot_PD.png"), volcano_plot, width = 8, height = 6, dpi = 300)
print(volcano_plot)


# ── 7. GENE ID CONVERSION (Symbol → Entrez) ───────────────────
sig_genes <- res_df %>% filter(padj < FDR_CUT, abs(log2FoldChange) > LFC_CUT)
entrez_ids <- bitr(sig_genes$gene, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
sig_genes <- merge(sig_genes, entrez_ids, by.x = "gene", by.y = "SYMBOL")


# ── 8. KEGG PATHWAY ENRICHMENT (needs internet) ───────────────
kegg_res <- enrichKEGG(gene = sig_genes$ENTREZID, organism = "hsa",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2)
kegg_plot <- dotplot(kegg_res, showCategory = 15) +
  labs(title = "KEGG Pathway Enrichment – Parkinson's Disease",
       subtitle = "Significantly dysregulated pathways (padj < 0.05)") +
  theme_classic(base_size = 12)
ggsave(out("KEGG_dotplot_PD.png"), kegg_plot, width = 9, height = 7, dpi = 300)
print(kegg_plot)


# ── 9. GO ENRICHMENT (Biological Process) ─────────────────────
go_res <- enrichGO(gene = sig_genes$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                   readable = TRUE)
go_plot <- dotplot(go_res, showCategory = 15) +
  labs(title = "GO Biological Process Enrichment – Parkinson's Disease",
       subtitle = "Top enriched biological processes") +
  theme_classic(base_size = 12)
ggsave(out("GO_dotplot_PD.png"), go_plot, width = 10, height = 7, dpi = 300)
print(go_plot)


# ── 10. GSEA (all genes ranked by log2 fold change) ───────────
all_entrez <- bitr(res_df$gene, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
res_ranked <- merge(res_df, all_entrez, by.x = "gene", by.y = "SYMBOL")
res_ranked <- res_ranked[!is.na(res_ranked$log2FoldChange), ]
gene_list <- res_ranked$log2FoldChange
names(gene_list) <- res_ranked$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

gsea_res <- gseKEGG(geneList = gene_list, organism = "hsa", minGSSize = 10,
                    maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsea_plot <- gseaplot2(gsea_res, geneSetID = 1:3,
                       title = "GSEA – Top Enriched KEGG Pathways in Parkinson's Disease")
ggsave(out("GSEA_plot_PD.png"), gsea_plot, width = 10, height = 7, dpi = 300)
print(gsea_plot)


# ── 11. METABOLIC PATHWAY CHECK (keyword filter) ──────────────
# Reported whatever the result; check the matched terms by eye.
metab_pat <- "metabol|glycoly|oxidative phosph|citrate|fatty acid|lipid|amino acid|glucose|pyruvate|biosynthesis|insulin|AMPK|PPAR"
kegg_df <- as.data.frame(kegg_res)
gsea_df <- as.data.frame(gsea_res)
kegg_metab <- kegg_df[grepl(metab_pat, kegg_df$Description, ignore.case = TRUE), ]
gsea_metab <- gsea_df[grepl(metab_pat, gsea_df$Description, ignore.case = TRUE), ]
write.csv(kegg_metab, out("KEGG_metabolic_PD.csv"), row.names = FALSE)
write.csv(gsea_metab, out("GSEA_metabolic_PD.csv"), row.names = FALSE)


# ── 12. HEATMAP OF TOP 30 DEGs ────────────────────────────────
top30 <- res_df %>% filter(padj < FDR_CUT) %>% arrange(padj) %>% head(30) %>% pull(gene)
mat <- log_expr[top30, ]
mat <- mat - rowMeans(mat)                      # centre each gene
annotation_col <- data.frame(Condition = sample_info$condition, row.names = colnames(mat))
ann_colors <- list(Condition = c(Control = "#4DBBD5", PD = "#E64B35"))
heat_args <- list(
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
# Saved with pheatmap's own file option (more reliable than png()/dev.off())
do.call(pheatmap, c(heat_args, list(filename = out("heatmap_top30_PD.png"),
                                    width = 10, height = 8)))
do.call(pheatmap, heat_args)                    # also draw in the Plots pane


# ── 13. RESULTS SUMMARY (written to a text file) ──────────────
n_up   <- sum(res_df$significance == "Upregulated")
n_down <- sum(res_df$significance == "Downregulated")
summary_lines <- c(
  "Results summary – GSE68719, BA9 prefrontal cortex, PD vs Control",
  paste("Samples: Control =", sum(sample_info$condition == "Control"),
        ", PD =", sum(sample_info$condition == "PD")),
  paste("Genes tested after filtering:", nrow(res_df)),
  paste0("Significant (padj < ", FDR_CUT, ", |log2FC| > ", LFC_CUT, "): ",
         n_up + n_down, " (up in PD: ", n_up, "; down in PD: ", n_down, ")"),
  paste("Correlation of log2FC with authors' published results:", round(auth_cor, 3)),
  paste("Top 10 genes by adjusted p-value:", paste(head(res_df$gene, 10), collapse = ", ")),
  paste("Top KEGG terms (over-representation):",
        paste(head(kegg_df$Description, 5), collapse = "; ")),
  paste("Top GO terms:", paste(head(as.data.frame(go_res)$Description, 5), collapse = "; ")),
  paste("Top GSEA KEGG terms:", paste(head(gsea_df$Description, 5), collapse = "; ")),
  paste("Metabolic-keyword KEGG terms found:", nrow(kegg_metab),
        "| GSEA terms:", nrow(gsea_metab)),
  "Limitations: no adjustment for age, post-mortem interval or RNA quality; bulk tissue, so",
  "cell-type composition can contribute to immune signals."
)
writeLines(summary_lines, out("results_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
writeLines(capture.output(sessionInfo()), out("session_info.txt"))
