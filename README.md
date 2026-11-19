# Bioinformatics Portfolio

RNA-seq differential expression and pathway analysis workflow, demonstrated on **public** Parkinson's disease data. My client and collaborator analyses are confidential and are not shown here.

**Author:** Hammad Ali Hassan · ORCID [0000-0002-9174-7673](https://orcid.org/0000-0002-9174-7673)

## Project: Transcriptomic signature of Parkinson's disease in prefrontal cortex

**Dataset:** [GSE68719](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68719) (NCBI GEO), post-mortem BA9 prefrontal cortex, Parkinson's disease (n = 29) vs control (n = 44).

**Script:** [`PD_analysis.R`](PD_analysis.R)

### Methods
- Expression input: DESeq2-normalized counts for protein-coding genes (GEO supplementary file). Because the values are already normalized, not raw counts, I used **limma** on log2(count + 1) with trend and robust empirical Bayes moderation. DESeq2 would not be appropriate for this input.
- Genes expressed in fewer than 10 samples were filtered out.
- Quality control: PCA.
- Significance: adjusted p < 0.05 and |log2 fold change| > 0.5.
- Enrichment: KEGG and GO (biological process) over-representation, plus KEGG GSEA on all genes ranked by log2 fold change (clusterProfiler).
- Validation: log2 fold changes compared with the original authors' published DESeq2 results.

### Results
<!-- Fill the numbers from results/results_summary.txt -->
- **N** genes passed filtering; **N** were significant (**N** up in PD, **N** down).
- Fold changes correlate with the original authors' results (r = **0.xx**). Strongly ranked genes such as HSPB1 and FXYD2 also appear among the authors' top hits.
- Enrichment points to **immune and inflammatory signalling**: cytokine-cytokine receptor interaction, MAPK signalling, T-cell differentiation, complement and coagulation (KEGG); cytokine production and lymphocyte activation and proliferation (GO); immune and infection gene sets in GSEA. This is consistent with neuroinflammation reported in Parkinson's disease.
- No metabolic pathway stood out.

| Volcano plot | Top 30 genes |
|---|---|
| ![Volcano plot](results/volcano_plot_PD.png) | ![Heatmap of top 30 genes](results/heatmap_top30_PD.png) |

| KEGG enrichment | GO enrichment |
|---|---|
| ![KEGG dotplot](results/KEGG_dotplot_PD.png) | ![GO dotplot](results/GO_dotplot_PD.png) |

![GSEA](results/GSEA_plot_PD.png)

### Limitations
- The model does not adjust for age, post-mortem interval, or RNA quality. The original authors did.
- Bulk tissue: cell-type composition (for example, more microglia) can contribute to immune signals.
- This is a workflow demonstration, so the findings are exploratory and not a claim about disease mechanism or treatment.

## How to run
1. Install R (4.4 or newer) and the packages listed at the top of `PD_analysis.R`.
2. Set the working directory to this folder and run `source("PD_analysis.R")`.
3. The script downloads the data into `data/` and writes all outputs to `results/`. KEGG steps need an internet connection.

## Files
| Path | Content |
|---|---|
| `PD_analysis.R` | Full analysis script |
| `results/` | Figures, full limma results table, `results_summary.txt`, `session_info.txt` |
