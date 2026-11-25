# Bioinformatics Portfolio

RNA-seq differential expression and pathway analysis workflow, demonstrated on **public** Parkinson's disease data. My client and collaborator analyses are confidential and are not shown here.

**Author:** Hammad Ali Hassan · ORCID [0000-0002-9174-7673](https://orcid.org/0000-0002-9174-7673)

Project                                 	What it is
1. Parkinson's disease transcriptomics	  RNA-seq differential expression and pathway analysis in R
2. DE Advisor	                            A local Shiny app that checks the data type, then runs the right test

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

------------------------------------------------------------------------------------------------------------------------------------------

2. DE Advisor: a local app that checks your data before testing

Folder: DE_Advisor/ · Download: see the Releases page

A common mistake in RNA-seq analysis is running DESeq2 on values that are not raw counts. DE Advisor is a local R Shiny app that inspects the uploaded expression table first (raw counts, normalized values, or log-scale), recommends DESeq2 or limma accordingly, and blocks DESeq2 on data that are not raw integer counts.

What it does

Detects the data type and explains the recommendation.
Defines groups from sample names or from a sample sheet, with optional covariates (version 1.2).
Runs DESeq2 or limma, and shows a results table, volcano plot, PCA, and heatmap.
Optional GO and KEGG enrichment for human gene symbols.
Writes a methods paragraph that records the test, filters, and thresholds.

Privacy: it runs on your own computer, so your data are not uploaded anywhere.

Validation: run on GSE68719 (above), it flagged the data as normalized, chose limma, and reproduced the top genes and group sizes (29 PD vs 44 control) from Project 1. The treat option, covariate models, and DESeq2 mode have not yet been validated on real data.

Setup and use (Windows)
Install R (version 4.4 or newer) from https://cran.r-project.org/bin/windows/base/ in the default location.
Extract the downloaded zip to a simple folder path (no apostrophes or special characters). Do not run the files from inside the zip.
First time only: double-click 1_Install_Once.bat and wait until it finishes (10-30 minutes, needs internet).
Every time: double-click 2_Start_DE_Advisor.bat. Keep the black window open while you use the app; the app opens in your browser. Close the black window to stop it.
Optional, for the GO/KEGG enrichment tab: run install_optional_packages.R once (open it in RStudio and click Source).

If R is installed in another location, or you are on Mac or Linux, open app.R in RStudio and click Run App. Tested on Windows with R 4.6.2.

Limitations
Two-group comparison only. Results are exploratory; check the model design before drawing conclusions.
Enrichment is limited to human gene symbols.
Files
Path	Content
PD_analysis.R	Parkinson's analysis script
results/	Figures, full limma results table, results_summary.txt, session_info.txt
DE_Advisor/	The local app and its launchers
