# ============================================================
# # DE Advisor
# Copyright (c) 2026 Hammad Ali Hassan
# All rights reserved.
#
# This source code is publicly available for inspection only.
# No permission is granted to copy, modify, redistribute,
# or commercially use this software without prior written
# permission from the copyright holder.
#
# DE Advisor 1.2 – check your expression data, then run the right test
# A local Shiny app: your data never leaves your computer.
# Author: Hammad Ali Hassan (ORCID 0000-0002-9174-7673)
#
# What it does
#   1. Reads an expression table (genes x samples) and a sample sheet.
#   2. Detects the data type (raw counts / normalized / log-scale) and
#      recommends DESeq2 or limma accordingly.
#   3. Runs the test, then shows results, volcano plot, PCA, heatmap and
#      (optional) GO / KEGG enrichment for human gene symbols.
#   4. Writes a plain-language methods paragraph you can paste in a report.
#
# Run:  open app.R in RStudio and press Ctrl+Shift+Enter (see README_app.md)
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})
options(shiny.maxRequestSize = 300 * 1024^2)

# ── helpers ──────────────────────────────────────────────────
read_any <- function(path, name) {
  df <- if (grepl("\\.csv(\\.gz)?$", name, ignore.case = TRUE)) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (anyDuplicated(names(df))) stop("The uploaded table contains duplicate column names. Rename the columns before importing.")
  if (!nrow(df) || !ncol(df)) stop("The uploaded table is empty.")
  df
}

detect_type <- function(x) {
  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if (!length(values)) stop("The expression matrix contains no finite numeric values.")
  neg <- any(values < 0)
  is_int <- all(abs(values - round(values)) < 1e-8)
  mx <- max(values); mn <- min(values); zero_fr <- mean(values == 0)
  # Value patterns cannot prove that integer-valued data are raw counts.
  if (neg) {
    type <- "scaled"; rec <- "limma"
    msg <- "Negative values are present. These may be log-transformed or scaled values; confirm the data type before analysis."
  } else if (is_int && mx > 100) {
    type <- "raw_candidate"; rec <- "deseq2"
    msg <- "Non-negative integer values are compatible with raw counts, but integer values alone cannot confirm this. Verify the data source."
  } else if (mn >= 0 && mx <= 30) {
    type <- "log_candidate"; rec <- "limma"
    msg <- "The value range may be compatible with log-scale expression, but the scale cannot be established from values alone. Verify the data source."
  } else {
    type <- "normalized_candidate"; rec <- "limma"
    msg <- "Values may represent normalized expression. Confirm whether they are raw counts, normalized counts, or transformed values."
  }
  list(type = type, rec = rec, msg = msg, neg = neg, is_int = is_int,
       range = range(values), zero_frac = zero_fr)
}

# GEO metadata lookup is advisory: accession metadata never overrides the
# user's confirmation of the actual file that was downloaded.
fetch_geo_metadata <- function(accession) {
  accession <- toupper(trimws(accession))
  if (!grepl("^GSE[0-9]+$", accession)) {
    return(list(ok = FALSE, accession = accession,
                message = "Enter a GEO Series accession in the format GSE12345."))
  }
  # Curated record for the dataset used during DE Advisor development.
  if (identical(accession, "GSE68719")) {
    return(list(ok = TRUE, accession = accession,
      title = "mRNA-Seq expression and MS3 proteomics profiling of human post-mortem BA9 brain tissue for Parkinson Disease and neurologically normal individuals",
      organism = "Homo sapiens", platform = "GPL11154 · Illumina HiSeq 2000",
      samples = 73, group_hint = "Sample names use C_ for neurologically normal controls and P_ for Parkinson disease samples.",
      files = "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz (submitter-provided DESeq2-normalized counts); GSE68719_mlpd_DESeq2_diffexp_pmi_age_rin_default.txt.gz (published differential-expression results)",
      recommendation = "For the confirmed normalized-count file: select Normalized counts and use limma on log2(value + 1) values with trend-aware empirical Bayes moderation. Do not use DESeq2 or voom as if this were a raw-count matrix.",
      caveat = "The published differential-expression file name indicates a model incorporating post-mortem interval (PMI), age and RNA integrity number (RIN). This DE Advisor release currently fits a two-group model only and does not adjust for those covariates. Check the study design and covariates before interpreting a reanalysis.",
      source = "NCBI GEO Series GSE68719"))
  }
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession,
                "&targ=self&form=text&view=brief")
  lines <- tryCatch(readLines(url, warn = FALSE, encoding = "UTF-8"),
                    error = function(e) character())
  if (!length(lines)) return(list(ok = FALSE, accession = accession,
    message = "GEO lookup could not be completed. Check your internet connection or enter the accession again. You can still use the app, but verify the file type from the GEO record."))
  get_field <- function(key) {
    hit <- grep(paste0("^!Series_", key, "\\s*="), lines, value = TRUE)
    if (!length(hit)) return(NA_character_)
    sub(paste0("^!Series_", key, "\\s*=\\s*"), "", hit[1])
  }
  samples <- grep("^!Series_sample_id\\s*=", lines, value = TRUE)
  list(ok = TRUE, accession = accession,
       title = get_field("title"), organism = get_field("organism"),
       platform = paste(unique(sub("^!Series_platform_id\\s*=\\s*", "", grep("^!Series_platform_id\\s*=", lines, value = TRUE))), collapse = ", "),
       samples = length(samples), files = NA_character_,
       recommendation = "GEO Series metadata were retrieved, but the uploaded supplementary file has not been identified automatically. Inspect the exact file description and confirm whether it contains raw counts, normalized values, or transformed expression before selecting a method.",
       caveat = "An accession identifies a study, not necessarily the provenance or scale of the particular file uploaded. GEO supplementary files may contain different data types.",
       source = paste0("NCBI GEO Series ", accession))
}

# ── UI ───────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      :root { --da-ink:#172b4d; --da-muted:#64748b; --da-teal:#087f8c; --da-teal-dark:#06636d; --da-bg:#f4f7fb; --da-border:#e2e8f0; }
      body { background:var(--da-bg); color:var(--da-ink); font-family:Inter,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif; }
      .container-fluid { max-width:1600px; padding:22px 28px 36px; }
      .da-hero { background:linear-gradient(120deg,#102a43 0%,#164e63 58%,#087f8c 100%); color:#fff; border-radius:18px; padding:25px 30px; margin:0 0 20px; box-shadow:0 8px 24px rgba(15,42,67,.12); }
      .da-brand { font-size:28px; font-weight:750; letter-spacing:-.6px; margin:0 0 5px; }
      .da-subtitle { color:#d8edf0; font-size:14px; margin:0; }
      .da-local { display:inline-block; margin-top:14px; padding:6px 10px; border:1px solid rgba(255,255,255,.25); border-radius:20px; color:#e6fffb; font-size:12px; }
      .da-section-title { font-size:15px; font-weight:700; color:var(--da-ink); margin:2px 0 12px; }
      .sidebar { background:#fff; border:1px solid var(--da-border); border-radius:15px; padding:20px 18px; box-shadow:0 3px 12px rgba(15,42,67,.045); }
      .main-panel { background:#fff; border:1px solid var(--da-border); border-radius:15px; padding:20px 22px; min-height:650px; box-shadow:0 3px 12px rgba(15,42,67,.045); }
      .control-label { color:#334155; font-weight:650; font-size:12px; }
      .help-block { color:var(--da-muted); font-size:12px; }
      .form-control { border:1px solid #cbd5e1; border-radius:8px; box-shadow:none; min-height:38px; }
      .form-control:focus { border-color:var(--da-teal); box-shadow:0 0 0 3px rgba(8,127,140,.12); }
      .btn-primary { background:var(--da-teal); border-color:var(--da-teal); border-radius:8px; font-weight:650; padding:10px 16px; }
      .btn-primary:hover,.btn-primary:focus { background:var(--da-teal-dark); border-color:var(--da-teal-dark); }
      .btn-default { border-radius:7px; }
      .nav-tabs { border-bottom:1px solid var(--da-border); margin-bottom:18px; }
      .nav-tabs>li>a { color:#52647a; font-weight:600; border:0!important; border-bottom:3px solid transparent!important; border-radius:0; padding:12px 15px; }
      .nav-tabs>li>a:hover { background:transparent; color:var(--da-teal); }
      .nav-tabs>li.active>a,.nav-tabs>li.active>a:hover,.nav-tabs>li.active>a:focus { color:var(--da-teal-dark); background:transparent; border-bottom:3px solid var(--da-teal)!important; }
      .shiny-output-error { color:#b42318; }
      .sidebar, .col-sm-8 { background:#fff; border:1px solid var(--da-border); border-radius:15px; padding:20px 22px; box-shadow:0 3px 12px rgba(15,42,67,.045); }
      .well { background:#fff; border:1px solid var(--da-border); border-radius:10px; box-shadow:none; }
      hr { border-top:1px solid var(--da-border); }
      @media (max-width: 767px) { .container-fluid { padding:12px; } .da-hero { padding:20px; } .da-brand { font-size:24px; } .main-panel { padding:14px; min-height:0; } }
    "))
  ),
  div(class = "da-hero",
      h1(class = "da-brand", "DE Advisor 1.2"),
      p(class = "da-subtitle", "Differential expression analysis with transparent, evidence-aware guidance"),
      span(class = "da-local", "LOCAL ANALYSIS  ·  Your expression files stay on this computer")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      div(class = "da-section-title", "01  ·  Prepare your data"),
      fileInput("expr", "Expression table (CSV/TSV, genes × samples)"),
      textInput("geo_accession", "Optional GEO Series accession", placeholder = "e.g. GSE68719"),
      actionButton("check_geo", "Check GEO study", class = "btn-default"),
      uiOutput("geo_ui"),
      uiOutput("gene_col_ui"),
      selectInput("data_type", "Expression data type (verify before analysis)",
                  choices = c("Automatic assessment (confirmation required)" = "auto",
                              "Raw unnormalized counts" = "raw",
                              "Normalized counts / TPM / FPKM / CPM" = "normalized",
                              "Log-transformed expression" = "log",
                              "Scaled / z-score expression" = "scaled"),
                  selected = "auto"),
      uiOutput("duplicate_ui"),
      div(class = "da-section-title", "02  ·  Define your experiment"),
      radioButtons("group_mode", "How will you define the groups?",
                   c("From sample names" = "names", "From a sample sheet" = "file")),
      conditionalPanel(
        "input.group_mode == 'names'",
        helpText("Samples whose names contain the comparison text form the comparison group; all others form the reference group."),
        textInput("case_pat", "Comparison group: sample names containing", "P_"),
        textInput("case_lab", "Comparison group label", "PD"),
        textInput("ref_lab", "Reference group label", "Control")
      ),
      conditionalPanel(
        "input.group_mode == 'file'",
        fileInput("meta", "Sample sheet (CSV/TSV)"),
        uiOutput("meta_ui"),
        uiOutput("group_ui"),
        uiOutput("covariate_ui")
      ),
      hr(),
      div(class = "da-section-title", "03  ·  Analysis settings"),
      radioButtons("method", "Statistical method",
                   c("Automatic recommendation" = "auto",
                     "DESeq2 · raw counts only" = "deseq2",
                     "limma · normalized / log values" = "limma")),
      numericInput("lfc", "Absolute log2 fold-change cutoff", 0.5, min = 0, step = 0.1),
      numericInput("fdr", "Adjusted p-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
      numericInput("min_n", "Minimum samples expressing a gene", 3, min = 1),
      numericInput("min_log_expr", "Minimum log2(expression + 1) for normalized-data filtering", 1, min = 0, step = 0.25),
      helpText("Filtering is applied across all samples before testing. For normalized counts, a gene must meet this threshold in at least the selected number of samples."),
      radioButtons("limma_test", "limma significance test",
                   c("Empirical Bayes + post-hoc fold-change cutoff (v1.1-compatible)" = "ebayes",
                     "Treat: test against the minimum fold change" = "treat"),
                   selected = "treat"),
      helpText("Treat tests whether the true effect exceeds the selected log2 fold-change threshold. It is not equivalent to filtering eBayes results after testing."),
      checkboxInput("pca_scale", "Scale genes to unit variance for PCA", TRUE),
      helpText("Scaling gives each gene equal variance in PCA; without scaling, high-variance genes contribute more."),
      sliderInput("n_top", "Genes displayed in heatmap", 10, 100, 30, step = 5),
      hr(),
      div(class = "da-section-title", "04  ·  Optional enrichment"),
      checkboxInput("do_enrich", "Run enrichment (currently human gene symbols)", FALSE),
      radioButtons("enrich_db", NULL, c("GO biological process" = "go",
                                        "KEGG (requires internet)" = "kegg"), inline = TRUE),
      actionButton("run", "Run differential expression analysis", class = "btn-primary", width = "100%")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Data check", br(), h4("Data quality and analysis readiness"),
                 p(class = "help-block", "Review the data assessment and sample-group checks before running the analysis."),
                 verbatimTextOutput("check_txt"),
                 hr(), h4("Duplicate gene identifiers"),
                 p(class = "help-block", "Inspect repeated identifiers before choosing how they should be handled."),
                 tableOutput("duplicate_tbl"), downloadButton("dl_duplicates", "Download duplicate rows (CSV)")),
        tabPanel("Results", br(), downloadButton("dl_res", "Download results (CSV)"),
                 p(class = "help-block", "P-values are displayed in scientific notation. The CSV preserves full available numeric precision; machine-underflow values may still be represented as zero."),
                 br(), tableOutput("res_tbl")),
        tabPanel("Volcano", br(), plotOutput("volcano", height = "520px"),
                 downloadButton("dl_volcano", "Download PNG")),
        tabPanel("PCA", br(), plotOutput("pca", height = "480px"),
                 downloadButton("dl_pca", "Download PNG")),
        tabPanel("Heatmap", br(), plotOutput("heat", height = "620px"),
                 downloadButton("dl_heat", "Download PNG")),
        tabPanel("Enrichment", br(), plotOutput("enrich_plot", height = "560px")),
        tabPanel("Methods", br(), h4("Analysis methods and settings"),
                 verbatimTextOutput("methods_txt"),
                 downloadButton("dl_methods", "Download methods text"))
      )
    )
  )
)
# ── Server ───────────────────────────────────────────────────
server <- function(input, output, session) {

  geo_record <- eventReactive(input$check_geo, {
    req(input$geo_accession)
    fetch_geo_metadata(input$geo_accession)
  }, ignoreInit = TRUE)

  output$geo_ui <- renderUI({
    g <- geo_record()
    if (is.null(g)) return(helpText("GEO lookup is optional. Enter a GSE accession and select Check GEO study to retrieve study context."))
    if (!isTRUE(g$ok)) return(div(class = "well", strong("GEO lookup"), p(g$message)))
    tagList(
      div(class = "well",
        strong(paste("GEO study:", g$accession)),
        p(g$title),
        p(paste("Organism:", ifelse(is.na(g$organism) || !nzchar(g$organism), "Not available in brief record", g$organism))),
        p(paste("Platform:", ifelse(is.na(g$platform) || !nzchar(g$platform), "Not available in brief record", g$platform))),
        p(paste("Samples listed:", ifelse(is.null(g$samples) || is.na(g$samples), "Not available", g$samples))),
        if (!is.null(g$group_hint)) p(strong("Sample naming clue: "), g$group_hint),
        if (!is.null(g$files) && !is.na(g$files)) p(paste("Relevant GEO files:", g$files)),
        p(strong("Workflow guidance: "), g$recommendation),
        p(strong("Important limitation: "), g$caveat),
        p(class = "help-block", paste("Source:", g$source)),
        if (identical(g$accession, "GSE68719")) p(class = "help-block", "For this study, the known supplementary normalized-count file is submitter-provided. Confirm that this is the file you uploaded." )
      )
    )
  })

  expr_raw <- reactive({
    req(input$expr)
    read_any(input$expr$datapath, input$expr$name)
  })

  output$gene_col_ui <- renderUI({
    df <- expr_raw()
    nonnum <- names(df)[!sapply(df, is.numeric)]
    validate(need(length(nonnum) > 0,
                  "The table needs a text column with gene IDs or symbols."))
    selectInput("gene_col", "Gene ID / symbol column", nonnum,
                selected = if ("symbol" %in% nonnum) "symbol" else nonnum[1])
  })

  expr_input <- reactive({
    df <- expr_raw()
    req(input$gene_col)
    validate(need(input$gene_col %in% names(df), "Select a valid gene identifier column."))
    ids <- trimws(as.character(df[[input$gene_col]]))
    numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, input$gene_col)
    validate(need(length(numeric_cols) >= 4,
                  "Need at least 4 numeric expression columns. Check whether sample columns were imported as numeric."))
    mat <- as.matrix(df[, numeric_cols, drop = FALSE])
    storage.mode(mat) <- "double"
    validate(need(nrow(mat) > 0, "The expression table has no gene rows."),
             need(all(!is.na(ids) & nzchar(ids)), "Missing or blank gene identifiers were found. Correct them before analysis."),
             need(all(is.finite(mat)), "Missing, infinite, or invalid expression values were found. Correct them before analysis."),
             need(!is.null(colnames(mat)) && all(nzchar(colnames(mat))), "Every expression sample column must have a name."),
             need(!anyDuplicated(colnames(mat)), "Duplicate sample column names were found. Rename them before analysis."))
    list(ids = ids, mat = mat)
  })

  duplicate_rows <- reactive({
    z <- expr_input(); ids <- z$ids
    idx <- which(duplicated(ids) | duplicated(ids, fromLast = TRUE))
    if (!length(idx)) return(data.frame())
    out <- data.frame(row_in_file = idx, gene_id = ids[idx],
                      mean_expression = rowMeans(z$mat[idx, , drop = FALSE]),
                      z$mat[idx, , drop = FALSE], check.names = FALSE)
    out[order(out$gene_id, out$row_in_file), , drop = FALSE]
  })

  output$duplicate_ui <- renderUI({
    d <- duplicate_rows()
    if (!nrow(d)) return(helpText("No duplicate gene identifiers detected."))
    tagList(
      div(class = "da-section-title", "Data-integrity check · duplicate identifiers"),
      p(sprintf("%d rows belong to repeated gene identifiers (%d unique identifiers).", nrow(d), length(unique(d$gene_id)))),
      selectInput("duplicate_policy", "How should DE Advisor handle these rows?",
        choices = c("Stop and let me resolve them in my file" = "stop",
                    "Keep the row with the highest mean expression (legacy v1.0 behavior)" = "highest",
                    "Sum rows (raw counts only; confirm rows are non-overlapping features)" = "sum",
                    "Average rows (normalized / transformed data; use only if justified)" = "mean"),
        selected = "stop"),
      helpText("No rows are changed silently. The selected rule is recorded in the methods output. Summing is only enabled for raw counts; averaging is only enabled for non-raw data." )
    )
  })

  output$duplicate_tbl <- renderTable({
    d <- duplicate_rows()
    if (!nrow(d)) return(data.frame(Status = "No duplicate gene identifiers detected."))
    head(d, 100)
  }, digits = 4, rownames = FALSE)

  output$dl_duplicates <- downloadHandler(
    filename = function() "duplicate_gene_rows.csv",
    content = function(file) {
      d <- duplicate_rows()
      if (!nrow(d)) d <- data.frame(Status = "No duplicate gene identifiers detected.")
      write.csv(d, file, row.names = FALSE)
    }
  )

  expr_mat <- reactive({
    z <- expr_input(); ids <- z$ids; mat <- z$mat
    dup <- unique(ids[duplicated(ids)])
    if (!length(dup)) { rownames(mat) <- ids; return(mat) }
    policy <- input$duplicate_policy %||% "stop"
    validate(need(policy != "stop", paste0("Duplicate gene identifiers were found (", length(dup), " identifiers; examples: ",
      paste(head(dup, 5), collapse = ", "), "). Review the Duplicate gene identifiers section and select a handling rule, or resolve them in your source file.")))
    if (policy == "sum") {
      validate(need((input$data_type %||% "auto") == "raw",
        "Summing duplicate rows is available only when you explicitly select Raw unnormalized counts."))
      agg <- rowsum(mat, group = ids, reorder = FALSE)
    } else if (policy == "mean") {
      validate(need((input$data_type %||% "auto") != "raw",
        "Averaging duplicate rows is not offered for raw counts. Use a justified count-summing strategy or resolve the source file."))
      spl <- split(seq_along(ids), ids)
      agg <- do.call(rbind, lapply(spl, function(ii) colMeans(mat[ii, , drop = FALSE])))
      rownames(agg) <- names(spl)
    } else {
      # Explicit compatibility option: retain the highest mean-expression row per identifier.
      ord <- order(-rowMeans(mat), seq_along(ids))
      keep <- !duplicated(ids[ord])
      selected <- ord[keep]
      agg <- mat[selected, , drop = FALSE]
      rownames(agg) <- ids[selected]
    }
    agg
  })

  meta_raw <- reactive({
    req(input$meta)
    read_any(input$meta$datapath, input$meta$name)
  })

  output$meta_ui <- renderUI({
    m <- meta_raw()
    tagList(
      selectInput("sample_col", "Sample ID column", names(m), selected = names(m)[1]),
      selectInput("group_col", "Group column", names(m),
                  selected = names(m)[min(2, ncol(m))])
    )
  })

  output$group_ui <- renderUI({
    req(input$group_col)
    m <- meta_raw()
    g <- unique(trimws(as.character(m[[input$group_col]])))
    g <- g[!is.na(g) & nzchar(g)]
    tagList(
      selectInput("ref", "Reference group (e.g. control)", g, selected = g[1]),
      selectInput("case", "Comparison group (e.g. disease)", g, selected = g[min(2, length(g))])
    )
  })

  output$covariate_ui <- renderUI({
    req(input$sample_col, input$group_col)
    m <- meta_raw()
    candidates <- setdiff(names(m), c(input$sample_col, input$group_col))
    tagList(
      selectInput("covariate_cols", "Optional covariates to adjust for", choices = candidates,
                  selected = character(0), multiple = TRUE),
      helpText("Choose only measured, scientifically justified variables (for example age, sex, PMI, RIN or batch). Leave blank for a group-only model. Numeric columns with more than five distinct values are modelled continuously; text or low-cardinality columns are categorical factors. Review coding before analysis.")
    )
  })

  matched <- reactive({
    sample_names <- colnames(expr_mat())
    if (input$group_mode == "names") {
      req(nzchar(input$case_pat), nzchar(input$case_lab), nzchar(input$ref_lab))
      ref <- trimws(input$ref_lab); case <- trimws(input$case_lab)
      validate(need(nzchar(ref) && nzchar(case) && ref != case,
                    "Enter two different, non-empty group labels."))
      is_case <- grepl(input$case_pat, sample_names, fixed = TRUE)
      validate(need(any(is_case) && any(!is_case),
                    "The sample-name pattern must identify samples in both comparison groups."))
      return(list(samples = sample_names,
                  group = factor(ifelse(is_case, case, ref), levels = c(ref, case)),
                  covariates = data.frame(row.names = sample_names),
                  covariate_names = character(0), covariate_model_names = character(0),
                  n_meta = length(sample_names), ref = ref, case = case))
    }
    req(input$sample_col, input$group_col, input$ref, input$case)
    m <- meta_raw()
    s <- trimws(as.character(m[[input$sample_col]]))
    g <- trimws(as.character(m[[input$group_col]]))
    dup <- unique(s[duplicated(s)])
    validate(need(input$sample_col != input$group_col, "Select different columns for sample IDs and group labels."),
             need(all(!is.na(s) & nzchar(s)), "The sample sheet contains blank sample IDs."),
             need(all(!is.na(g) & nzchar(g)), "The sample sheet contains blank group labels."),
             need(!length(dup), paste0("The sample sheet contains duplicated sample IDs: ", paste(head(dup, 5), collapse = ", "))),
             need(input$ref != input$case, "Reference and comparison groups must be different."))
    missing_meta <- setdiff(sample_names, s)
    extra_meta <- setdiff(s, sample_names)
    validate(need(length(missing_meta) == 0,
                  paste0("Expression samples missing from the sample sheet: ", paste(head(missing_meta, 8), collapse = ", "))),
             need(length(extra_meta) == 0,
                  paste0("Sample-sheet IDs absent from the expression table: ", paste(head(extra_meta, 8), collapse = ", "))))
    ord <- match(sample_names, s)
    g <- g[ord]
    validate(need(all(g %in% c(input$ref, input$case)),
                  "Every expression sample must belong to one of the two selected groups. Other groups cannot be silently excluded."))
    cov_names <- setdiff(input$covariate_cols %||% character(0), "group")
    cov_df <- data.frame(row.names = sample_names)
    if (length(cov_names)) {
      validate(need(all(cov_names %in% names(m)), "One or more selected covariates are no longer present in the sample sheet."))
      model_names <- make.names(cov_names, unique = TRUE)
      for (j in seq_along(cov_names)) {
        nm <- cov_names[j]
        model_nm <- model_names[j]
        v <- m[[nm]][ord]
        validate(need(!any(is.na(v)) && all(nzchar(trimws(as.character(v)))),
                      paste0("Covariate '", nm, "' contains missing or blank values. Resolve them or remove this covariate.")))
        if (is.numeric(v) && length(unique(v)) > 5) {
          validate(need(all(is.finite(v)), paste0("Covariate '", nm, "' contains non-finite values.")))
          cov_df[[model_nm]] <- as.numeric(v)
        } else {
          cov_df[[model_nm]] <- droplevels(factor(trimws(as.character(v))))
          validate(need(nlevels(cov_df[[model_nm]]) >= 2,
                        paste0("Covariate '", nm, "' has fewer than two observed levels.")))
        }
      }
    }
    list(samples = sample_names, group = factor(g, levels = c(input$ref, input$case)),
         covariates = cov_df, covariate_names = cov_names,
         covariate_model_names = names(cov_df),
         n_meta = length(s), ref = input$ref, case = input$case)
  })

  dtype <- reactive({
    detected <- detect_type(expr_input()$mat)
    selected <- input$data_type %||% "auto"
    if (selected == "auto") return(detected)
    list(type = selected, rec = if (selected == "raw") "deseq2" else "limma",
         msg = paste("User-selected data type:", selected), neg = detected$neg,
         is_int = detected$is_int, range = detected$range, zero_frac = detected$zero_frac)
  })

  output$check_txt <- renderText({
    rawz <- expr_input(); mat <- rawz$mat; dt <- dtype()
    lines <- c(
      sprintf("Genes: %d   Samples (numeric columns): %d", nrow(mat), ncol(mat)),
      paste("First sample columns:", paste(head(colnames(mat), 6), collapse = ", ")),
      sprintf("Value range: %s to %s", format(dt$range[1], digits = 6), format(dt$range[2], digits = 6)),
      sprintf("All whole numbers: %s   Negative values: %s   Share of zeros: %.1f%%",
              dt$is_int, dt$neg, 100 * dt$zero_frac),
      "",
      paste("Assessment:", dt$msg),
      paste("Selected data type:", dt$type),
      paste("Recommended method:", if (dt$rec == "deseq2") "DESeq2" else "limma"),
      "Important: value patterns cannot prove whether data are raw or normalized. Verify the source and select the correct type above."
    )
    g <- tryCatch(geo_record(), error = function(e) NULL)
    if (!is.null(g) && isTRUE(g$ok)) {
      lines <- c(lines, "", paste("GEO accession:", g$accession),
                 paste("GEO study:", g$title),
                 paste("GEO organism:", g$organism),
                 paste("GEO sample count:", g$samples),
                 if (!is.null(g$group_hint)) paste("GEO sample naming clue:", g$group_hint),
                 paste("GEO workflow guidance:", g$recommendation),
                 paste("GEO caveat:", g$caveat))
    }
    mt <- tryCatch(matched(), error = function(e) structure(list(message = conditionMessage(e)), class = "match_error"))
    if (inherits(mt, "match_error")) {
      lines <- c(lines, "", paste("Sample/group validation:", mt$message))
    } else {
      tb <- table(mt$group)
      lines <- c(lines, "",
                 sprintf("Samples considered: %d   assigned to the two groups: %d",
                         mt$n_meta, length(mt$samples)),
                 paste0("Group sizes: ", paste(names(tb), tb, sep = " = ", collapse = ", ")))
      if (any(tb < 3)) lines <- c(lines, "WARNING: fewer than 3 samples in a group; results may have limited power.")
    }
    dups <- unique(expr_input()$ids[duplicated(expr_input()$ids)])
    policy <- input$duplicate_policy %||% "stop"
    lines <- c(lines, "", paste0("Duplicate gene identifiers: ", length(dups), " unique IDs affected"),
               paste0("Duplicate handling: ", if (!length(dups)) "not needed" else switch(policy,
                 stop = "not resolved; analysis is blocked", highest = "highest mean-expression row retained",
                 sum = "rows summed (raw counts only)", mean = "rows averaged (non-raw data)", policy)),
               paste0("Covariates selected: ", if (length(mt$covariate_names)) paste(mt$covariate_names, collapse = ", ") else "none (group-only model)"))
    paste(lines, collapse = "\n")
  })

  # ── analysis ──
  analysis <- eventReactive(input$run, {
    mt <- matched(); dt <- dtype()
    validate(need(input$data_type != "auto",
                  "Please confirm the expression data type above. Automatic assessment is only a heuristic and cannot establish whether data are raw or normalized."),
             need(!is.null(input$lfc) && is.finite(input$lfc) && input$lfc >= 0,
                  "Enter a valid non-negative log2 fold-change cutoff."),
             need(!is.null(input$fdr) && is.finite(input$fdr) && input$fdr > 0 && input$fdr <= 1,
                  "Enter an adjusted p-value cutoff greater than 0 and no greater than 1."),
             need(!is.null(input$min_n) && input$min_n >= 1 && input$min_n == as.integer(input$min_n),
                  "The minimum number of samples for filtering must be a positive whole number."),
             need(!is.null(input$min_log_expr) && is.finite(input$min_log_expr) && input$min_log_expr >= 0,
                  "Enter a valid non-negative expression-filter threshold."))
    tb <- table(mt$group)
    validate(need(mt$ref != mt$case, "Choose two different groups."),
             need(all(tb >= 2),
                  paste0("Each group needs at least 2 matched samples. Found: ",
                         paste(names(tb), tb, sep = " = ", collapse = ", "),
                         ". The expression table has ", ncol(expr_mat()),
                         " sample columns, starting with: ",
                         paste(head(colnames(expr_mat()), 4), collapse = ", "),
                         ". Check that you uploaded the normalized-counts file (one column per sample).")))
    x <- expr_mat()[, mt$samples, drop = FALSE]
    method <- if (input$method == "auto") dt$rec else input$method
    design_data <- data.frame(group = mt$group, mt$covariates, check.names = FALSE)
    design_formula <- stats::reformulate(c(mt$covariate_model_names, "group"))
    design <- stats::model.matrix(design_formula, data = design_data)
    validate(need(qr(design)$rank == ncol(design),
                  "The selected covariates and group are confounded or redundant (the design matrix is not full rank). Remove or reconsider covariates."),
             need(nrow(design) > ncol(design),
                  "There are not enough residual degrees of freedom for the selected covariates and group. Reduce the model."))

    withProgress(message = "Running analysis...", value = 0.3, {
      if (method == "deseq2") {
        validate(need(dt$type == "raw",
                      "DESeq2 requires data explicitly identified as raw counts. Select 'Raw unnormalized counts' only if this is what your file contains."))
        validate(need(requireNamespace("DESeq2", quietly = TRUE),
                      "Install DESeq2: BiocManager::install('DESeq2')"))
        validate(need(all(is.finite(x)) && all(x >= 0), "DESeq2 requires finite, non-negative raw counts."),
                 need(max(x) <= .Machine$integer.max, "At least one count exceeds R's supported integer range for this DESeq2 input path."),
                 need(all(abs(x - round(x)) < 1e-8),
                      "DESeq2 requires integer raw counts. Values will not be rounded automatically. Verify the source file and select the correct expression data type."))
        cnt <- x
        storage.mode(cnt) <- "integer"
        keep <- rowSums(cnt >= 10) >= input$min_n
        cnt <- cnt[keep, , drop = FALSE]
        validate(need(nrow(cnt) > 0, "No genes passed the DESeq2 count filter. Reduce the filter only if scientifically justified."))
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = cnt, colData = data.frame(group = mt$group, mt$covariates,
            row.names = colnames(cnt), check.names = FALSE),
          design = design_formula)
        dds <- DESeq2::DESeq(dds)
        r <- DESeq2::results(dds, contrast = c("group", mt$case, mt$ref))
        res <- data.frame(gene = rownames(r), baseMean = r$baseMean,
                          log2FoldChange = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
        mat <- tryCatch(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE)),
                        error = function(e) log2(DESeq2::counts(dds, normalized = TRUE) + 1))
      } else {
        validate(need(requireNamespace("limma", quietly = TRUE),
                      "Install limma: BiocManager::install('limma')"),
                 need(dt$type != "raw",
                      "This release does not apply ordinary limma directly to raw counts. Use DESeq2; limma-voom is planned for a later release."))
        validate(need(dt$type != "normalized" || all(x >= 0),
                      "Normalized non-negative expression values are required for the log2(value + 1) transformation."))
        do_log <- dt$type == "normalized"
        lx <- if (do_log) log2(x + 1) else x
        keep <- if (do_log) rowSums(lx > input$min_log_expr) >= input$min_n else apply(lx, 1, var, na.rm = TRUE) > 0
        keep[is.na(keep)] <- FALSE
        lx <- lx[keep, , drop = FALSE]
        validate(need(nrow(lx) > 0, "No genes passed the expression filter. Review the data type and filtering threshold."))
        fit <- limma::lmFit(lx, design)
        term_index <- match("group", attr(stats::terms(design_formula), "term.labels"))
        group_coef <- which(attr(design, "assign") == term_index)
        validate(need(length(group_coef) == 1,
                      "Could not identify a single group coefficient. Check the group labels and covariate design."))
        if (identical(input$limma_test, "treat")) {
          fit <- limma::treat(fit, lfc = input$lfc, trend = TRUE, robust = TRUE)
        } else {
          fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
        }
        tt <- limma::topTable(fit, coef = group_coef, number = Inf, sort.by = "P")
        res <- data.frame(gene = rownames(tt), AveExpr = tt$AveExpr,
                          log2FoldChange = tt$logFC, pvalue = tt$P.Value, padj = tt$adj.P.Val)
        mat <- lx
      }
    })
    res <- res[order(res$padj), ]
    res$significance <- "Not significant"
    if (method == "limma" && identical(input$limma_test, "treat")) {
      res$significance[!is.na(res$padj) & res$padj < input$fdr & res$log2FoldChange > 0] <- "Up"
      res$significance[!is.na(res$padj) & res$padj < input$fdr & res$log2FoldChange < 0] <- "Down"
    } else {
      res$significance[!is.na(res$padj) & res$padj < input$fdr & res$log2FoldChange >  input$lfc] <- "Up"
      res$significance[!is.na(res$padj) & res$padj < input$fdr & res$log2FoldChange < -input$lfc] <- "Down"
    }
    list(res = res, mat = mat, group = mt$group, method = method, type = dt$type,
         ref = mt$ref, case = mt$case, lfc = input$lfc, fdr = input$fdr,
         covariate_names = mt$covariate_names,
         limma_test = if (method == "limma") input$limma_test else NA_character_,
         n_up = sum(res$significance == "Up"), n_down = sum(res$significance == "Down"))
  })

  output$res_tbl <- renderTable({
    z <- head(analysis()$res, 100)
    for (nm in intersect(c("pvalue", "padj"), names(z))) {
      z[[nm]] <- ifelse(is.na(z[[nm]]), NA_character_,
                        ifelse(z[[nm]] == 0, "< 2.2e-308 (underflow)",
                               formatC(z[[nm]], format = "e", digits = 3)))
    }
    z
  }, digits = 7, rownames = FALSE)
  output$dl_res <- downloadHandler(
    filename = function() "de_results.csv",
    content  = function(file) {
      old <- options(digits = 17)
      on.exit(options(old), add = TRUE)
      write.csv(analysis()$res, file, row.names = FALSE, na = "")
    })

  # ── plots ──
  volcano_obj <- reactive({
    a <- analysis(); r <- a$res; r <- r[!is.na(r$padj) & is.finite(r$log2FoldChange), ]
    if (nrow(r)) {
      positive <- r$padj[r$padj > 0]
      floor_p <- if (length(positive)) min(positive) / 10 else .Machine$double.xmin
      r$plot_padj <- pmax(r$padj, floor_p)
    } else r$plot_padj <- numeric(0)
    lab <- head(r[r$significance != "Not significant", ], 12)
    p <- ggplot(r, aes(log2FoldChange, -log10(plot_padj), color = significance)) +
      geom_point(alpha = 0.6, size = 1.4) +
      scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", `Not significant` = "grey70")) +
      geom_vline(xintercept = c(-a$lfc, a$lfc), linetype = "dashed", linewidth = 0.4) +
      geom_hline(yintercept = -log10(a$fdr), linetype = "dashed", linewidth = 0.4) +
      labs(title = paste(a$case, "vs", a$ref), subtitle = paste("Method:", toupper(a$method)),
           x = "log2 fold change", y = "-log10 adjusted p-value", color = "") +
      theme_classic(base_size = 13) + theme(legend.position = "top")
    if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab) > 0) {
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3,
                                        max.overlaps = 20, show.legend = FALSE)
    }
    p
  })
  output$volcano <- renderPlot(volcano_obj())
  output$dl_volcano <- downloadHandler(filename = function() "volcano.png",
    content = function(file) ggsave(file, volcano_obj(), width = 8, height = 6, dpi = 300))

  pca_obj <- reactive({
    a <- analysis(); m <- a$mat
    m <- m[apply(m, 1, function(z) is.finite(var(z)) && var(z) > 0), , drop = FALSE]
    validate(need(nrow(m) >= 2 && ncol(m) >= 3,
                  "PCA requires at least two variable genes and three samples."))
    p <- prcomp(t(m), scale. = isTRUE(input$pca_scale))
    pct <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
    df <- data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], Group = a$group)
    ggplot(df, aes(PC1, PC2, color = Group)) + geom_point(size = 3, alpha = 0.85) +
      labs(title = "PCA of analysed samples", x = paste0("PC1 (", pct[1], "%)"),
           y = paste0("PC2 (", pct[2], "%)")) + theme_classic(base_size = 13)
  })
  output$pca <- renderPlot(pca_obj())
  output$dl_pca <- downloadHandler(filename = function() "pca.png",
    content = function(file) ggsave(file, pca_obj(), width = 7, height = 5.5, dpi = 300))

  draw_heat <- function(file = NA, width = 7, height = 7) {
    a <- analysis()
    top <- head(a$res$gene[!is.na(a$res$padj)], input$n_top)
    m <- a$mat[top, , drop = FALSE]
    m <- m - rowMeans(m)
    ann <- data.frame(Group = a$group, row.names = colnames(m))
    pheatmap(m, annotation_col = ann, show_colnames = FALSE, fontsize_row = 8,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             main = paste("Top", length(top), "genes by adjusted p-value"),
             filename = file, width = width, height = height)
  }
  output$heat <- renderPlot(draw_heat())
  output$dl_heat <- downloadHandler(filename = function() "heatmap.png",
    content = function(file) draw_heat(file, 10, 8))

  # ── enrichment (optional) ──
  enrich_res <- reactive({
    req(input$do_enrich)
    a <- analysis()
    validate(need(all(vapply(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"),
                             requireNamespace, logical(1), quietly = TRUE)),
                  "Install clusterProfiler, org.Hs.eg.db and enrichplot via BiocManager."))
    sig <- a$res$gene[a$res$significance != "Not significant"]
    validate(need(length(sig) >= 5, "Fewer than 5 significant genes: not enough for enrichment."))
    ids <- clusterProfiler::bitr(sig, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db::org.Hs.eg.db)
    withProgress(message = "Running enrichment...", value = 0.5, {
      if (input$enrich_db == "go") {
        clusterProfiler::enrichGO(ids$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                  ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
      } else {
        clusterProfiler::enrichKEGG(ids$ENTREZID, organism = "hsa",
                                    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      }
    })
  })
  output$enrich_plot <- renderPlot({
    er <- enrich_res()
    validate(need(!is.null(er) && nrow(as.data.frame(er)) > 0, "No significant terms."))
    enrichplot::dotplot(er, showCategory = 15)
  })

  # ── methods text ──
  methods_text <- reactive({
    a <- analysis()
    ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "n/a")
    tool <- if (a$method == "deseq2") {
      sprintf("Differential expression was tested with DESeq2 (v%s) on raw counts (genes with counts >= 10 in at least %d samples were kept).",
              ver("DESeq2"), input$min_n)
    } else {
      filter_desc <- if (a$type == "normalized") {
        sprintf("genes with log2(value + 1) > %s in at least %d samples", input$min_log_expr, input$min_n)
      } else {
        "genes with non-zero variance across samples"
      }
      if (identical(a$limma_test, "treat")) {
        sprintf("Differential expression was tested with limma (v%s) using treat with trend-aware, robust empirical Bayes moderation on %s values; the null hypothesis was that the absolute log2 fold change does not exceed %s. %s were retained.",
                ver("limma"), if (a$type == "normalized") "log2(value + 1)" else "the supplied", a$lfc, filter_desc)
      } else {
        sprintf("Differential expression was tested with limma (v%s) using eBayes moderation (trend and robust) on %s values; significance used an adjusted p-value cutoff plus a post-hoc absolute log2 fold-change cutoff of %s. %s were retained.",
                ver("limma"), if (a$type == "normalized") "log2(value + 1)" else "the supplied", a$lfc, filter_desc)
      }
    }
    g <- tryCatch(geo_record(), error = function(e) NULL)
    provenance <- if (!is.null(g) && isTRUE(g$ok)) paste0(" GEO accession recorded: ", g$accession,
      ". GEO metadata were used as study context only; the uploaded file type was confirmed by the user and was not inferred from the accession alone.") else ""
    significance_desc <- if (a$method == "limma" && identical(a$limma_test, "treat")) {
      paste0("Genes with treat-adjusted p < ", a$fdr,
             " were called significant against the minimum absolute log2 fold-change threshold of ", a$lfc)
    } else {
      paste0("Genes with adjusted p < ", a$fdr, " and |log2 fold change| > ", a$lfc,
             " were called significant")
    }
    paste0("The expression table was classified as '", a$type, "' data based on the selected user setting or heuristic assessment. ", tool,
           provenance,
           " The comparison was ", a$case, " versus ", a$ref, " (", sum(a$group == a$case), " vs ",
           sum(a$group == a$ref), " samples). ", significance_desc, " (", a$n_up, " up, ",
           a$n_down, " down). Covariates included in the model: ",
           if (length(a$covariate_names)) paste(a$covariate_names, collapse = ", ") else "none (group-only model)",
           ". PCA gene-wise unit-variance scaling was ", if (isTRUE(input$pca_scale)) "enabled" else "disabled",
           ". Duplicate gene handling: ", {
             ndup <- length(unique(expr_input()$ids[duplicated(expr_input()$ids)]))
             if (!ndup) "none required" else switch(input$duplicate_policy,
               highest = "highest mean-expression row retained per identifier (legacy v1.0 rule)",
               sum = "counts summed by identifier", mean = "expression values averaged by identifier",
               stop = "unresolved", "not specified")
           }, ". Analysis performed with DE Advisor in R ", getRversion(), ".")
  })
  output$methods_txt <- renderText({ strwrap(methods_text(), width = 90) |> paste(collapse = "\n") })
  output$dl_methods <- downloadHandler(filename = function() "methods.txt",
    content = function(file) writeLines(methods_text(), file))
}

shinyApp(ui, server)
