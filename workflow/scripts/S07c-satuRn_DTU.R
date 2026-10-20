## _________________________________________________
##
## Differential transcript usage with satuRn
##
## Aim: Test differential usage of the ENDome bins between sample groups
##      (e.g. PDD vs Control) within each stratum (e.g. brain region) and cell
##      type, from the pseudobulked scUTRquant counts of one ENDome (07b).
##
## Project: ENDome generation - UTR Usage
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-06
##
## Latest Version: v3.0 (2026-10-16)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - Input formats, outputs and the reasoning behind each step are described
##      in docs/UTR_quantification_and_usage.md.
##    - One fit per stratum x cell type, with every group of the contrasts.
##      All contrasts are tested from that fit.
##    - Per fit: samples with >= min_cells nuclei -> groups with >=
##      min_samples_per_group samples -> DRIMSeq::dmFilter -> satuRn::fitDTU ->
##      satuRn::testDTU -> stageR (gene-level screening + bin confirmation).
##    - p_value: "empirical" (satuRn's recommendation; falls back to the raw
##      p-values when satuRn cannot compute them, see the p_type column), "raw",
##      or "none" (no stageR).
##    - Numeric covariates are centred and scaled; any other covariate is a
##      factor. Covariates that are constant within a fit are dropped.
##
## Changelog:
##    - v3.0 (2026-10-16): Simplified. stratify_by is a single column name; the
##      sample table is read as TSV; the constant "fit" column is removed.
##    - v2.0 (2026-10-15): Dataset-independent sample table.
##    - v1.1 (2026-10-07): Results carry the bin annotations of 07b.
##    - v1.0 (2026-10-06): Initial version.
##
## Contact: guillermorocamora@gmail.com
## _________________________________________________

############################################################################## #
# ---- 0. Setup ----

#----------------------------------------------------------------------------- #
## 0.0 Define Snakemake interactive parameters ----
if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots = c(
      input = 'list',
      output = 'list',
      params = 'list',
      wildcards = 'list',
      log = 'list',
      threads = 'numeric',
      scriptdir = 'character')
  )
  ## Inputs linked by the `test_data` rule; outputs use the "interactive"
  ## prefix so they never overwrite them
  test_dir <- "data/test_data/Wood.control.st_ref"
  test_endome <- "ref_pc.w500.3p"
  test_pseudobulk <- paste0(test_endome, ".default")

  snakemake <- Snakemake(
    input = list(
      pseudobulk = file.path(test_dir, paste0("test.", test_pseudobulk, ".pseudobulk.rds")),
      n_cells = file.path(test_dir, paste0("test.", test_pseudobulk, ".pseudobulk_n_cells.tsv")),
      sample_table = "data/Hardy_snRNAseq/prepared/Hardy.sample_table.tsv"
    ),
    output = list(
      results = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.tsv")),
      diagnostics = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.diagnostics.tsv")),
      results_rds = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.rds")),
      diagnostics_rds = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.diagnostics.rds")),
      diagplots = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.diagplots"))
    ),
    params = list(
      exclude_samples = c("gosot_ACG", "gosot_IPL"),
      stratify_by = "region",
      celltypes = NULL,
      contrasts = c("PDD-Control", "PD-Control", "PDD-PD"),
      covariates = c("sex", "rin_bxp", "deletion_length", "insertion_length", "uniquely_mapped_percent"),
      min_cells = 5,
      min_genes = 10,
      min_samples_per_group = 3,
      dmfilter = list(min_samps_gene_prop = 0.75, min_gene_expr = 10, min_feature_expr = 10),
      p_value = "raw",
      alpha = 0.05
    ),
    wildcards = list(
      endome_name = paste0(basename(test_dir), ".", test_endome)
    ),
    log = list(),
    threads = 8,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages

shhh({
  library(conflicted)
  library(tidyverse)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(satuRn)
  library(DRIMSeq)
  library(DEXSeq)
  library(stageR)
  library(limma)
  library(BiocParallel)
})

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflict_prefer_all("tidyr", quiet = TRUE)
conflicted::conflict_prefer("counts", "DRIMSeq", quiet = TRUE)
conflicted::conflict_prefer("samples", "DRIMSeq", quiet = TRUE)

## One BLAS thread per worker: the parallelism is across fits
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
smk_inputs <- snakemake@input
smk_outputs <- snakemake@output
params <- snakemake@params
wc <- snakemake@wildcards

exclude_samples <- as.character(unlist(params$exclude_samples))
stratify_by <- params$stratify_by          # a sample-table column, or NULL for no strata
celltypes <- params$celltypes              # NULL = every cell type of the pseudobulk
contrasts <- params$contrasts
covariates <- as.character(unlist(params$covariates))
min_cells <- params$min_cells
min_genes <- params$min_genes
min_samples_per_group <- params$min_samples_per_group
dmfilter <- params$dmfilter
p_value <- match.arg(params$p_value, c("raw", "empirical", "none"))
alpha <- params$alpha

bpparam <- if (interactive()) {
  BiocParallel::SerialParam(stop.on.error = FALSE)
} else {
  BiocParallel::MulticoreParam(workers = snakemake@threads, stop.on.error = FALSE)
}

dir.create(dirname(smk_outputs$results), showWarnings = FALSE, recursive = TRUE)
dir.create(smk_outputs$diagplots, showWarnings = FALSE, recursive = TRUE)

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Parse contrasts written as "<group1>-<group2>"
#'
#' @param contrasts Character vector, e.g. c("PDD-Control", "PD-Control")
#' @return A tibble with contrast (e.g. "PDDvsControl"), group1 and group2
parse_contrasts <- function(contrasts) {
  parts <- strsplit(contrasts, "-", fixed = TRUE)
  if (any(lengths(parts) != 2)) stop("Contrasts must be '<group1>-<group2>': ", paste(contrasts, collapse = ", "))
  tibble::tibble(
    contrast = vapply(parts, paste, character(1), collapse = "vs"),
    group1 = vapply(parts, `[`, character(1), 1),
    group2 = vapply(parts, `[`, character(1), 2)
  )
}

#' Make a string safe for file names (cell types such as "CUX2_L2/3")
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

#' Empirical null estimates (delta, sigma) of one contrast
#'
#' Only used for the diagnostics table: satuRn computes them internally but
#' does not return them. Reproduces the steps of satuRn:::p.adjust_empirical.
#'
#' @param pval Raw p-values from testDTU
#' @param tstat t-statistics from testDTU
#' @return A named numeric vector c(delta, sigma), NA if the estimation fails
estimate_empirical_null <- function(pval, tstat) {
  tryCatch({
    zvalues <- qnorm(pval / 2) * sign(tstat)
    zvalues_mid <- zvalues[abs(zvalues) < 10]
    zvalues_mid <- zvalues_mid[!is.na(zvalues_mid)]
    N <- length(zvalues_mid)
    b <- 4.3 * exp(-0.26 * log(N, 10))
    med <- median(zvalues_mid)
    sc <- diff(quantile(zvalues_mid)[c(2, 4)]) / (2 * qnorm(0.75))
    mlests <- satuRn:::locfdr_locmle(zvalues_mid, xlim = c(med, b * sc))
    lo <- min(zvalues_mid)
    up <- max(zvalues_mid)
    breaks <- seq(lo, up, length = 120)
    zzz <- pmax(pmin(zvalues_mid, up), lo)
    x <- (breaks[-1] + breaks[-length(breaks)]) / 2
    X <- cbind(1, poly(x, df = 7))
    y <- hist(zzz, breaks = breaks, plot = FALSE)$counts
    f <- glm(y ~ poly(x, df = 7), poisson)$fit
    ml.out <- satuRn:::locfdr_locmle(zvalues_mid, xlim = c(mlests[1], b * mlests[2]),
                                     d = mlests[1], s = mlests[2],
                                     Cov.in = list(x = x, X = X, f = f, sw = 0))
    c(delta = unname(ml.out$mle[1]), sigma = unname(ml.out$mle[2]))
  }, error = function(e) c(delta = NA_real_, sigma = NA_real_))
}

#' Gene-level q-values for the stageR screening stage
#'
#' satuRn vignette recipe: per-gene minimum p-value, then DEXSeq's
#' perGeneQValue. Bins with NA p-values are ignored; genes without any valid
#' p-value get 1.
#'
#' @param pvals Bin-level p-values
#' @param gene_id Gene of each bin
#' @return A named vector of gene q-values
compute_gene_qvalues <- function(pvals, gene_id) {
  gene_id <- factor(gene_id)
  gene_split <- split(seq_along(gene_id), gene_id)
  p_gene <- vapply(gene_split, function(i) suppressWarnings(min(pvals[i], na.rm = TRUE)), numeric(1))
  p_gene[!is.finite(p_gene)] <- 1
  theta <- unique(sort(p_gene))
  q <- DEXSeq:::perGeneQValueExact(p_gene, theta, gene_split)
  q_screen <- pmin(1, q[match(p_gene, theta)])
  names(q_screen) <- names(gene_split)
  q_screen
}

#' Diagnostics rows of one fit: one row per contrast
make_diagnostics <- function(fit, contrast_names, status, ...) {
  tibble::tibble(endome = wc$endome_name, stratum = fit$stratum, celltype = fit$celltype,
                 contrast = contrast_names, status = status, ...)
}

############################################################################## #
# ---- 1. Load inputs ----

#----------------------------------------------------------------------------- #
## 1.1 Pseudobulk and nuclei per pseudobulk (07b) ----
message("Loading the pseudobulk: ", smk_inputs$pseudobulk)
pseudobulk <- readRDS(smk_inputs$pseudobulk)
n_cells <- readr::read_tsv(smk_inputs$n_cells, col_types = readr::cols(.default = readr::col_character(), n_cells = readr::col_double())) %>%
  dplyr::select(sample_id, celltype, n_cells)

if (is.null(celltypes)) celltypes <- assayNames(pseudobulk)
missing_celltypes <- setdiff(celltypes, assayNames(pseudobulk))
if (length(missing_celltypes) > 0) {
  warning("Cell types not in the pseudobulk (skipped): ", paste(missing_celltypes, collapse = ", "))
  celltypes <- intersect(celltypes, assayNames(pseudobulk))
}

#----------------------------------------------------------------------------- #
## 1.2 Sample table ----
contrast_tbl <- parse_contrasts(contrasts)
all_groups <- unique(c(contrast_tbl$group1, contrast_tbl$group2))

samples <- readr::read_tsv(smk_inputs$sample_table, col_types = readr::cols(
  .default = readr::col_guess(), sample_id = readr::col_character(), group = readr::col_character()))

missing_cols <- setdiff(c("sample_id", "group", stratify_by, covariates), colnames(samples))
if (length(missing_cols) > 0) stop("Columns missing from the sample table: ", paste(missing_cols, collapse = ", "))
missing_groups <- setdiff(all_groups, samples$group)
if (length(missing_groups) > 0) stop("Contrast groups not in the sample table: ", paste(missing_groups, collapse = ", "))

## stratum: the value of the stratification column ("all" without stratify_by)
samples <- samples %>%
  dplyr::mutate(stratum = if (is.null(stratify_by)) "all" else as.character(.data[[stratify_by]])) %>%
  dplyr::filter(group %in% all_groups, !sample_id %in% exclude_samples, sample_id %in% colnames(pseudobulk))
if (anyNA(samples[, covariates])) stop("Missing values in the covariates of the sample table")

strata_n <- table(samples$stratum)
message("Samples: ", nrow(samples), " (", paste(names(strata_n), strata_n, sep = "=", collapse = ", "), ")")

#----------------------------------------------------------------------------- #
## 1.3 Fits to run: every stratum x cell type ----
fit_grid <- tidyr::expand_grid(stratum = sort(unique(samples$stratum)), celltype = celltypes)
message("Fits to run: ", nrow(fit_grid))

############################################################################## #
# ---- 2. satuRn DTU ----

#' Run one satuRn fit (one stratum x cell type) and test every contrast
#'
#' @param i Row of fit_grid
#' @return A list with `results` (tibble or NULL) and `diagnostics` (tibble)
run_fit <- function(i) {
  # i <- 1
  start_time <- Sys.time()
  fit <- fit_grid[i, ]

  #----------------------------------------------------------------------------- #
  ## 2.1 Samples of this fit ----
  ## Samples with >= min_cells nuclei of this cell type, then groups with >=
  ## min_samples_per_group such samples. Contrasts with a dropped group are
  ## reported as skipped.
  fit_samples <- samples %>%
    dplyr::filter(stratum == fit$stratum) %>%
    dplyr::inner_join(n_cells %>% dplyr::filter(celltype == fit$celltype) %>% dplyr::select(-celltype), by = "sample_id") %>%
    dplyr::filter(n_cells >= min_cells)

  group_sizes <- table(factor(fit_samples$group, levels = all_groups))
  group_sizes_txt <- paste(names(group_sizes), group_sizes, sep = "=", collapse = ";")
  kept_groups <- names(group_sizes)[group_sizes >= min_samples_per_group]
  fit_samples <- fit_samples %>% dplyr::filter(group %in% kept_groups)
  cts <- assay(pseudobulk, fit$celltype)[, fit_samples$sample_id, drop = FALSE]

  testable <- contrast_tbl$group1 %in% kept_groups & contrast_tbl$group2 %in% kept_groups
  skipped_diag <- make_diagnostics(fit, contrast_tbl$contrast[!testable],
    paste0("skipped: < ", min_samples_per_group, " samples with >= ", min_cells, " nuclei in a group"),
    group_sizes = group_sizes_txt)
  fit_contrasts <- contrast_tbl[testable, ]
  if (nrow(fit_contrasts) == 0) return(list(results = NULL, diagnostics = skipped_diag))

  skip_fit <- function(status) {
    list(results = NULL, diagnostics = dplyr::bind_rows(skipped_diag,
      make_diagnostics(fit, fit_contrasts$contrast, status, group_sizes = group_sizes_txt)))
  }

  #----------------------------------------------------------------------------- #
  ## 2.2 Filter bins with DRIMSeq::dmFilter ----
  ## DRIMSeq renames non-syntactic sample columns (e.g. "001" -> "X001"), so it
  ## gets internal IDs (s1, s2, ...), mapped back to sample_id in 2.4
  n_bins_input <- sum(rowSums(cts) > 0)
  dm_ids <- paste0("s", seq_len(nrow(fit_samples)))
  d <- DRIMSeq::dmDSdata(
    counts = data.frame(gene_id = rowData(pseudobulk)$gene_id, feature_id = rownames(cts),
                        stats::setNames(as.data.frame(as.matrix(cts)), dm_ids), check.names = FALSE, stringsAsFactors = FALSE),
    samples = data.frame(sample_id = dm_ids, group = fit_samples$group, stringsAsFactors = FALSE)
  )
  ## dmFilter fails when no gene passes
  d <- tryCatch(
    DRIMSeq::dmFilter(d,
      min_samps_gene_expr = round(nrow(fit_samples) * dmfilter$min_samps_gene_prop),
      min_samps_feature_expr = min(group_sizes[kept_groups]),
      min_gene_expr = dmfilter$min_gene_expr,
      min_feature_expr = dmfilter$min_feature_expr),
    error = function(e) NULL)
  if (is.null(d)) return(skip_fit("skipped: no genes left after filtering"))

  ## satuRn needs >= 2 bins per gene
  filt <- DRIMSeq::counts(d) %>%
    dplyr::mutate(gene_id = as.character(gene_id), feature_id = as.character(feature_id)) %>%
    dplyr::group_by(gene_id) %>% dplyr::filter(dplyr::n() >= 2) %>% dplyr::ungroup()
  n_genes_filt <- dplyr::n_distinct(filt$gene_id)
  if (n_genes_filt < min_genes) {
    return(skip_fit(paste0("skipped: ", n_genes_filt, " genes after filtering (< min_genes = ", min_genes, ")")))
  }

  #----------------------------------------------------------------------------- #
  ## 2.3 Design matrix: ~ 0 + group + covariates ----
  coldata <- fit_samples %>% dplyr::mutate(group = factor(group, levels = kept_groups))
  used_covariates <- character(0)
  for (cov in covariates) {
    values <- coldata[[cov]]
    if (dplyr::n_distinct(values) < 2) next
    coldata[[cov]] <- if (is.numeric(values)) as.numeric(scale(values)) else factor(values)
    used_covariates <- c(used_covariates, cov)
  }

  design_formula <- as.formula(paste("~ 0 + group", paste(c("", used_covariates), collapse = " + ")))
  coldata <- as.data.frame(coldata)
  rownames(coldata) <- coldata$sample_id
  design <- model.matrix(design_formula, data = coldata)
  if (qr(design)$rank < ncol(design)) return(skip_fit("skipped: rank-deficient design"))

  #----------------------------------------------------------------------------- #
  ## 2.4 Fit satuRn ----
  count_mat <- as.matrix(filt[, dm_ids])
  dimnames(count_mat) <- list(filt$feature_id, fit_samples$sample_id)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = count_mat),
    colData = S4Vectors::DataFrame(coldata),
    rowData = S4Vectors::DataFrame(isoform_id = filt$feature_id, gene_id = filt$gene_id, row.names = filt$feature_id)
  )
  S4Vectors::metadata(se)$formula <- design_formula
  se <- satuRn::fitDTU(object = se, formula = design_formula, parallel = FALSE, verbose = FALSE)

  #----------------------------------------------------------------------------- #
  ## 2.5 Test the contrasts ----
  L <- limma::makeContrasts(
    contrasts = paste0("group", fit_contrasts$group1, "-group", fit_contrasts$group2),
    levels = colnames(design)
  )
  colnames(L) <- fit_contrasts$contrast

  ## satuRn 1.18 always draws its diagnostic plots, so they go to a PDF per fit
  fit_label <- paste(fit$stratum, fit$celltype, sep = " | ")
  grDevices::pdf(file.path(smk_outputs$diagplots, paste0(safe_name(paste(fit$stratum, fit$celltype, sep = ".")), ".pdf")),
                 width = 8, height = 6)
  graphics::plot.new()
  graphics::text(0.5, 0.5, paste0(wc$endome_name, "\n", fit_label, "\n", group_sizes_txt), cex = 1.1)
  se <- tryCatch(
    satuRn::testDTU(object = se, contrasts = L, diagplot1 = TRUE, diagplot2 = TRUE, sort = FALSE),
    finally = grDevices::dev.off())

  #----------------------------------------------------------------------------- #
  ## 2.6 Per contrast: p-values, gene q-values and stageR ----
  per_contrast <- purrr::map(fit_contrasts$contrast, function(contrast_name) {
    # contrast_name <- fit_contrasts$contrast[[1]]
    res <- as.data.frame(rowData(se)[[paste0("fitDTUResult_", contrast_name)]])
    res$isoform_id <- rowData(se)$isoform_id
    res$gene_id <- rowData(se)$gene_id

    ## satuRn 1.18 only computes empirical p-values with >= 500 tested bins
    has_empirical <- any(!is.na(res$empirical_pval))
    p_type <- if (p_value == "empirical" && !has_empirical) "raw" else p_value
    res$p_used <- switch(p_type, empirical = res$empirical_pval, raw = res$pval, none = NA_real_)
    null_est <- if (has_empirical) estimate_empirical_null(res$pval, res$t) else c(delta = NA_real_, sigma = NA_real_)

    if (p_type != "none" && any(!is.na(res$p_used))) {
      q_screen <- compute_gene_qvalues(res$p_used, res$gene_id)
      stage_obj <- stageR::stageRTx(
        pScreen = q_screen,
        pConfirmation = matrix(res$p_used, ncol = 1, dimnames = list(res$isoform_id, "transcript")),
        pScreenAdjusted = TRUE,
        tx2gene = data.frame(transcript = res$isoform_id, gene = res$gene_id)
      )
      stage_obj <- stageR::stageWiseAdjustment(stage_obj, method = "dtu", alpha = alpha, allowNA = TRUE)
      padj <- suppressMessages(stageR::getAdjustedPValues(stage_obj, order = FALSE, onlySignificantGenes = FALSE)) %>%
        dplyr::select(isoform_id = txID, stageR_gene = gene, stageR_transcript = transcript)
      res <- res %>%
        dplyr::mutate(gene_qvalue = base::unname(q_screen[gene_id])) %>%
        dplyr::left_join(padj, by = "isoform_id")
    } else {
      res <- res %>% dplyr::mutate(gene_qvalue = NA_real_, stageR_gene = NA_real_, stageR_transcript = NA_real_)
    }

    results <- res %>%
      dplyr::transmute(
        endome = wc$endome_name, stratum = fit$stratum, celltype = fit$celltype,
        contrast = contrast_name, gene_id, isoform_id,
        estimates, se, df, t, pval, FDR = regular_FDR, empirical_pval, empirical_FDR,
        p_type = p_type, p_used, gene_qvalue, stageR_gene, stageR_transcript
      )

    diagnostics <- make_diagnostics(fit, contrast_name, "ok",
      group_sizes = group_sizes_txt,
      covariates_used = paste(used_covariates, collapse = ";"),
      covariates_dropped = paste(setdiff(covariates, used_covariates), collapse = ";"),
      n_bins_input = n_bins_input,
      n_genes_tested = dplyr::n_distinct(results$gene_id),
      n_bins_tested = nrow(results),
      n_bins_valid_p = sum(!is.na(results$pval)),
      p_type = p_type,
      null_delta = null_est[["delta"]],
      null_sigma = null_est[["sigma"]],
      frac_raw_p05 = mean(results$pval < 0.05, na.rm = TRUE),
      frac_emp_p05 = mean(results$empirical_pval < 0.05, na.rm = TRUE),
      frac_used_p05 = mean(results$p_used < 0.05, na.rm = TRUE),
      n_sig_genes = dplyr::n_distinct(results$gene_id[!is.na(results$stageR_gene) & results$stageR_gene < alpha]),
      n_sig_bins = sum(!is.na(results$stageR_transcript) & results$stageR_transcript < alpha)
    )
    list(results = results, diagnostics = diagnostics)
  })

  minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  list(
    results = dplyr::bind_rows(purrr::map(per_contrast, "results")),
    diagnostics = dplyr::bind_rows(skipped_diag,
      dplyr::bind_rows(purrr::map(per_contrast, "diagnostics")) %>% dplyr::mutate(minutes = minutes))
  )
}

#----------------------------------------------------------------------------- #
## 2.7 Run all fits ----
fit_outputs <- BiocParallel::bplapply(seq_len(nrow(fit_grid)), run_fit, BPPARAM = bpparam)

## A fit that fails does not stop the others: its error goes to the diagnostics
ok <- BiocParallel::bpok(fit_outputs)
if (any(!ok)) {
  message(sum(!ok), " fit(s) failed: ", paste(which(!ok), collapse = ", "))
  fit_outputs[!ok] <- purrr::map(which(!ok), function(i) {
    list(results = NULL, diagnostics = make_diagnostics(fit_grid[i, ], contrast_tbl$contrast,
      paste0("error: ", conditionMessage(fit_outputs[[i]]))))
  })
}

############################################################################## #
# ---- 3. Outputs ----
results <- dplyr::bind_rows(purrr::map(fit_outputs, "results")) %>% tibble::as_tibble()
diagnostics <- dplyr::bind_rows(purrr::map(fit_outputs, "diagnostics"))

## Bin annotations of 07b (coordinates, end rank, merged transcripts)
if (nrow(results) > 0) {
  bin_annots <- as.data.frame(rowData(pseudobulk)) %>%
    tibble::as_tibble() %>%
    dplyr::select(-dplyr::any_of("gene_id")) %>%
    dplyr::rename(isoform_id = transcript_id)
  results <- results %>% dplyr::left_join(bin_annots, by = "isoform_id")
} else {
  warning("No fit was run; see the diagnostics for the reasons")
  results <- tibble::tibble(endome = character(0), stratum = character(0))   # header-only file
}

## The stratum column is named after the stratification column (e.g. region)
name_stratum <- function(df) {
  if (is.null(stratify_by)) dplyr::select(df, -stratum) else dplyr::rename(df, !!stratify_by := stratum)
}

readr::write_tsv(name_stratum(results), smk_outputs$results)
readr::write_tsv(name_stratum(diagnostics), smk_outputs$diagnostics)
saveRDS(name_stratum(results), smk_outputs$results_rds)
saveRDS(name_stratum(diagnostics), smk_outputs$diagnostics_rds)

message(sum(diagnostics$status == "ok"), " of ", nrow(diagnostics), " fit x contrast tests ran; ",
        sum(diagnostics$n_sig_genes, na.rm = TRUE), " significant genes in total (stageR, alpha = ", alpha, ")")
