## _________________________________________________
##
## Differential transcript usage with satuRn
##
## Aim: Test differential usage of the ENDome bins between clinical groups
##      (e.g. PDD vs Control) within each region and cell type, from the
##      pseudobulked scUTRquant counts of one ENDome (06d).
##
## Project: ENDome generation - UTR Quantification
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-06
##
## Latest Version: v1.1 (2026-10-07)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - Follows Aine Fairbrother-Browne's DRIMSeq pipeline (04d; manuscript
##      methods "Transcript-level quantification and differential transcript
##      usage analysis"), replacing the DTU engine with satuRn:
##        - Same unit (sample_id = donor x region), same covariates, one model
##          per region x cell type, same DRIMSeq::dmFilter filter
##          (min_samps_gene_expr = 75% of the samples, min_samps_feature_expr =
##          smallest group, min_gene_expr = min_feature_expr = 10) and stageR
##          two-stage FDR control.
##        - satuRn::fitDTU (quasi-binomial GLM) + satuRn::testDTU replace
##          dmPrecision / dmFit / dmTest. Gene-level q-values for the stageR
##          screening stage come from DEXSeq's perGeneQValue, as in the satuRn
##          vignette (docs/satuRn/).
##    - fit_mode:
##        - "joint" (default): one fit per region x cell type with every group
##          in the contrasts, all contrasts tested from that fit (satuRn's
##          contrast system). All donors inform the dispersion and covariate
##          estimates.
##        - "pairwise": one fit per contrast with only its two groups, as in
##          Aine's pipeline.
##    - P-values: satuRn's empirical p-values (vignette recommendation). In
##      satuRn 1.18 they are only computed for contrasts with >= 500 tested
##      bins; below that, p_fallback = "raw" uses the raw p-values (p_type
##      column) and "none" skips the stageR step.
##    - Numeric covariates are centred and scaled (as in the manuscript's DE
##      model). This does not change the group contrasts. Covariates that are
##      constant within a fit are dropped and reported.
##    - Cell type vs cell type contrasts are not implemented yet; they need a
##      different design (donor as a blocking factor).
##
## Changelog:
##    - v1.1 (2026-10-07): Results carry the bin annotations of the pseudobulk
##      rowData (06d v1.1).
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
  ## dataset.group.merge_method folder created by the `test_data` rule (its
  ## name is the prefix wildcard), and the ENDome (orf_filter, width, txEnd) to
  ## test. The pseudobulk and nucleus counts are the 06d outputs linked by
  ## `test_data`. The outputs use the "interactive" prefix so they never
  ## overwrite the linked results.
  test_dir <- "data/test_data/Ebbert.control.st_ref"
  test_orf_filter <- "ref_pc"
  test_width <- "500"
  test_txEnd <- "3p"
  test_pseudobulk_method <- "default"
  test_endome <- paste0(test_orf_filter, ".w", test_width, ".", test_txEnd)
  test_pseudobulk <- paste0(test_endome, ".", test_pseudobulk_method)

  snakemake <- Snakemake(
    input = list(
      pseudobulk = file.path(test_dir, paste0("test.", test_pseudobulk, ".pseudobulk.rds")),
      n_cells = file.path(test_dir, paste0("test.", test_pseudobulk, ".pseudobulk_n_cells.tsv")),
      sample_covariates = "data/Hardy_snRNAseq/annotations/Hardy_final_covariates_weighted_and_combined_batches.csv"
    ),
    output = list(
      results = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.tsv.gz")),
      diagnostics = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.diagnostics.tsv")),
      diagplots = file.path(test_dir, paste0("interactive.", test_pseudobulk, ".satuRn.condition.diagplots"))
    ),
    params = list(
      exclude_samples = c("gosot"),
      regions = c("ACG", "IPL"),
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

exclude_samples <- params$exclude_samples
regions <- params$regions
celltypes <- params$celltypes
contrasts <- params$contrasts
covariates <- params$covariates
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

#' Empirical null estimates (delta, sigma) used by satuRn::testDTU
#'
#' Reproduces the estimation steps of satuRn:::p.adjust_empirical
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
#' perGeneQValue. Bins with NA p-values are ignored within a gene (the vignette
#' sets the whole gene to 1 if any bin is NA); genes with no valid p-value get 1.
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

#' One row per contrast describing a fit (or why it was skipped)
make_diagnostics <- function(fit, contrast_tbl, status, ...) {
  tibble::tibble(
    endome = wc$endome_name, region = fit$region, celltype = fit$celltype, fit = fit$fit,
    contrast = contrast_tbl$contrast, status = status, ...
  )
}

############################################################################## #
# ---- 1. Load inputs ----

#----------------------------------------------------------------------------- #
## 1.1 Pseudobulk and nuclei per pseudobulk (06d) ----
message("Loading the pseudobulk: ", smk_inputs$pseudobulk)
pseudobulk <- readRDS(smk_inputs$pseudobulk)
n_cells <- readr::read_tsv(smk_inputs$n_cells) %>% dplyr::select(sample_id, celltype, n_cells)

### Filter by celltypes if provided
if (is.null(celltypes)) celltypes <- assayNames(pseudobulk)
missing_celltypes <- setdiff(celltypes, assayNames(pseudobulk))
if (length(missing_celltypes) > 0) {
  warning("Cell types not in the pseudobulk (skipped): ", paste(missing_celltypes, collapse = ", "))
  celltypes <- intersect(celltypes, assayNames(pseudobulk))
}

#----------------------------------------------------------------------------- #
## 1.2 Sample covariates ----
contrast_tbl_all <- parse_contrasts(contrasts)
all_groups <- unique(c(contrast_tbl_all$group1, contrast_tbl_all$group2))

samples <- readr::read_csv(smk_inputs$sample_covariates) %>%
  dplyr::rename(sample_id = sample, group = group) %>%
  dplyr::mutate(region = stringr::str_extract(sample_id, "[^_]+$")) %>%
  dplyr::filter(group %in% all_groups, region %in% regions)

### Remove outlier/excluded samples
if (length(exclude_samples) > 0) {
  samples <- samples %>% dplyr::filter(!grepl(paste(exclude_samples, collapse = "|"), sample_id))
}

### Ensure that covariates are provided
missing_covariates <- setdiff(covariates, colnames(samples))
if (length(missing_covariates) > 0) stop("Covariates not in the covariate table: ", paste(missing_covariates, collapse = ", "))
if (anyNA(samples[, covariates])) stop("Missing values in the covariates; impute or drop them before 06e")

### Ensure samples match
not_in_pseudobulk <- setdiff(samples$sample_id, colnames(pseudobulk))
if (length(not_in_pseudobulk) > 0) {
  message("Samples in the covariate table but not in the pseudobulk (ignored): ", paste(not_in_pseudobulk, collapse = ", "))
}
samples <- samples %>% dplyr::filter(sample_id %in% colnames(pseudobulk))
message("Samples: ", nrow(samples), " (", paste(names(table(samples$region)), table(samples$region), sep = "=", collapse = ", "), ")")

#----------------------------------------------------------------------------- #
## 1.3 Fits to run ----
fit_grid <- tidyr::expand_grid(region = regions, celltype = celltypes) %>% 
  dplyr::mutate(fit = "joint", contrasts = list(contrast_tbl_all))
message("Fits to run: ", nrow(fit_grid))

############################################################################## #
# ---- 2. satuRn DTU ----

#' Run one satuRn fit (one region x cell type [x contrast])
#'
#' @param i Row of fit_grid
#' @return A list with `results` (tibble or NULL) and `diagnostics` (tibble)
run_fit <- function(i) {
  # i <- 1
  start_time <- Sys.time()
  fit <- fit_grid[i, ]
  contrast_tbl <- fit$contrasts[[1]]
  fit_groups <- unique(c(contrast_tbl$group1, contrast_tbl$group2))

  #----------------------------------------------------------------------------- #
  ## 2.1 Samples of this fit ----
  fit_samples <- samples %>%
    dplyr::filter(region == fit$region, group %in% fit_groups) %>%
    dplyr::inner_join(n_cells %>% dplyr::filter(celltype == fit$celltype) %>% dplyr::select(-celltype), by = "sample_id") %>%
    dplyr::filter(n_cells >= min_cells)

  cts <- assay(pseudobulk, fit$celltype)[, fit_samples$sample_id, drop = FALSE]
  fit_samples <- fit_samples %>% dplyr::filter(colSums(cts) > 0)
  cts <- cts[, fit_samples$sample_id, drop = FALSE]

  ## Drop groups with too few samples
  group_sizes <- table(factor(fit_samples$group, levels = fit_groups))
  kept_groups <- names(group_sizes)[group_sizes >= min_samples_per_group]
  fit_samples <- fit_samples %>% dplyr::filter(group %in% kept_groups)
  cts <- cts[, fit_samples$sample_id, drop = FALSE]
  group_sizes <- group_sizes[kept_groups]
  group_sizes_txt <- paste(names(group_sizes), group_sizes, sep = "=", collapse = ";")

  testable <- contrast_tbl$group1 %in% kept_groups & contrast_tbl$group2 %in% kept_groups
  skipped_diag <- if (any(!testable)) {
    make_diagnostics(fit, contrast_tbl[!testable, ], paste0("skipped: < ", min_samples_per_group, " samples with >= ", min_cells, " nuclei in a group"),
                     group_sizes = group_sizes_txt)
  } else NULL
  contrast_tbl <- contrast_tbl[testable, ]
  if (nrow(contrast_tbl) == 0) return(list(results = NULL, diagnostics = skipped_diag))

  #----------------------------------------------------------------------------- #
  ## 2.2 Filter bins with DRIMSeq::dmFilter ----
  n_bins_input <- sum(rowSums(cts) > 0)
  d <- DRIMSeq::dmDSdata(
    counts = data.frame(gene_id = rowData(pseudobulk)$gene_id, feature_id = rownames(cts),
                        as.matrix(cts), check.names = FALSE, stringsAsFactors = FALSE),
    samples = data.frame(sample_id = fit_samples$sample_id, group = fit_samples$group, stringsAsFactors = FALSE)
  )
  d <- tryCatch(
    DRIMSeq::dmFilter(d,
      min_samps_gene_expr = round(nrow(fit_samples) * dmfilter$min_samps_gene_prop),
      min_samps_feature_expr = min(group_sizes),
      min_gene_expr = dmfilter$min_gene_expr,
      min_feature_expr = dmfilter$min_feature_expr),
    error = function(e) NULL)

  filt <- if (is.null(d)) NULL else DRIMSeq::counts(d)
  if (!is.null(filt)) {
    ## satuRn needs >= 2 bins per gene
    filt <- filt %>%
      dplyr::mutate(gene_id = as.character(gene_id), feature_id = as.character(feature_id)) %>%
      dplyr::group_by(gene_id) %>% dplyr::filter(dplyr::n() >= 2) %>% dplyr::ungroup()
  }
  if (is.null(filt) || nrow(filt) == 0) {
    return(list(results = NULL, diagnostics = dplyr::bind_rows(skipped_diag,
      make_diagnostics(fit, contrast_tbl, "skipped: no genes left after filtering",
                       group_sizes = group_sizes_txt, n_bins_input = n_bins_input))))
  }
  n_genes_filt <- dplyr::n_distinct(filt$gene_id)
  if (n_genes_filt < min_genes) {
    return(list(results = NULL, diagnostics = dplyr::bind_rows(skipped_diag,
      make_diagnostics(fit, contrast_tbl,
                      paste0("skipped: ", n_genes_filt, " genes after filtering (< min_genes = ", min_genes, ")"),
                      group_sizes = group_sizes_txt, n_bins_input = n_bins_input))))
  }

  #----------------------------------------------------------------------------- #
  ## 2.3 Design matrix ----
  coldata <- fit_samples %>% dplyr::mutate(group = factor(group, levels = kept_groups))
  used_covariates <- character(0)
  for (cov in covariates) {
    values <- coldata[[cov]]
    if (dplyr::n_distinct(values) < 2) next
    if (is.numeric(values)) {
      coldata[[cov]] <- as.numeric(scale(values))
    } else {
      coldata[[cov]] <- factor(values)
    }
    used_covariates <- c(used_covariates, cov)
  }
  dropped_covariates <- setdiff(covariates, used_covariates)

  design_formula <- as.formula(paste("~ 0 + group", paste(c("", used_covariates), collapse = " + ")))
  coldata <- as.data.frame(coldata)
  rownames(coldata) <- coldata$sample_id
  design <- model.matrix(design_formula, data = coldata)
  if (qr(design)$rank < ncol(design)) {
    return(list(results = NULL, diagnostics = dplyr::bind_rows(skipped_diag,
      make_diagnostics(fit, contrast_tbl, "skipped: rank-deficient design", group_sizes = group_sizes_txt))))
  }

  #----------------------------------------------------------------------------- #
  ## 2.4 Fit satuRn ----
  count_mat <- as.matrix(filt[, fit_samples$sample_id])
  rownames(count_mat) <- filt$feature_id
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
    contrasts = paste0("group", contrast_tbl$group1, "-group", contrast_tbl$group2),
    levels = colnames(design)
  )
  colnames(L) <- contrast_tbl$contrast

  ## satuRn 1.18 draws its diagnostic plots regardless of diagplot1/diagplot2
  fit_id <- safe_name(paste(fit$region, fit$celltype, fit$fit, sep = "."))
  grDevices::pdf(file.path(smk_outputs$diagplots, paste0(fit_id, ".pdf")), width = 8, height = 6)
  graphics::plot.new()
  graphics::text(0.5, 0.5, paste0(wc$endome_name, "\n", fit$region, " | ", fit$celltype, " | ", fit$fit,
                                  "\n", group_sizes_txt), cex = 1.1)
  se <- tryCatch(
    satuRn::testDTU(object = se, contrasts = L, diagplot1 = TRUE, diagplot2 = TRUE, sort = FALSE),
    finally = grDevices::dev.off())

  #----------------------------------------------------------------------------- #
  ## 2.6 Per contrast: p-values, gene q-values and stageR ----
  per_contrast <- purrr::map(contrast_tbl$contrast, function(contrast_name) {
    # contrast_name <- contrast_tbl$contrast[[1]]
    res <- as.data.frame(rowData(se)[[paste0("fitDTUResult_", contrast_name)]])
    res$isoform_id <- rowData(se)$isoform_id
    res$gene_id <- rowData(se)$gene_id

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
        endome = wc$endome_name, region = fit$region, celltype = fit$celltype, fit = fit$fit,
        contrast = contrast_name, gene_id, isoform_id,
        estimates, se, df, t, pval, FDR = regular_FDR, empirical_pval, empirical_FDR,
        p_type = p_type, p_used, gene_qvalue, stageR_gene, stageR_transcript
      )

    diagnostics <- make_diagnostics(fit, contrast_tbl[contrast_tbl$contrast == contrast_name, ], "ok",
      group_sizes = group_sizes_txt,
      covariates_used = paste(used_covariates, collapse = ";"),
      covariates_dropped = paste(dropped_covariates, collapse = ";"),
      n_bins_input = n_bins_input,
      n_genes_tested = dplyr::n_distinct(results$gene_id),
      n_bins_tested = nrow(results),
      n_bins_valid_p = sum(!is.na(results$pval)),
      p_type = p_type,
      null_delta = null_est[["delta"]],
      null_sigma = null_est[["sigma"]],
      frac_raw_p05 = mean(results$pval < 0.05, na.rm = TRUE),
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

## Fits that failed unexpectedly are reported in the diagnostics, not dropped
ok <- BiocParallel::bpok(fit_outputs)
if (any(!ok)) {
  message(sum(!ok), " fit(s) failed: ", paste(which(!ok), collapse = ", "))
  fit_outputs[!ok] <- purrr::map(which(!ok), function(i) {
    list(results = NULL, diagnostics = make_diagnostics(fit_grid[i, ], fit_grid$contrasts[[i]],
      paste0("error: ", conditionMessage(fit_outputs[[i]]))))
  })
}

############################################################################## #
# ---- 3. Outputs ----
results <- dplyr::bind_rows(purrr::map(fit_outputs, "results")) %>% tibble::as_tibble()
diagnostics <- dplyr::bind_rows(purrr::map(fit_outputs, "diagnostics"))

## Bin annotations from the pseudobulk rowData (06d: ENDome coordinates, end
## rank, merged transcripts), when present
bin_annots <- as.data.frame(rowData(pseudobulk)) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("gene_id")) %>%
  dplyr::rename(isoform_id = transcript_id)
if (ncol(bin_annots) > 1 && nrow(results) > 0) {
  results <- results %>% dplyr::left_join(bin_annots, by = "isoform_id")
}

readr::write_tsv(results, smk_outputs$results)
readr::write_tsv(diagnostics, smk_outputs$diagnostics)

message("Done. ", sum(diagnostics$status == "ok"), " of ", nrow(diagnostics), " fit x contrast tests ran; ",
        sum(diagnostics$n_sig_genes, na.rm = TRUE), " significant genes in total (stageR, alpha = ", alpha, ")")
