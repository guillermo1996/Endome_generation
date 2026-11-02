#!/usr/bin/env Rscript

## _________________________________________________
##
## Process the scUTRquant output
##
## Aim: Match the nuclei of one scUTRquant SingleCellExperiment to the cell
## annotation table, compute their effective number of transcripts and flag
## the nuclei that pass QC
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S07a-Process_scUTRquant
## Date Created: 2026-10-06
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v3.1 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v3.1 (2026-10-28): Standardization of the script.
##    - v3.0 (2026-10-16): Cell annotation table read as TSV.
##    - v2.0 (2026-10-15): Dataset-independent cell annotation table.
##    - v1.0 (2026-10-06): Initial version.
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

  rs <- "results/Wood.control"
  endome_name <- "Wood.control.st_ref.pc_novel.w500.3p"
  pseudobulk_method <- "default"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      sce = file.path(rs, "06-UTR_Quantification/data/sce", endome_name, "hardy_snRNAseq.txs.Rds"),
      cell_annotations = "data/Hardy_snRNAseq/prepared/Hardy.cell_annotations.tsv.gz"
    ),
    output = list(
      cells = file.path(rs, "temp", sprintf("%s.%s.cells.tsv", endome_name, pseudobulk_method)),
      qc = file.path(rs, "temp", sprintf("%s.%s.cells_qc.tsv", endome_name, pseudobulk_method))
    ),
    params = list(
      exclude_samples = c("gosot_ACG", "gosot_IPL"),
      effective_txs_mode = "fixed",
      min_effective_txs = 300,
      mad_below_median = 3
    ),
    wildcards = list(
      endome_name = endome_name,
      pseudobulk_method = pseudobulk_method
    ),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(Matrix)
  library(SingleCellExperiment)
  library(conflicted)
  library(tidyverse)
})

## Package options
options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

## Solve possible conflicts between packages
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflict_prefer_all("tidyr", quiet = TRUE)

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
smk_inputs <- snakemake@input
smk_outputs <- snakemake@output
params <- snakemake@params
wc <- snakemake@wildcards

exclude_samples <- as.character(unlist(params$exclude_samples))
effective_txs_mode <- match.arg(params$effective_txs_mode, c("fixed", "mad"))
min_effective_txs <- params$min_effective_txs
mad_below_median <- params$mad_below_median

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Effective number of transcripts per nucleus
#'
#' 2^H, with H the Shannon entropy (bits) of the bin counts of the nucleus.
#' Same as compute_effective_count in the scUTRquant repository.
#'
#' @param cts Bins x nuclei count matrix.
#'
#' @returns Numeric vector with one value per nucleus (NaN for empty nuclei).
compute_effective_count <- function(cts) {
  cts <- as(cts, "CsparseMatrix")
  S <- Matrix::colSums(cts)
  xlogx <- cts
  xlogx@x <- cts@x * log2(cts@x)
  xlogx@x[is.nan(xlogx@x)] <- 0   # explicitly stored zeros (0 * -Inf)
  base::unname(2^(log2(S) - Matrix::colSums(xlogx) / S))
}

############################################################################## #
# ---- 1. Load the inputs ----

#----------------------------------------------------------------------------- #
## 1.1 Cell annotation table ----
message("Loading the cell annotations: ", smk_inputs$cell_annotations)
cell_annots <- readr::read_tsv(smk_inputs$cell_annotations, col_types = readr::cols(.default = readr::col_character()))

missing_cols <- setdiff(c("library_id", "barcode", "sample_id", "celltype"), colnames(cell_annots))
if (length(missing_cols) > 0) stop("Columns missing from the cell annotation table: ", paste(missing_cols, collapse = ", "), call. = FALSE)

### scUTRquant names each nucleus "<library_id>_<barcode>", without the "-1"
### suffix of the Cell Ranger barcodes
cell_annots <- cell_annots %>%
  dplyr::filter(!is.na(celltype), !sample_id %in% exclude_samples) %>%
  dplyr::mutate(barcode = sub("-\\d+$", "", barcode)) %>%
  dplyr::mutate(cell_id = paste0(library_id, "_", barcode), .before = 1)
annotation_cols <- setdiff(colnames(cell_annots), c("cell_id", "library_id", "barcode", "sample_id"))

#----------------------------------------------------------------------------- #
## 1.2 scUTRquant SCE ----
message("Loading the scUTRquant SCE: ", smk_inputs$sce)
sce <- readRDS(smk_inputs$sce)

### scUTRquant keeps every barcode above min_umis (mostly empty droplets). If
### no nucleus matches, the library_id or barcode format differs.
keep_cells <- colnames(sce) %in% cell_annots$cell_id
if (!any(keep_cells)) {
  stop("No annotated nucleus is in the SCE. SCE nucleus: '", colnames(sce)[1],
       "'; annotation table nucleus: '", cell_annots$cell_id[1], "'", call. = FALSE)
}
message("  ", sum(keep_cells), " of ", nrow(cell_annots), " annotated nuclei found in the SCE (",
        ncol(sce), " barcodes in the SCE)")

#----------------------------------------------------------------------------- #
## 1.3 Nucleus metrics ----
cts <- counts(sce)[, keep_cells, drop = FALSE]
cells <- tibble::tibble(
  cell_id = colnames(cts),
  total_umis = Matrix::colSums(cts),
  effective_txs = compute_effective_count(cts)
) %>%
  dplyr::inner_join(cell_annots, by = "cell_id", relationship = "one-to-one")
rm(sce, cts); invisible(gc())

############################################################################## #
# ---- 2. Effective-transcript QC ----

### "fixed": a single threshold. "mad": median - k * MAD of log(effective_txs)
### per cell type.
if (effective_txs_mode == "fixed") {
  cells <- cells %>%
    dplyr::mutate(effective_txs_threshold = min_effective_txs)
} else {
  cells <- cells %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(effective_txs_threshold = exp(
      median(log(effective_txs), na.rm = TRUE) - mad_below_median * mad(log(effective_txs), na.rm = TRUE))) %>%
    dplyr::ungroup()
}
cells <- cells %>%
  dplyr::mutate(pass_qc = !is.na(effective_txs) & effective_txs >= effective_txs_threshold)
message("  ", sum(cells$pass_qc), " of ", nrow(cells), " nuclei pass the QC (mode: ", effective_txs_mode, ")")

############################################################################## #
# ---- 3. Write the outputs ----

#----------------------------------------------------------------------------- #
## 3.1 Nucleus table ----
cells %>%
  dplyr::mutate(endome = wc$endome_name) %>%
  dplyr::select(endome, cell_id, library_id, barcode, sample_id, dplyr::all_of(annotation_cols),
                total_umis, effective_txs, effective_txs_threshold, pass_qc) %>%
  readr::write_tsv(smk_outputs$cells)

#----------------------------------------------------------------------------- #
## 3.2 QC summary per library and cell type ----
cells %>%
  dplyr::group_by(sample_id, library_id, celltype) %>%
  dplyr::summarise(
    n_annotated = dplyr::n(),
    n_pass = sum(pass_qc),
    frac_removed = 1 - n_pass / n_annotated,
    median_total_umis = median(total_umis),
    median_effective_txs = median(effective_txs, na.rm = TRUE),
    effective_txs_threshold = dplyr::first(effective_txs_threshold),
    .groups = "drop"
  ) %>%
  dplyr::mutate(endome = wc$endome_name, .before = 1) %>%
  readr::write_tsv(smk_outputs$qc)
