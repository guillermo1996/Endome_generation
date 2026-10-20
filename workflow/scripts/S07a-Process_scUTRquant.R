## _________________________________________________
##
## Process scUTRquant output: nucleus annotation and QC
##
## Aim: Match the nuclei of one scUTRquant SingleCellExperiment (one ENDome) to
##      the cell annotation table, compute the effective number of transcripts
##      of each nucleus and flag the nuclei that pass QC. The output is a
##      nucleus table (no counts) that 07b uses to build the pseudobulk.
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
##    - QC (effective_txs_mode):
##        - "fixed": keep nuclei with effective_txs >= min_effective_txs.
##        - "mad": keep nuclei above median - mad_below_median * MAD of
##          log(effective_txs), computed per cell type.
##
## Changelog:
##    - v3.0 (2026-10-16): Simplified. Reads the cell annotation table as TSV.
##    - v2.0 (2026-10-15): Dataset-independent cell annotation table.
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
  test_dir <- "data/test_data/Ebbert.control.st_ref"
  test_endome <- "ref_pc.w500.3p"

  snakemake <- Snakemake(
    input = list(
      sce = file.path(test_dir, paste0("test.", test_endome, ".hardy_snRNAseq.txs.Rds")),
      cell_annotations = "data/Hardy_snRNAseq/prepared/Hardy.cell_annotations.tsv.gz"
    ),
    output = list(
      cells = file.path(test_dir, paste0("interactive.", test_endome, ".cells.tsv")),
      qc = file.path(test_dir, paste0("interactive.", test_endome, ".cells_qc.tsv"))
    ),
    params = list(
      exclude_samples = c("gosot_ACG", "gosot_IPL"),
      effective_txs_mode = "fixed",
      min_effective_txs = 300,
      mad_below_median = 3,
      cell_type_level = "annotation_level_2"
    ),
    wildcards = list(
      endome_name = paste0(basename(test_dir), ".", test_endome)
    ),
    log = list(),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages

shhh({
  library(conflicted)
  library(tidyverse)
  library(Matrix)
  library(SingleCellExperiment)
})

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

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
cell_type_level <- params$cell_type_level

dir.create(dirname(smk_outputs$cells), showWarnings = FALSE, recursive = TRUE)

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Effective number of transcripts per nucleus
#'
#' 2^H, with H the Shannon entropy (bits) of the nucleus' bin counts. Same
#' formula as compute_effective_count in the scUTRquant repository.
#'
#' @param cts A bins x nuclei count matrix
#' @return A numeric vector, one value per nucleus (NaN for empty nuclei)
compute_effective_count <- function(cts) {
  cts <- as(cts, "CsparseMatrix")
  S <- Matrix::colSums(cts)
  xlogx <- cts
  xlogx@x <- cts@x * log2(cts@x)
  xlogx@x[is.nan(xlogx@x)] <- 0   # explicitly stored zeros (0 * -Inf)
  base::unname(2^(log2(S) - Matrix::colSums(xlogx) / S))
}

############################################################################## #
# ---- 1. Load inputs ----

#----------------------------------------------------------------------------- #
## 1.1 Cell annotation table ----
message("Loading the cell annotations: ", smk_inputs$cell_annotations)
cell_annots <- readr::read_tsv(smk_inputs$cell_annotations, col_types = readr::cols(.default = readr::col_character()))

missing_cols <- setdiff(c("library_id", "barcode", "sample_id", "celltype"), colnames(cell_annots))
if (length(missing_cols) > 0) stop("Columns missing from the cell annotation table: ", paste(missing_cols, collapse = ", "))

## scUTRquant names each nucleus "<library_id>_<barcode>". Cell Ranger barcodes
## carry a "-1" suffix that scUTRquant does not use, so it is removed.
cell_annots <- cell_annots %>%
  dplyr::filter(!is.na(celltype), !sample_id %in% exclude_samples) %>%
  dplyr::mutate(barcode = sub("-\\d+$", "", barcode)) %>%
  dplyr::mutate(cell_id = paste0(library_id, "_", barcode), .before = 1)
annotation_cols <- setdiff(colnames(cell_annots), c("cell_id", "library_id", "barcode", "sample_id"))

#----------------------------------------------------------------------------- #
## 1.2 scUTRquant SCE ----
message("Loading the scUTRquant SCE: ", smk_inputs$sce)
sce <- readRDS(smk_inputs$sce)

## scUTRquant keeps every barcode above min_umis (mostly empty droplets): only
## the annotated nuclei are kept. If none match, the library_id or barcode
## format of the annotation table differs from scUTRquant's.
keep_cells <- colnames(sce) %in% cell_annots$cell_id
if (!any(keep_cells)) {
  stop("No annotated nucleus is in the SCE. SCE nucleus: '", colnames(sce)[1],
       "'; annotation table nucleus: '", cell_annots$cell_id[1], "'")
}
message("  ", sum(keep_cells), " of ", nrow(cell_annots), " annotated nuclei found in the SCE (",
        ncol(sce), " barcodes in the SCE)")

#----------------------------------------------------------------------------- #
## 1.3 Per-nucleus metrics ----
cts <- counts(sce)[, keep_cells, drop = FALSE]
cells <- tibble::tibble(
  cell_id = colnames(cts),
  total_umis = Matrix::colSums(cts),
  effective_txs = compute_effective_count(cts)
) %>%
  ## one-to-one: fails if a nucleus is annotated twice
  dplyr::inner_join(cell_annots, by = "cell_id", relationship = "one-to-one")
rm(sce, cts); invisible(gc())

############################################################################## #
# ---- 2. Effective-transcript QC ----
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
# ---- 3. Outputs ----

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
