## _________________________________________________
##
## Process scUTRquant output: nucleus annotation and QC
##
## Aim: Annotate the nuclei of one scUTRquant SingleCellExperiment (one ENDome)
##      with the Hardy cell-type annotations, compute the per-nucleus effective
##      number of transcripts, and flag the nuclei that pass QC. The output is a
##      nucleus table (no counts) that 06d uses to build the processed SCE and
##      the pseudobulk.
##
## Project: ENDome generation - UTR Quantification
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-06
##
## Latest Version: v1.0 (2026-10-06)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - Follows Aine Fairbrother-Browne's processing (04c-process_output_scUTRquant.R
##      and section 2 of 04d-run_DRIMSeq_DUTR_ALLcovs.R; 41/Data/Aine-DRIMSeq):
##        1. Effective number of transcripts per nucleus (Shannon entropy,
##           compute_effective_count from the scUTRquant GitHub repository).
##        2. scUTRquant's annotation columns are dropped and Melissa's annotations
##           (annotation_all_levels_all_nuclei_updated.csv, the cell-type names
##           used in the manuscript; the older file only differs in the names of
##           3 astrocyte types at level 2 and in coarser glial states at level 3)
##           are re-joined with the sample -> S-number mapping from SAMPLE.csv,
##           including Aine's library-label fixes (borah_IPL, nobin_ACG and the
##           zupam_IPL B1/B2 swap). This mapping is kept unmodified on purpose:
##           docs/scUTRquant_cell_selection_and_hardy_annotations.md describes
##           the remaining B1/B2 mismatches (6 pairs), to be revisited separately.
##        3. Nuclei without annotation (or labelled "drop") and excluded samples
##           (the outlier donor, gosot) are removed.
##    - Differences with Aine's processing:
##        - The effective-transcript threshold is a parameter. "fixed" keeps
##          effective_txs >= min_effective_txs (manuscript: 300). "mad" keeps
##          nuclei above median - nmads * MAD of log(effective_txs) within each
##          cell type, for ENDomes whose bin sets make a fixed threshold
##          incomparable (e.g. 5' vs 3' ENDomes).
##        - No size factors or normalised counts: 06e models within-sample
##          proportions with satuRn, which takes raw counts.
##        - No count matrix is written. The output is the nucleus table only,
##          so the multi-GB SCE is not duplicated on disk.
##
## Changelog:
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
  test_orf_filter <- "ref_pc"
  test_width <- "500"
  test_txEnd <- "3p"
  test_dataset_name <- "hardy_snRNAseq"
  test_endome <- paste0(test_orf_filter, ".w", test_width, ".", test_txEnd, ".", test_dataset_name)

  snakemake <- Snakemake(
    input = list(
      sce = file.path(test_dir, paste0("test.", test_endome, ".txs.Rds")),
      annotations = "data/Hardy_snRNAseq/annotations/annotation_all_levels_all_nuclei_updated.csv",
      sample_metadata = "data/Hardy_snRNAseq/annotations/SAMPLE.csv"
    ),
    output = list(
      cells = file.path(test_dir, paste0("interactive.", test_endome, ".cells.tsv")),
      qc = file.path(test_dir, paste0("interactive.", test_endome, ".cells_qc.tsv"))
    ),
    params = list(
      exclude_samples = c("gosot"),
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

exclude_samples <- params$exclude_samples
effective_txs_mode <- match.arg(params$effective_txs_mode, c("fixed", "mad"))
min_effective_txs <- params$min_effective_txs
mad_below_median <- params$mad_below_median
cell_type_level <- params$cell_type_level

dir.create(dirname(smk_outputs$cells), showWarnings = FALSE, recursive = TRUE)

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Effective number of transcripts per nucleus
#'
#' 2^H, with H the Shannon entropy (bits) of the nucleus' transcript counts.
#' Same formula as compute_effective_count in the scUTRquant GitHub repository
#' (used by Aine as computeEffectiveTxCount in 04c).
#'
#' @param cts A features x nuclei count matrix (dgCMatrix)
#' @return A numeric vector, one value per nucleus (NaN for empty nuclei)
compute_effective_count <- function(cts) {
  cts <- as(cts, "CsparseMatrix")
  S <- Matrix::colSums(cts)
  xlogx <- cts
  xlogx@x <- cts@x * log2(cts@x)
  xlogx@x[is.nan(xlogx@x)] <- 0   # explicitly stored zeros (0 * -Inf)
  base::unname(2^(log2(S) - Matrix::colSums(xlogx) / S))
}

#' Melissa's nucleus annotations, with Aine's library-label fixes
#'
#' Index format: "<barcode>-<sample_id>[B1|B2]". Copied from 04d (section 2).
#'
#' @param path Path to annotation_all_levels_all_nuclei[_updated].csv
#' @param exclude_samples Sample stubs to remove (e.g. the outlier donor)
#' @return A tibble with index, sample_id_batch and annotation_level_*
read_annotations <- function(path, exclude_samples) {
  annots <- readr::read_csv(path)
  if (length(exclude_samples) > 0) {
    annots <- annots %>% dplyr::filter(!grepl(paste(exclude_samples, collapse = "|"), index))
  }
  annots %>%
    dplyr::filter(annotation_level_0 != "drop") %>%
    tidyr::separate(index, into = c("barcode", "sample_id_batch"), sep = "-", remove = FALSE) %>%
    dplyr::select(-barcode) %>%
    ## AFB: the unreplicated libraries of these two samples are their B2 library
    dplyr::mutate(
      sample_id_batch = dplyr::case_when(
        sample_id_batch == "borah_IPL" ~ "borah_IPLB2",
        sample_id_batch == "nobin_ACG" ~ "nobin_ACGB2",
        TRUE ~ sample_id_batch),
      index = dplyr::case_when(
        sample_id_batch %in% c("borah_IPLB2", "nobin_ACGB2") ~ paste0(index, "B2"),
        TRUE ~ index)
    ) %>%
    ## AFB: zupam_IPL B1/B2 labels were swapped
    dplyr::mutate(
      index = dplyr::case_when(
        grepl("zupam_IPLB1", index) ~ gsub("B1", "B2", index),
        grepl("zupam_IPLB2", index) ~ gsub("B2", "B1", index),
        TRUE ~ index),
      sample_id_batch = dplyr::case_when(
        grepl("zupam_IPLB1", sample_id_batch) ~ gsub("B1", "B2", sample_id_batch),
        grepl("zupam_IPLB2", sample_id_batch) ~ gsub("B2", "B1", sample_id_batch),
        TRUE ~ sample_id_batch)
    )
}

#' Library (S-number) -> sample_id / sample_id_batch mapping from SAMPLE.csv
#'
#' Copied from 04d (section 2): sample_id_batch = sample_id + "B<replicate>" for
#' replicated samples, sample_id otherwise.
#'
#' @param path Path to the CRN SAMPLE.csv
#' @return A tibble with sample_id, sample_id_batch and sample_id_snum
read_library_mapping <- function(path) {
  readr::read_csv(path) %>%
    dplyr::select(sample_id, replicate, file_name) %>%
    dplyr::mutate(s_number = stringr::str_match(file_name, "_(S\\d+)_")[, 2]) %>%
    dplyr::filter(!is.na(s_number)) %>%
    dplyr::distinct(sample_id, replicate, s_number) %>%
    dplyr::mutate(
      sample_id_snum = paste0(sample_id, "_", s_number),
      sample_id_batch = dplyr::if_else(!is.na(replicate), paste0(sample_id, "B", gsub("rep", "", replicate)), sample_id)
    ) %>%
    dplyr::distinct(sample_id, sample_id_batch, sample_id_snum)
}

############################################################################## #
# ---- 1. Load scUTRquant output ----
message("Loading the scUTRquant SCE: ", smk_inputs$sce)
sce <- readRDS(smk_inputs$sce)
message("  ", nrow(sce), " bins x ", ncol(sce), " nuclei")

#----------------------------------------------------------------------------- #
## 1.1 Per-nucleus metrics ----
message("Computing total UMIs and effective number of transcripts...")
cells <- tibble::tibble(
  cell_id = colnames(sce),
  total_umis = Matrix::colSums(counts(sce)),
  effective_txs = compute_effective_count(counts(sce))
)

## cell_id = "<sample_id>_S<n>_<barcode>", e.g. babom_ACG_S20_AAACCCAAGACGGAAA
cells <- cells %>%
  tidyr::extract(cell_id, into = c("library_sample", "s_num", "barcode"),
                 regex = "^(.+)_S(\\d+)_([ACGTN]+)$", remove = FALSE) %>%
  dplyr::mutate(sample_id_snum = paste0(library_sample, "_S", s_num)) %>%
  dplyr::select(-library_sample, -s_num)

n_unparsed <- sum(is.na(cells$barcode))
if (n_unparsed > 0) warning(n_unparsed, " nuclei with an unexpected cell_id format were dropped")
cells <- cells %>% dplyr::filter(!is.na(barcode))

rm(sce); invisible(gc())

############################################################################## #
# ---- 2. Annotate nuclei (Aine's mapping) ----

#----------------------------------------------------------------------------- #
## 2.1 Library mapping and annotations ----
message("Joining the library mapping and Melissa's annotations...")
library_mapping <- read_library_mapping(smk_inputs$sample_metadata)
annotations <- read_annotations(smk_inputs$annotations, exclude_samples)

cells <- cells %>%
  dplyr::left_join(library_mapping, by = "sample_id_snum") %>%
  dplyr::mutate(index = paste0(barcode, "-", sample_id_batch)) %>%
  dplyr::left_join(annotations, by = c("index", "sample_id_batch"))

#----------------------------------------------------------------------------- #
## 2.2 Keep annotated nuclei ----
annotation_cols <- grep("^annotation_level_", colnames(cells), value = TRUE)
n_libraries_unmapped <- cells %>% dplyr::filter(is.na(sample_id)) %>% dplyr::distinct(sample_id_snum) %>% nrow()
if (n_libraries_unmapped > 0) warning(n_libraries_unmapped, " libraries are missing from SAMPLE.csv")

cells <- cells %>%
  dplyr::filter(!is.na(sample_id), dplyr::if_all(dplyr::all_of(annotation_cols), ~ !is.na(.x)))
if (length(exclude_samples) > 0) {
  cells <- cells %>% dplyr::filter(!grepl(paste(exclude_samples, collapse = "|"), sample_id))
}
message("  ", nrow(cells), " annotated nuclei in ", dplyr::n_distinct(cells$sample_id), " samples (",
        dplyr::n_distinct(cells$sample_id_batch), " libraries)")

############################################################################## #
# ---- 3. Effective-transcript QC ----
if (effective_txs_mode == "fixed") {
  cells <- cells %>%
    dplyr::mutate(effective_txs_threshold = min_effective_txs)
} else {
  cells <- cells %>%
    dplyr::group_by(.data[[cell_type_level]]) %>%
    dplyr::mutate(effective_txs_threshold = exp(
      median(log(effective_txs), na.rm = TRUE) - mad_below_median * mad(log(effective_txs), na.rm = TRUE))) %>%
    dplyr::ungroup()
}
cells <- cells %>%
  dplyr::mutate(pass_qc = !is.na(effective_txs) & effective_txs >= effective_txs_threshold)
message("  ", sum(cells$pass_qc), " nuclei pass the effective-transcript QC (mode: ", effective_txs_mode, ")")

############################################################################## #
# ---- 4. Outputs ----

#----------------------------------------------------------------------------- #
## 4.1 Nucleus table ----
cells %>%
  dplyr::mutate(endome = wc$endome_name) %>%
  dplyr::select(endome, cell_id, sample_id, sample_id_batch, barcode, dplyr::all_of(annotation_cols),
                total_umis, effective_txs, effective_txs_threshold, pass_qc) %>%
  readr::write_tsv(smk_outputs$cells)

#----------------------------------------------------------------------------- #
## 4.2 QC summary per library and cell type ----
## Fraction of annotated nuclei removed by the effective-transcript threshold,
## to compare thresholds across cell types and ENDomes.
cells %>%
  dplyr::group_by(sample_id, sample_id_batch, celltype = .data[[cell_type_level]]) %>%
  dplyr::summarise(
    n_annotated = dplyr::n(),
    n_pass = sum(pass_qc),
    frac_removed = 1 - n_pass / n_annotated,
    median_total_umis = median(total_umis),
    median_effective_txs = median(effective_txs, na.rm = TRUE),
    effective_txs_threshold = dplyr::first(effective_txs_threshold),
    .groups = "drop"
  ) %>%
  dplyr::mutate(endome = wc$endome_name, annotation_level = cell_type_level, .before = 1) %>%
  readr::write_tsv(smk_outputs$qc)

