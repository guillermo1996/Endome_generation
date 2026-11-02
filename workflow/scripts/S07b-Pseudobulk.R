#!/usr/bin/env Rscript

## _________________________________________________
##
## Processed scUTRquant SCE and pseudobulk
##
## Aim: Build the processed SingleCellExperiment of one ENDome (QC-passing
## nuclei, bins annotated with their ENDome coordinates) and sum its raw counts
## by sample and cell type
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S07b-Pseudobulk
## Date Created: 2026-10-06
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v2.2 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - No minimum number of nuclei is applied here. S07c applies min_cells
##      with the n_cells table.
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v2.2 (2026-10-28): Standardization of the script.
##    - v2.1 (2026-10-20): is_mt also matches the Ensembl "MT".
##    - v2.0 (2026-10-16): The ENDome GTF and merge table are required.
##    - v1.1 (2026-10-07): Processed SCE output, bin annotations and removal of
##      bins without counts.
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
      cells = file.path(rs, "07-UTR_Usage/cells", sprintf("%s.%s.cells.tsv", endome_name, pseudobulk_method)),
      gtf = file.path(rs, "05-Truncation/txendcutr", sprintf("%s.txendcutr.gtf", endome_name)),
      merge_tsv = file.path(rs, "05-Truncation/txendcutr", sprintf("%s.txendcutr.merge.tsv", endome_name))
    ),
    output = list(
      sce = file.path(rs, "temp", sprintf("%s.%s.processed.txs.rds", endome_name, pseudobulk_method)),
      pseudobulk = file.path(rs, "temp", sprintf("%s.%s.pseudobulk.rds", endome_name, pseudobulk_method)),
      n_cells = file.path(rs, "temp", sprintf("%s.%s.pseudobulk_n_cells.tsv", endome_name, pseudobulk_method))
    ),
    params = list(),
    wildcards = list(
      endome_name = endome_name,
      pseudobulk_method = pseudobulk_method
    ),
    threads = 4,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(Matrix)
  library(SingleCellExperiment)
  library(dreamlet)
  library(BiocParallel)
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
wc <- snakemake@wildcards

bpparam <- BiocParallel::MulticoreParam(workers = snakemake@threads)

### The ENDome name ends in .w{width}.{txEnd}
tx_end <- stringr::str_extract(wc$endome_name, "(3p|5p)$")

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Bin coordinates from the ENDome GTF
#'
#' @param path Path to the ENDome GTF.
#' @param tx_end "3p" or "5p".
#'
#' @returns Tibble with `transcript_id`, `bin_chrom`, `bin_start`, `bin_end`,
#'   `bin_strand` and `end_coord`, the position of the end that defines the bin
#'   (3' end for 3p ENDomes, 5' end for 5p ENDomes).
readBinCoordinates <- function(path, tx_end) {
  rtracklayer::import(path, format = "gtf", feature.type = "transcript") %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      transcript_id,
      bin_chrom = as.character(seqnames),
      bin_start = start,
      bin_end = end,
      bin_strand = as.character(strand),
      end_coord = dplyr::case_when(
        tx_end == "3p" & bin_strand == "-" ~ bin_start,
        tx_end == "3p" ~ bin_end,
        tx_end == "5p" & bin_strand == "-" ~ bin_end,
        TRUE ~ bin_start)
    ) %>%
    dplyr::distinct(transcript_id, .keep_all = TRUE)
}

#' Transcripts collapsed into each bin
#'
#' @param path Path to the txendcutr merge table (tx_in, tx_out, gene_out).
#'
#' @returns Tibble with `transcript_id` (= tx_out), `n_merged_txs` and `merged_txs`.
readBinMembers <- function(path) {
  readr::read_tsv(path) %>%
    dplyr::group_by(transcript_id = tx_out) %>%
    dplyr::summarise(n_merged_txs = dplyr::n(), merged_txs = paste(sort(tx_in), collapse = ";"), .groups = "drop")
}

############################################################################## #
# ---- 1. Load the inputs ----
cells <- readr::read_tsv(smk_inputs$cells, col_types = readr::cols(
    .default = readr::col_character(),
    total_umis = readr::col_double(),
    effective_txs = readr::col_double(),
    effective_txs_threshold = readr::col_double(),
    pass_qc = readr::col_logical())) %>%
  dplyr::filter(pass_qc)
message("QC-passing nuclei: ", nrow(cells), " in ", dplyr::n_distinct(cells$sample_id), " samples")

message("Loading the scUTRquant SCE: ", smk_inputs$sce)
sce <- readRDS(smk_inputs$sce)

############################################################################## #
# ---- 2. Processed SCE ----

#----------------------------------------------------------------------------- #
## 2.1 QC-passing nuclei ----
sce <- sce[, cells$cell_id]
colData(sce) <- S4Vectors::DataFrame(
  cells %>% dplyr::select(-endome, -pass_qc) %>% as.data.frame(),
  row.names = cells$cell_id
)

#----------------------------------------------------------------------------- #
## 2.2 Bins with counts ----
keep_bins <- Matrix::rowSums(counts(sce)) > 0
message("Bins with counts in the retained nuclei: ", sum(keep_bins), " of ", length(keep_bins))
sce <- sce[keep_bins, ]

#----------------------------------------------------------------------------- #
## 2.3 Bin annotations ----

### end_rank orders the bins of a gene by the position of their end, in the
### direction of transcription (1 = most upstream). The coordinates are named
### bin_* because seqnames/start/end/strand clash with SummarizedExperiment.
message("Adding the bin annotations from the ENDome GTF and merge table...")
bin_annots <- tibble::tibble(
  transcript_id = rowData(sce)$transcript_id,
  gene_id = rowData(sce)$gene_id
) %>%
  dplyr::left_join(readBinCoordinates(smk_inputs$gtf, tx_end), by = "transcript_id") %>%
  dplyr::mutate(signed_end = ifelse(bin_strand == "-", -1, 1) * end_coord) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(
    n_bins_gene = dplyr::n(),
    end_rank = rank(signed_end, ties.method = "min", na.last = "keep"),
    is_most_upstream_end = end_rank == min(end_rank, na.rm = TRUE),
    is_most_downstream_end = end_rank == max(end_rank, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-signed_end) %>%
  dplyr::mutate(is_mt = bin_chrom %in% c("chrM", "MT")) %>%
  dplyr::left_join(readBinMembers(smk_inputs$merge_tsv), by = "transcript_id")

rowData(sce) <- S4Vectors::DataFrame(as.data.frame(bin_annots), row.names = rownames(sce))
S4Vectors::metadata(sce)$endome <- wc$endome_name
S4Vectors::metadata(sce)$tx_end <- tx_end

############################################################################## #
# ---- 3. Pseudobulk ----

### Only sample_id and celltype are kept in the colData, so that no nucleus
### column ends up in the pseudobulk. One assay per cell type, one column per
### sample.
sce_pb <- sce
colData(sce_pb) <- S4Vectors::DataFrame(
  sample_id = sce$sample_id,
  celltype = sce$celltype,
  row.names = colnames(sce)
)
pseudobulk <- dreamlet::aggregateToPseudoBulk(
  sce_pb,
  assay = "counts",
  fun = "sum",
  sample_id = "sample_id",
  cluster_id = "celltype",
  scale = FALSE,
  verbose = TRUE,
  checkValues = TRUE,
  BPPARAM = bpparam
)
rm(sce_pb)
S4Vectors::metadata(pseudobulk)$endome <- wc$endome_name
S4Vectors::metadata(pseudobulk)$annotation_level <- "unknown"
message("  ", length(assayNames(pseudobulk)), " cell types x ", ncol(pseudobulk), " samples")

############################################################################## #
# ---- 4. Write the outputs ----
message("Writing the processed SCE: ", smk_outputs$sce)
saveRDS(sce, smk_outputs$sce)
saveRDS(pseudobulk, smk_outputs$pseudobulk)

cells %>%
  dplyr::count(sample_id, celltype, name = "n_cells") %>%
  dplyr::mutate(endome = wc$endome_name, .before = 1) %>%
  readr::write_tsv(smk_outputs$n_cells)
