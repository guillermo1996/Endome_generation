## _________________________________________________
##
## Processed scUTRquant SCE and pseudobulk
##
## Aim: Build the processed SingleCellExperiment of one ENDome (QC-passing
##      nuclei from 07a, bins annotated with their ENDome coordinates) and sum
##      its raw counts by sample and cell type.
##
## Project: ENDome generation - UTR Usage
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-06
##
## Latest Version: v2.0 (2026-10-16)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - Input formats, outputs and the reasoning behind each step are described
##      in docs/UTR_quantification_and_usage.md.
##    - Raw counts are summed (not normalised counts): satuRn (07c) models each
##      bin's share of its gene's counts within a sample.
##    - No minimum number of nuclei is applied here: 07c applies min_cells using
##      the n_cells table written next to the pseudobulk.
##
## Changelog:
##    - v2.0 (2026-10-16): Simplified. The GTF and merge table are required.
##    - v1.1 (2026-10-07): Writes the processed SCE; bin annotations from the
##      ENDome GTF and merge table; drops bins without counts.
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
      cells = file.path(test_dir, paste0("test.", test_endome, ".cells.tsv")),
      gtf = file.path(test_dir, paste0("test.", test_endome, ".txendcutr.gtf")),
      merge_tsv = file.path(test_dir, paste0("test.", test_endome, ".txendcutr.merge.tsv"))
    ),
    output = list(
      sce = file.path(test_dir, paste0("interactive.", test_endome, ".processed.txs.rds")),
      pseudobulk = file.path(test_dir, paste0("interactive.", test_endome, ".pseudobulk.rds")),
      n_cells = file.path(test_dir, paste0("interactive.", test_endome, ".pseudobulk_n_cells.tsv"))
    ),
    wildcards = list(
      endome_name = paste0(basename(test_dir), ".", test_endome)
    ),
    log = list(),
    threads = 4,
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
  library(dreamlet)
  library(BiocParallel)
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

bpparam <- BiocParallel::MulticoreParam(workers = snakemake@threads)

## Truncation site of the ENDome (the name ends in .w{width}.{txEnd})
tx_end <- stringr::str_extract(wc$endome_name, "(3p|5p)$")

dir.create(dirname(smk_outputs$pseudobulk), showWarnings = FALSE, recursive = TRUE)

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Bin coordinates from the ENDome GTF
#'
#' end_coord is the position of the end that defines the bins: the 3' end for
#' 3p ENDomes, the 5' end for 5p ENDomes.
#'
#' @param path Path to the ENDome GTF
#' @param tx_end "3p" or "5p"
#' @return A tibble with transcript_id, bin_chrom, bin_start, bin_end, bin_strand, end_coord
read_bin_coordinates <- function(path, tx_end) {
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

#' ENDome transcripts collapsed into each bin (txendcutr merge table)
#'
#' @param path Path to {endome}.txendcutr.merge.tsv (tx_in, tx_out, gene_out)
#' @return A tibble with transcript_id (= tx_out), n_merged_txs and merged_txs
read_bin_members <- function(path) {
  readr::read_tsv(path) %>%
    dplyr::group_by(transcript_id = tx_out) %>%
    dplyr::summarise(n_merged_txs = dplyr::n(), merged_txs = paste(sort(tx_in), collapse = ";"), .groups = "drop")
}

############################################################################## #
# ---- 1. Load inputs ----
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
## 2.1 Nuclei: the QC-passing nuclei, with the 07a table as colData ----
sce <- sce[, cells$cell_id]
colData(sce) <- S4Vectors::DataFrame(
  cells %>% dplyr::select(-endome, -pass_qc) %>% as.data.frame(),
  row.names = cells$cell_id
)

#----------------------------------------------------------------------------- #
## 2.2 Bins: drop bins without counts ----
keep_bins <- Matrix::rowSums(counts(sce)) > 0
message("Bins with counts in the retained nuclei: ", sum(keep_bins), " of ", length(keep_bins))
sce <- sce[keep_bins, ]

#----------------------------------------------------------------------------- #
## 2.3 Bin annotations ----
## end_rank orders the bins of a gene by the position of their end, in the
## direction of transcription (1 = most upstream; for 3p ENDomes, the shortest
## 3' UTR). Coordinates are named bin_* because seqnames/start/end/strand clash
## with SummarizedExperiment.
message("Adding bin annotations from the ENDome GTF and merge table")
bin_annots <- tibble::tibble(
  transcript_id = rowData(sce)$transcript_id,
  gene_id = rowData(sce)$gene_id
) %>%
  dplyr::left_join(read_bin_coordinates(smk_inputs$gtf, tx_end), by = "transcript_id") %>%
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
  dplyr::mutate(is_mt = bin_chrom == "chrM") %>%
  dplyr::left_join(read_bin_members(smk_inputs$merge_tsv), by = "transcript_id")

rowData(sce) <- S4Vectors::DataFrame(as.data.frame(bin_annots), row.names = rownames(sce))
S4Vectors::metadata(sce)$endome <- wc$endome_name
S4Vectors::metadata(sce)$tx_end <- tx_end

############################################################################## #
# ---- 3. Pseudobulk ----
## Only sample_id and the cell type are passed to aggregateToPseudoBulk, so no
## per-nucleus column ends up in the pseudobulk colData. Result: one assay per
## cell type, one column per sample.
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
# ---- 4. Outputs ----
message("Writing the processed SCE: ", smk_outputs$sce)
saveRDS(sce, smk_outputs$sce)
saveRDS(pseudobulk, smk_outputs$pseudobulk)

## Nuclei per sample and cell type (07c drops pseudobulks below min_cells)
cells %>%
  dplyr::count(sample_id, celltype, name = "n_cells") %>%
  dplyr::mutate(endome = wc$endome_name, .before = 1) %>%
  readr::write_tsv(smk_outputs$n_cells)
