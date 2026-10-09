## _________________________________________________
##
## Processed scUTRquant SCE and pseudobulk
##
## Aim: Build the processed SingleCellExperiment of one ENDome (QC-passing,
##      annotated nuclei from 06c; bins annotated with their ENDome coordinates)
##      and sum its raw counts by sample (donor x region) and cell type.
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
##    - Processed SCE (equivalent of the annotated SCE of Aine's 04c/04d and
##      sections 1-2 of 41/Scripts/02-satuRn_DTU/T01-Hardy_satuRn.R):
##        - Nuclei: the QC-passing nuclei of 06c, with the full nucleus table
##          (sample, library, annotation levels, total_umis, effective_txs) as
##          colData.
##        - Bins: bins with no counts in the retained nuclei are dropped (T01
##          1.2 drops them before the nucleus QC; no result depends on it).
##        - rowData: transcript_id and gene_id from scUTRquant, plus, when the
##          inputs are given:
##            - gtf (the ENDome GTF scUTRquant was built from,
##              05-Truncation/txendcutr/{endome}.txendcutr.gtf): bin coordinates
##              and the rank of the ENDome-defining end within its gene. This
##              generalises Aine's utr_rank / is_shortest_utr / is_longest_utr to
##              5' ENDomes: the end is the 3' end for 3p ENDomes and the 5' end
##              for 5p ENDomes, ranked in the direction of transcription
##              (1 = most upstream; for 3p ENDomes, the shortest 3' UTR).
##              Column names avoid seqnames/start/end/strand, which clash with
##              SummarizedExperiment (Aine removed them later for that reason).
##              Bins on non-standard chromosomes are kept, not filtered.
##            - merge_tsv (05-Truncation/txendcutr/{endome}.txendcutr.merge.tsv):
##              the ENDome transcripts collapsed into each bin.
##    - Pseudobulk: same function and grouping as Aine's pipeline and the
##      manuscript (dreamlet::aggregateToPseudoBulk, nuclei summed by sample_id,
##      so replicate libraries B1/B2 are summed together, and cell type). The
##      bin annotations are carried to the pseudobulk rowData.
##    - Difference with Aine's pipeline: the RAW counts are summed, not the
##      scran-normalised counts. satuRn (06e) models each bin's count out of
##      its gene's total count within a sample, so it needs counts, and
##      within-sample proportions do not need library-size normalisation (see
##      docs/satuRn/README.md). No size factors / normcounts are computed.
##    - No minimum number of nuclei is applied here: the per-sample, per-cell
##      type nucleus counts are written next to the pseudobulk and 06e applies
##      min_cells.
##
## Changelog:
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
  test_orf_filter <- "ref_pc"
  test_width <- "500"
  test_txEnd <- "3p"
  test_dataset_name <- "hardy_snRNAseq"
  test_endome <- paste0(test_orf_filter, ".w", test_width, ".", test_txEnd)
  test_sce <- paste0(test_endome, ".", test_dataset_name)

  snakemake <- Snakemake(
    input = list(
      sce = file.path(test_dir, paste0("test.", test_sce, ".txs.Rds")),
      cells = file.path(test_dir, paste0("test.", test_endome, ".cells.tsv")),
      gtf = file.path(test_dir, paste0("test.", test_endome, ".txendcutr.gtf")),
      merge_tsv = file.path(test_dir, paste0("test.", test_endome, ".txendcutr.merge.tsv"))
    ),
    output = list(
      sce = file.path(test_dir, paste0("interactive.", test_sce, ".processed.txs.rds")),
      pseudobulk = file.path(test_dir, paste0("interactive.", test_sce, ".pseudobulk.rds")),
      n_cells = file.path(test_dir, paste0("interactive.", test_sce, ".pseudobulk_n_cells.tsv"))
    ),
    params = list(
      cell_type_level = "annotation_level_2"
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

cell_type_level <- params$cell_type_level
bpparam <- BiocParallel::MulticoreParam(workers = snakemake@threads)

## Truncation site of the ENDome ({orf_filter}.w{width}.{txEnd}): which end defines the bins
tx_end <- stringr::str_extract(wc$endome_name, "(3p|5p)$")

dir.create(dirname(smk_outputs$pseudobulk), showWarnings = FALSE, recursive = TRUE)

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Bin coordinates from the ENDome GTF
#'
#' Uses the transcript records of the truncated ENDome GTF (one per bin
#' representative, i.e. the SCE row names). end_coord is the genomic position
#' of the ENDome-defining end: the 3' end for 3p ENDomes, the 5' end for 5p
#' ENDomes. Following Aine's 04c (signed_end), it is signed by strand so that
#' ranking it orders the bins in the direction of transcription.
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

#----------------------------------------------------------------------------- #
## 1.1 QC-passing nuclei (06c) ----
cells <- readr::read_tsv(smk_inputs$cells) %>%
  dplyr::filter(pass_qc)
if (!cell_type_level %in% colnames(cells)) {
  stop("Annotation level '", cell_type_level, "' not found in the nucleus table. Columns: ",
       paste(colnames(cells), collapse = ", "))
}
message("QC-passing nuclei: ", nrow(cells), " in ", dplyr::n_distinct(cells$sample_id), " samples")

#----------------------------------------------------------------------------- #
## 1.2 scUTRquant SCE ----
message("Loading the scUTRquant SCE: ", smk_inputs$sce)
sce <- readRDS(smk_inputs$sce)

missing_cells <- setdiff(cells$cell_id, colnames(sce))
if (length(missing_cells) > 0) {
  stop(length(missing_cells), " nuclei of the nucleus table are not in the SCE. ",
       "Was 06c run on this ENDome? First missing: ", missing_cells[1])
}

############################################################################## #
# ---- 2. Processed SCE ----

#----------------------------------------------------------------------------- #
## 2.1 Nuclei: QC-passing nuclei with the 06c nucleus table ----
sce <- sce[, cells$cell_id]
colData(sce) <- S4Vectors::DataFrame(
  cells %>% dplyr::select(-dplyr::any_of(c("endome", "pass_qc"))) %>% as.data.frame(),
  row.names = cells$cell_id
)

#----------------------------------------------------------------------------- #
## 2.2 Bins: drop bins without counts ----
keep_bins <- Matrix::rowSums(counts(sce)) > 0
message("Bins with counts in the retained nuclei: ", sum(keep_bins), " of ", length(keep_bins))
sce <- sce[keep_bins, ]

#----------------------------------------------------------------------------- #
## 2.3 Bin annotations ----
bin_annots <- tibble::tibble(
  transcript_id = rowData(sce)$transcript_id,
  gene_id = rowData(sce)$gene_id
)

if (!is.null(smk_inputs$gtf)) {
  message("Adding bin coordinates from the ENDome GTF (", tx_end, " ends): ", smk_inputs$gtf)
  bin_annots <- bin_annots %>%
    dplyr::left_join(read_bin_coordinates(smk_inputs$gtf, tx_end), by = "transcript_id") %>%
    ## Rank of the ENDome-defining end within the gene, in the direction of
    ## transcription (Aine's utr_rank / is_shortest_utr / is_longest_utr)
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
    dplyr::mutate(is_mt = bin_chrom == "chrM")
  n_no_coords <- sum(is.na(bin_annots$end_coord))
  if (n_no_coords > 0) warning(n_no_coords, " bins are not in the ENDome GTF (NA coordinates)")
}

if (!is.null(smk_inputs$merge_tsv)) {
  message("Adding the merged ENDome transcripts of each bin: ", smk_inputs$merge_tsv)
  bin_annots <- bin_annots %>%
    dplyr::left_join(read_bin_members(smk_inputs$merge_tsv), by = "transcript_id")
}

rowData(sce) <- S4Vectors::DataFrame(as.data.frame(bin_annots), row.names = rownames(sce))
S4Vectors::metadata(sce)$endome <- wc$endome_name
S4Vectors::metadata(sce)$tx_end <- tx_end

############################################################################## #
# ---- 3. Pseudobulk ----
message("Pseudobulking raw counts by sample_id x ", cell_type_level, "...")
## Only sample_id and the cell type go to aggregateToPseudoBulk, so no
## per-nucleus column (library, annotation levels, QC metrics) ends up in the
## pseudobulk colData as if it were a sample-level value. The counts are shared,
## not copied.
sce_pb <- sce
colData(sce_pb) <- S4Vectors::DataFrame(
  sample_id = sce$sample_id,
  celltype = sce[[cell_type_level]],
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
S4Vectors::metadata(pseudobulk)$annotation_level <- cell_type_level
message("  ", length(assayNames(pseudobulk)), " cell types x ", ncol(pseudobulk), " samples")

############################################################################## #
# ---- 4. Outputs ----
message("Writing the processed SCE: ", smk_outputs$sce)
saveRDS(sce, smk_outputs$sce)
saveRDS(pseudobulk, smk_outputs$pseudobulk)

## Nuclei per sample and cell type (06e drops pseudobulks below min_cells)
cells %>%
  dplyr::count(sample_id, celltype = .data[[cell_type_level]], name = "n_cells") %>%
  dplyr::mutate(endome = wc$endome_name, annotation_level = cell_type_level, .before = 1) %>%
  readr::write_tsv(smk_outputs$n_cells)
