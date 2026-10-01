#!/usr/bin/env Rscript

## _________________________________________________
##
## Isomatch Classify Fix
##
## Make the isomatch output follow the StringTie + gffcompare
## conventions expected downstream (e.g. 04a joins GENCODE by transcript_id).
##
## Author: Guillermo Rocamora Pérez
## Date Created: 2026-07-15
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.2 (2026-09-29)
## _________________________________________________
##
## - Notes:
##
## - Changelog:
##   + v1.2 (2026-09-29): transcripts containing a reference transcript take
##   its id; added merge_transcript_id and the min_sample_cnt filter; exon
##   features only keep gene_id and transcript_id.
##   + v1.1 (2026-09-25): fixed syntax error; transcript_id is no longer
##   replaced by the reference id; added ref_transcript_id, isom_sample_cnt and
##   isom_ref_source.
##
## - Please contact guillermorocamora@gmail.com for further assistance.
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
      threads = 'numeric'
    )
  )
  ## dataset.group.merge_method folder created by the `test_data` rule. Outputs
  ## use the "interactive" prefix so they never overwrite the linked results.
  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  snakemake <- Snakemake(
      input = list(gtf = file.path(test_dir, "test.annotated.gtf.gz"),
                   present_absent = file.path(test_dir, "test.present_absent.tsv.gz"),
                   track = file.path(test_dir, "test.track.tsv.gz")),
      output = list(gtf = file.path(test_dir, "interactive.fixed.gtf")),
      params = list(ref_source = "ref_annotation.gtf", min_sample_cnt = 0),
      wildcards = list(dataset = "Ebbert", group = "control"),
      threads = 1
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
suppressMessages({
  library(rtracklayer)
  library(plyranges)
  library(tidyverse)
})

## Package options - these are some examples and do not modify the script usage.
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)

############################################################################## #
# ---- 1. Fix isomatch output ----

## Import the isomatch classify output, the merge present/absent table (one
## row per merged transcript, one 0/1 column per merge input file) and the
## merge track table (one row per input transcript and its merged transcript)
gr <- rtracklayer::import(snakemake@input$gtf)
present_absent <- readr::read_tsv(snakemake@input$present_absent)
track <- readr::read_tsv(snakemake@input$track)
ref_source <- snakemake@params$ref_source
min_sample_cnt <- as.integer(snakemake@params$min_sample_cnt)

## Novel transcripts carry "novel" or "NA" in the ISOM_REF_* attributes, both
## converted to NA
isomRefNA <- function(x) dplyr::if_else(x %in% c("novel", "NA"), NA_character_, x)

#----------------------------------------------------------------------------- #
## 1.1 Long-read support and filtering ----

## Long-read support per merged transcript. The reference annotation column is
## excluded from the sample count and reported on its own.
id_cols <- c("merged_tx_id", "merged_gene_id", "merged_exon_num", "src_tx_count_in_merged_group")
sample_cols <- setdiff(colnames(present_absent), c(id_cols, ref_source))

support <- present_absent %>%
  dplyr::transmute(
    transcript_id = merged_tx_id,
    isom_sample_cnt = rowSums(dplyr::across(dplyr::all_of(sample_cols)) > 0),
    isom_ref_source = if (nzchar(ref_source)) .data[[ref_source]] > 0 else NA
  )

## Keep the reference transcripts regardless of their long-read support, and
## the rest only when found in at least `min_sample_cnt` samples
kept_tx <- support %>%
  dplyr::filter(isom_ref_source %in% TRUE | isom_sample_cnt >= min_sample_cnt) %>%
  dplyr::pull(transcript_id)

message(sprintf("Kept %d of %d merged transcripts (reference or >= %d samples)",
                length(kept_tx), nrow(support), min_sample_cnt))
gr <- gr[gr$transcript_id %in% kept_tx]

#----------------------------------------------------------------------------- #
## 1.2 Reference ids ----

## Merged transcripts containing a reference transcript take its id. When they
## contain several, the one matched by classify goes first.
tx_ref_match <- mcols(gr) %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::transmute(transcript_id, ref_transcript_id = isomRefNA(ISOM_REF_TX_ID))

new_ids <- track %>%
  dplyr::filter(src_file_name == ref_source) %>%
  dplyr::select(transcript_id = merged_tx_id, new_transcript_id = src_tx_id) %>%
  dplyr::inner_join(tx_ref_match, by = "transcript_id") %>%
  dplyr::arrange(transcript_id, dplyr::desc(new_transcript_id %in% ref_transcript_id)) %>%
  dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
  dplyr::select(transcript_id, new_transcript_id)

message(sprintf("Renamed %d merged transcripts to their reference id", nrow(new_ids)))

#----------------------------------------------------------------------------- #
## 1.3 Fix the metadata columns ----

## Fix the metadata columns to match the output expected from a gffcompare
## pipeline. Every attribute other than gene_id and transcript_id is only set
## on the transcript features, stored as character so the GTF values stay
## quoted.
tx_cols <- c("ISOM_CATEGORY", "ISOM_SUBCATEGORY", "ref_transcript_id", "ref_gene_id", "gene_name",
             "isom_sample_cnt", "isom_ref_source")

metadata <- mcols(gr) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    ref_transcript_id = isomRefNA(ISOM_REF_TX_ID),
    ref_gene_id = isomRefNA(ISOM_REF_GENE_ID),
    gene_name = isomRefNA(ISOM_REF_GENE_NAME)
  ) %>%
  dplyr::select(-ISOM_REF_TX_ID, -ISOM_REF_GENE_ID, -ISOM_REF_GENE_NAME) %>%
  dplyr::left_join(support, by = "transcript_id") %>%
  dplyr::left_join(new_ids, by = "transcript_id") %>%
  dplyr::mutate(
    merge_transcript_id = transcript_id,
    transcript_id = dplyr::coalesce(new_transcript_id, transcript_id)
  ) %>%
  dplyr::select(-new_transcript_id) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(tx_cols), ~ dplyr::if_else(type == "transcript", as.character(.x), NA_character_))) %>%
  dplyr::relocate(merge_transcript_id, .after = dplyr::last_col()) %>%
  S4Vectors::DataFrame()
mcols(gr) <- metadata

## Export the modified GRanges
export(gr, snakemake@output$gtf, format = "gtf")