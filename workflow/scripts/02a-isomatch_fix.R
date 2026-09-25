#!/usr/bin/env Rscript

## _________________________________________________
##
## Isomatch Classify Fix
##
## Author: Guillermo Rocamora Pérez
## Date Created: 2026-07-15
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.1 (2026-09-25)
## _________________________________________________
##
## - Notes:
##   + transcript_id: isomatch's own id (ISOMT_*) is kept. Several merged
##   transcripts can match the same reference transcript (e.g. the reference
##   isoform and its alternative-end versions), so the reference id cannot be
##   used as transcript_id.
##   + ISOM_REF_TX_ID -> ref_transcript_id (NA for novel transcripts)
##   + ISOM_REF_GENE_ID -> ref_gene_id (NA when absent)
##   + ISOM_REF_GENE_NAME -> gene_name (NA when absent)
##   + isom_sample_cnt: number of long-read samples the merged transcript was
##   found in (transcript features only).
##   + isom_ref_source: TRUE when the reference annotation contributed to the
##   merged transcript. Only set when the merge included the reference.
##   + gene_id and every other attribute (ISOM_CATEGORY, ISOM_SUBCATEGORY, ...)
##   are left untouched.
##
## - Changelog:
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
  snakemake <- Snakemake(
      input = list(gtf = "data/test_data/test.isomatch_classify.annotated.gtf.gz",
                   present_absent = "data/test_data/test.isomatch_merge.present_absent.tsv.gz"),
      output = list(gtf = "data/test_data/test.isomatch_fix.annotated.gtf"),
      params = list(ref_source = "ref_annotation.gtf"),
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

## Import the isomatch classify output and the merge present/absent table (one
## row per merged transcript, one 0/1 column per merge input file)
gr <- rtracklayer::import(snakemake@input$gtf)
present_absent <- readr::read_tsv(snakemake@input$present_absent)
ref_source <- snakemake@params$ref_source

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

## Fix the metadata columns to match the output expected from a gffcompare
## pipeline. transcript_id is kept unique (isomatch id). Novel transcripts carry
## "novel" or "NA" in the ISOM_REF_* attributes, both converted to NA. The
## support columns are stored as character so the GTF values stay quoted.
isomRefNA <- function(x) dplyr::if_else(x %in% c("novel", "NA"), NA_character_, x)

metadata <- mcols(gr) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    ref_transcript_id = isomRefNA(ISOM_REF_TX_ID),
    ref_gene_id = isomRefNA(ISOM_REF_GENE_ID),
    gene_name = isomRefNA(ISOM_REF_GENE_NAME)
  ) %>%
  dplyr::select(-ISOM_REF_TX_ID, -ISOM_REF_GENE_ID, -ISOM_REF_GENE_NAME) %>%
  dplyr::left_join(support, by = "transcript_id") %>%
  dplyr::mutate(dplyr::across(c(isom_sample_cnt, isom_ref_source), ~ dplyr::if_else(type == "transcript", as.character(.x), NA_character_))) %>%
  S4Vectors::DataFrame()
mcols(gr) <- metadata

## Export the modified GRanges
export(gr, snakemake@output$gtf, format = "gtf")