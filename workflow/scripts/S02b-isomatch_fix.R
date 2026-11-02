#!/usr/bin/env Rscript

## _________________________________________________
##
## Isomatch classify fix
##
## Aim: Make the isomatch output follow the StringTie + gffcompare conventions
## expected downstream
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S02b-isomatch_fix
## Date Created: 2026-07-15
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.3 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.3 (2026-10-28): Standardization of the script.
##    - v1.2 (2026-09-29): Merged transcripts take the id of the reference
##      transcript they contain. Added merge_transcript_id and min_sample_cnt.
##    - v1.1 (2026-09-25): Added ref_transcript_id, isom_sample_cnt and
##      isom_ref_source.
##    - v1.0 (2026-07-15): Initial version.
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
  prefix <- "Wood.control.iso_ref"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(rs, "02-Transcriptome_Assembly/isomatch_classify", paste0(prefix, ".annotated.gtf.gz")),
      present_absent = file.path(rs, "02-Transcriptome_Assembly/isomatch_merge", paste0(prefix, ".present_absent.tsv.gz")),
      track = file.path(rs, "02-Transcriptome_Assembly/isomatch_merge", paste0(prefix, ".track.tsv.gz"))
    ),
    output = list(
      gtf = file.path(rs, "temp", paste0(prefix, ".fixed.gtf"))
    ),
    params = list(
      ref_source = "ref_annotation.gtf",
      min_sample_cnt = 0
    ),
    wildcards = list(
      dataset = "Wood",
      group = "control",
      merge_method = "iso_ref"
    ),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(rtracklayer)
  library(plyranges)
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

ref_source <- params$ref_source
min_sample_cnt <- as.integer(params$min_sample_cnt)

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Set the "novel" and "NA" values of the ISOM_REF_* attributes to NA
isomRefNA <- function(x) dplyr::if_else(x %in% c("novel", "NA"), NA_character_, x)

############################################################################## #
# ---- 1. Load the isomatch outputs ----
gr <- rtracklayer::import(smk_inputs$gtf)

### One row per merged transcript, one 0/1 column per merge input file
present_absent <- readr::read_tsv(smk_inputs$present_absent)

### One row per input transcript and the merged transcript it belongs to
track <- readr::read_tsv(smk_inputs$track)

############################################################################## #
# ---- 2. Fix the isomatch output ----

#----------------------------------------------------------------------------- #
## 2.1 Long-read support and filtering ----

### The reference annotation is not counted as a sample
id_cols <- c("merged_tx_id", "merged_gene_id", "merged_exon_num", "src_tx_count_in_merged_group")
sample_cols <- setdiff(colnames(present_absent), c(id_cols, ref_source))

support <- present_absent %>%
  dplyr::transmute(
    transcript_id = merged_tx_id,
    isom_sample_cnt = rowSums(dplyr::across(dplyr::all_of(sample_cols)) > 0),
    isom_ref_source = if (nzchar(ref_source)) .data[[ref_source]] > 0 else NA
  )

kept_tx <- support %>%
  dplyr::filter(isom_ref_source %in% TRUE | isom_sample_cnt >= min_sample_cnt) %>%
  dplyr::pull(transcript_id)

message(sprintf("Kept %d of %d merged transcripts (reference or >= %d samples)",
                length(kept_tx), nrow(support), min_sample_cnt))
gr <- gr[gr$transcript_id %in% kept_tx]

#----------------------------------------------------------------------------- #
## 2.2 Reference ids ----

### If a merged transcript contains several reference transcripts, the one
### matched by classify is preferred
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
## 2.3 Metadata columns ----

### Attributes other than gene_id and transcript_id are only set on the
### transcript features, as character so that they are quoted in the GTF
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

############################################################################## #
# ---- 3. Output the GTF ----
rtracklayer::export(gr, smk_outputs$gtf, format = "gtf")
