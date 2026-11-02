#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript filtration by ORF status
##
## Aim: Keep the transcripts that match the biotypes, ORFannotate types,
## structural categories and reference isoform status of an ORF filter preset
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S04b-ORF_Filtration
## Date Created: 2025-03-27
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.1 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - "all" in any of the list filters disables that filter.
##    - A GTF without `ref_gene_type` is treated as the reference annotation:
##      gene_type and transcript_type are used as the reference biotypes and
##      the ORFannotate and structural category filters are disabled.
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.1 (2026-10-28): Standardization of the script.
##    - v1.0 (2025-03-27): Initial version.
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
  prefix <- "Wood.control.st_ref"
  orf_filter <- "ref_pc"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(rs, "04-ORF_Identification/ORF_Category", paste0(prefix, ".orf_annotated.gtf"))
    ),
    output = list(
      gtf_filter = file.path(rs, "temp", sprintf("%s.%s.orf_filter.gtf", prefix, orf_filter))
    ),
    params = list(
      main_config = orf_filter,
      valid_ref_gene_type = c("protein_coding"),
      valid_ref_tx_type = c("protein_coding"),
      valid_orfannotate_type = c("coding"),
      valid_structural_category = c("all"),
      ref_isoform_filter = "TRUE"
    ),
    wildcards = list(
      prefix = prefix,
      orf_filter = orf_filter
    ),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(GenomicRanges)
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

main_config <- params$main_config
valid_ref_gene_type <- params$valid_ref_gene_type
valid_ref_tx_type <- params$valid_ref_tx_type
valid_orfannotate_type <- params$valid_orfannotate_type
valid_structural_category <- params$valid_structural_category
ref_isoform_filter <- as.logical(params$ref_isoform_filter)

############################################################################## #
# ---- 1. Load the transcriptome ----
gtf <- rtracklayer::import(smk_inputs$gtf)

############################################################################## #
# ---- 2. Filter the transcripts ----

### ref_isoform is written as "TRUE"/"FALSE" in the GTF
valid_transcripts <- gtf %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::filter(
    (all(valid_ref_gene_type == "all") | ref_gene_type %in% valid_ref_gene_type),
    (all(valid_ref_tx_type == "all") | ref_tx_type %in% valid_ref_tx_type),
    (all(valid_orfannotate_type == "all") | orfannotate_type %in% valid_orfannotate_type),
    (all(valid_structural_category == "all") | structural_category %in% valid_structural_category),
    (is.na(ref_isoform_filter) | as.logical(ref_isoform) %in% ref_isoform_filter)
  ) %>%
  dplyr::pull(transcript_id)

if(length(valid_transcripts) == 0) stop("No transcripts pass the filters.", call. = FALSE)
message(sprintf("Preset %s: %d transcripts kept", main_config, length(valid_transcripts)))

### Gene features are not needed downstream
gtf <- gtf %>% dplyr::filter(type != "gene", transcript_id %in% valid_transcripts)

############################################################################## #
# ---- 3. Output the GTF ----
rtracklayer::export(gtf, smk_outputs$gtf_filter)
