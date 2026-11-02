#!/usr/bin/env Rscript

## _________________________________________________
##
## Isomatch sample filter
##
## Aim: Remove the low-abundance transcripts of a StringTie sample assembly
## before it is indexed and merged by isomatch
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S02a-isomatch_filter_sample
## Date Created: 2026-09-29
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.1 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - Reproduces the input filter of `stringtie --merge` (-F min_fpkm and
##      -T min_tpm, both 1.0 by default in StringTie 3.0.0).
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.1 (2026-10-28): Standardization of the script.
##    - v1.0 (2026-09-29): Initial version.
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
  merge_method <- "iso_ref"
  sample <- "LIB01_PDC029_F_CTX"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(rs, "02-Transcriptome_Assembly/StringTie_Sample_Assembly", paste0(sample, ".gtf"))
    ),
    output = list(
      gtf = file.path(rs, "temp", sprintf("%s.%s.filter.gtf", merge_method, sample))
    ),
    params = list(
      min_fpkm = 1,
      min_tpm = 1
    ),
    wildcards = list(
      merge_method = merge_method,
      sample = sample
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

min_fpkm <- as.numeric(params$min_fpkm)
min_tpm <- as.numeric(params$min_tpm)

############################################################################## #
# ---- 1. Filter the sample assembly ----
gr <- rtracklayer::import(smk_inputs$gtf)

### FPKM and TPM are only set on the transcript features
kept_tx <- mcols(gr) %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(FPKM = as.numeric(FPKM), TPM = as.numeric(TPM)) %>%
  dplyr::filter(FPKM >= min_fpkm, TPM >= min_tpm) %>%
  dplyr::pull(transcript_id)

message(sprintf("Kept %d of %d transcripts in sample %s (FPKM >= %s and TPM >= %s)",
                length(kept_tx), sum(gr$type == "transcript"), wc$sample, min_fpkm, min_tpm))

############################################################################## #
# ---- 2. Output the GTF ----
gr <- gr[gr$transcript_id %in% kept_tx]
rtracklayer::export(gr, smk_outputs$gtf, format = "gtf")
