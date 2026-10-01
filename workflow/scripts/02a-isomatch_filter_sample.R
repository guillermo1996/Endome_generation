#!/usr/bin/env Rscript

## _________________________________________________
##
## Isomatch Sample Filter
##
## Removes low-abundance transcripts from a StringTie sample assembly
## before it is indexed and merged by isomatch. It reproduces the input
## filter of `stringtie --merge` (-F min_fpkm, -T min_tpm; both 1.0 by
## default in StringTie 3.0.0)
##
## Author: Guillermo Rocamora Pérez
## Date Created: 2026-09-29
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2026-09-29)
## _________________________________________________
##
## - Notes:
##
## - Changelog:
##   + v1.0 (2026-09-29): first version.
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
  ## dataset.group.merge_method folder created by the `test_data` rule. Sample
  ## assemblies are not linked by it: copy one there as test.stringtie_sample.gtf
  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  snakemake <- Snakemake(
      input = list(gtf = file.path(test_dir, "test.stringtie_sample.gtf")),
      output = list(gtf = file.path(test_dir, "interactive.isomatch_filter_sample.gtf")),
      params = list(min_fpkm = 1, min_tpm = 1),
      wildcards = list(sample = "test"),
      threads = 1
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
suppressMessages({
  library(rtracklayer)
  library(tidyverse)
})

## Package options - these are some examples and do not modify the script usage.
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")

############################################################################## #
# ---- 1. Filter the sample assembly ----

## Import the StringTie sample assembly. FPKM and TPM are only set on the
## transcript features.
gr <- rtracklayer::import(snakemake@input$gtf)
sample_name <- as.character(snakemake@wildcards$sample)
min_fpkm <- as.numeric(snakemake@params$min_fpkm)
min_tpm <- as.numeric(snakemake@params$min_tpm)

## Transcripts passing both thresholds
kept_tx <- mcols(gr) %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(FPKM = as.numeric(FPKM), TPM = as.numeric(TPM)) %>%
  dplyr::filter(FPKM >= min_fpkm, TPM >= min_tpm) %>%
  dplyr::pull(transcript_id)

message(sprintf("Kept %d of %d transcripts in sample %s (FPKM >= %s and TPM >= %s)",
                length(kept_tx), sum(gr$type == "transcript"), sample_name, min_fpkm, min_tpm))

## Keep the transcripts and their exons, and export the filtered assembly
gr <- gr[gr$transcript_id %in% kept_tx]
export(gr, snakemake@output$gtf, format = "gtf")
