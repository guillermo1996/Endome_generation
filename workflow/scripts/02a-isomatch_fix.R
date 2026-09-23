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
## Latest version: v1.0 (2026-07-15)
## _________________________________________________
##
## - Notes:
##   + transcript_id: replaced with ISOM_REF_TX_ID when it is a real reference
##   match; otherwise isomatch's own id is kept.
##   + ISOM_REF_GENE_ID -> ref_gene_id
##   + ISOM_REF_GENE_NAME -> gene_name
##   + gene_id and every other attribute (ISOM_CATEGORY, ISOM_SUBCATEGORY, ...) 
##   are left untouched.
## 
## - Changelog:
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
      input = list(gtf = "data/test_data/test.isomatch_classify.annotated.gtf.gz"),
      output = list(gtf = "data/test_data/test.isomatch_fix.annotated.gtf"),
      params = list(),
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

## Import the isomatch classify output
gr <- rtracklayer::import(snakemake@input$gtf)

## Fix the metadata columns to match the output expected from a gffcompare pipeline
metadata <- mcols(gr) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(transcript_id = dplyr::case_when(
    ISOM_REF_TX_ID != "novel" &&  ~ ISOM_REF_TX_ID,
    .default = transcript_id
  )) %>% 
  dplyr::rename(ref_gene_id = ISOM_REF_GENE_ID, gene_name = ISOM_REF_GENE_NAME) %>% 
  S4Vectors::DataFrame()
mcols(gr) <- metadata

## Export the modified GRanges
export(gr, snakemake@output$gtf, format = "gtf")