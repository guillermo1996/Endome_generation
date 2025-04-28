#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript filtration by ORF status
##
## Aim: Filter the transcripts prior to truncation based on the ORF status
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-03-27
##
## Copyright (c) Guillermo Rocamora Pérez, 2025
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-03-27)
##
## _________________________________________________
##
## - Notes:
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
    input = list(gtf = "/home/grocamora/RytenLab-Research/38-Endome_generation/results/ORF_Category/sq3.annotated_orf.gtf"),
    output = list(gtf = "/home/grocamora/RytenLab-Research/38-Endome_generation/results/ORF_Filter/sq3.annotated_orf.filter.gtf"),
    params = list(keep_known_orf = FALSE,
                  keep_novel_orf = FALSE,
                  keep_not_orf = FALSE,
                  keep_hasCDS = TRUE),
    threads = 1
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

## Package options - these are some examples and do not modify the script usage.
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)

#----------------------------------------------------------------------------- #
## 0.2 Script Paths ----
main_path <- "/home/grocamora/RytenLab-Research/38-Endome_generation" # here::here()

### Input Paths
input_gtf_path <- snakemake@input$gtf

### Output Paths
output_gtf_path <- snakemake@output$gtf

#### Create output directory
dir.create(dirname(output_gtf_path), showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.3 Script Parameters ----
#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----
#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----

############################################################################## #
# ---- 1. Load transcriptome ----
gtf <- rtracklayer::import.gff(input_gtf_path)


############################################################################## #
# ---- 2. Split transcripts by ORF categories ----
transcripts <- gtf[gtf$type == "transcript"]

tx_known <- transcripts$transcript_id[transcripts$ORF_cat == "known-ORF"]
tx_novel <- transcripts$transcript_id[transcripts$ORF_cat == "novel-ORF"]
tx_not <- transcripts$transcript_id[transcripts$ORF_cat == "not-ORF"]

tx_hasCDS <- transcripts$transcript_id[transcripts$hasCDS == "true"]

## GR: Define the transcripts to keep based on the config file of the pipeline
tx_keep <- c()

if(snakemake@params$keep_known_orf) tx_keep <- c(tx_keep, tx_known)
if(snakemake@params$keep_novel_orf) tx_keep <- c(tx_keep, tx_novel)
if(snakemake@params$keep_not_orf) tx_keep <- c(tx_keep, tx_not)
if(snakemake@params$keep_hasCDS) tx_keep <- c(tx_keep, tx_hasCDS)

tx_keep <- unique(tx_keep)

############################################################################## #
# ---- 3. Filter the transcriptome ----

## GR: Genes are removed from the transcriptome as they will not be needed in
## future modules.
gtf <- gtf %>% plyranges::filter(type != "gene" & transcript_id %in% tx_keep)

############################################################################## #
# ---- 4. Output the GTF ----
gtf %>% rtracklayer::export(output_gtf_path)
