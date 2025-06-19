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
      params = 'list'
    )
  )
  snakemake <- Snakemake(
    input = list(gtf = "~/RytenLab-Research/38-Endome_generation/results_k15/04-ORF_Filtration/ORF_Category/sq3.annotated_orf.gtf"),
    output = list(gtf = "~/RytenLab-Research/38-Endome_generation/results_k15/04-ORF_Filtration/ORF_Filter/sq3.annotated_orf.filter.gtf"),
    params = list(
      valid_ref_gene_type = c("protein_coding"),
      valid_ref_tx_type = c("protein_coding"),
      valid_sq3_type = c("known-ORF", "novel-ORF"),
      hasCDS_filter = TRUE,
      in_ref_filter = TRUE
    )
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
# ---- 2. Filter the GTF based on configuration ----

#----------------------------------------------------------------------------- #
## 2.1 Load filters ----
valid_ref_gene_type <- snakemake@params$valid_ref_gene_type
valid_ref_tx_type <- snakemake@params$valid_ref_tx_type
valid_sq3_type <- snakemake@params$valid_sq3_type
hasCDS_filter <- snakemake@params$hasCDS_filter
in_ref_filter <- snakemake@params$in_ref_filter


#----------------------------------------------------------------------------- #
## 2.2 Apply the filterst o extract the valid transcripts ----

### GR: Note that for each of the `_type` filters, the value "all" can be
### provided to ignore any filtering (which will leave NAs). The filter `in_ref`
### removes every entry with NA in the `ref_gene_type` and `ref_tx_type` columns
valid_transcripts <- gtf[gtf$type == "transcript"] %>% plyranges::filter(
  (all(valid_ref_gene_type == "all") | ref_gene_type %in% valid_ref_gene_type),
  (all(valid_ref_tx_type == "all") | ref_tx_type %in% valid_ref_tx_type),
  (all(valid_sq3_type == "all") | sq3_type %in% valid_sq3_type),
  hasCDS == hasCDS_filter,
  in_ref = in_ref_filter,
)

############################################################################## #
# ---- 3. Filter the transcriptome ----

## GR: Genes are removed from the transcriptome as they will not be needed in
## future modules.
gtf <- gtf %>% plyranges::filter(type != "gene" & transcript_id %in% valid_transcripts$transcript_id)

############################################################################## #
# ---- 4. Output the GTF ----
gtf %>% rtracklayer::export(output_gtf_path)
