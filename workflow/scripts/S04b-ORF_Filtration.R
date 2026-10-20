#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript filtration by ORF status
##
## Author: Guillermo Rocamora Pérez
## Date Created: 2025-03-27
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-03-27)
## _________________________________________________
##
## - Notes:
##   + Keeps only transcripts matching, per `params`: valid_ref_gene_type,
##   valid_ref_tx_type (reference gene/transcript biotypes),
##   valid_orfannotate_type (ORFannotate coding_class) and ref_isoform_filter.
##   + If the input GTF has no "ref_gene_type" column (i.e. it's the raw
##   reference annotation, not 04a's output), gene_type/transcript_type are
##   copied into ref_gene_type/ref_tx_type and the ORFannotate-specific
##   filters are disabled, so the same filtering logic applies to both.
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

  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  test_orf_filter <- "ref_pc"

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(test_dir, "test.orf_annotated.gtf")
    ),
    output = list(
      gtf_filter = file.path(test_dir, paste0("interactive.", test_orf_filter, ".orf_filter.gtf"))
    ),
    params = list(
      main_config = test_orf_filter,
      valid_ref_gene_type = c("protein_coding"),
      valid_ref_tx_type = c("protein_coding"),
      valid_orfannotate_type = c("coding"),
      valid_structural_category = c("all"),
      ref_isoform_filter = "TRUE"
    )
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(plyranges)
  library(tidyverse)
})

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
output_gtf_path <- snakemake@output$gtf_filter

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
gtf <- rtracklayer::import(input_gtf_path)

############################################################################## #
# ---- 2. Filter the GTF based on configuration ----

#----------------------------------------------------------------------------- #
## 2.1 Load filters ----
main_config <- snakemake@params$main_config
valid_ref_gene_type <- snakemake@params$valid_ref_gene_type
valid_ref_tx_type <- snakemake@params$valid_ref_tx_type
valid_orfannotate_type <- snakemake@params$valid_orfannotate_type
valid_structural_category <- snakemake@params$valid_structural_category
ref_isoform_filter <- as.logical(snakemake@params$ref_isoform_filter)

#----------------------------------------------------------------------------- #
## 2.2 Apply the filters and extract the valid transcripts ----

### GR: If "ref_gene_type" is not found, we are dealing with the reference
### annotation. Modify their columns so that the same logic can be applied on
### later steps
if(!"ref_gene_type" %in% colnames(mcols(gtf))){
  if("coding" %in% valid_orfannotate_type){
    valid_ref_gene_type = c("protein_coding")
    valid_ref_tx_type = c("protein_coding")
  }
  
  mcols(gtf)["ref_gene_type"] <- mcols(gtf)["gene_type"]
  mcols(gtf)["ref_tx_type"] <- mcols(gtf)["transcript_type"]
  mcols(gtf)["ref_isoform"] <- "TRUE"
  valid_orfannotate_type <- "all"
  mcols(gtf)["orfannotate_type"] <- ""
  valid_structural_category <- "all"
  mcols(gtf)["structural_category"] <- ""
}

### GR: Note that for each of the list filters, the value "all" can be
### provided to ignore any filtering, and ref_isoform_filter = "" (NA) ignores
### the ref_isoform filter. ref_isoform is written as "TRUE"/"FALSE" in the GTF.
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

if(identical(valid_transcripts, character(0))) stop("No transcripts pass the requirements.")
message(sprintf("Preset %s: %d transcripts kept", main_config, length(valid_transcripts)))

############################################################################## #
# ---- 3. Filter the transcriptome ----

## GR: Genes are removed from the transcriptome as they will not be needed in
## future modules.
gtf <- gtf %>% dplyr::filter(type != "gene", transcript_id %in% valid_transcripts)

############################################################################## #
# ---- 4. Output the GTF ----
gtf %>% rtracklayer::export(output_gtf_path)
