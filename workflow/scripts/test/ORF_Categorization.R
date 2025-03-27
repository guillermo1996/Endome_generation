#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript classification by ORF status
##
## Aim: add a column in the GTF to represent the
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-03-26
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-03-26)
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
    input = list(gtf = "/home/grocamora/RytenLab-Research/38-Endome_generation/results/AGAT/sq3.annotated.gtf",
               sq3_class = "/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt"),
    output = list(gtf = "/home/grocamora/RytenLab-Research/38-Endome_generation/results/ORF_Categorization/sq3.annotated_orf.gtf"),
    threads=1
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
input_sq3_class_path <- snakemake@input$sq3_class

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
# ---- 1. Load transcriptome and SQ3 Classification ----
gtf <- rtracklayer::import.gff(input_gtf_path)
sq3_class <- readr::read_delim(input_sq3_class_path)


############################################################################## #
# ---- 2. Generate set of categories ----

#----------------------------------------------------------------------------- #
## 2.1 Assign categories to the transcripts ----
sq3_isoforms <- sq3_class %>%
  dplyr::mutate(ORF_cat = dplyr::case_when(
    structural_category == "full-splice_match" & coding == "coding" ~ "known-ORF",
    structural_category != "full-splice_match" & coding == "coding" ~ "novel-ORF",
    coding != "coding" ~ "not-ORF",
    .default = "Unknown"
  )) %>%
  dplyr::select(isoform, ORF_cat, ORF_length)

#----------------------------------------------------------------------------- #
## 2.2 Define the transcript lists ----
tx_known <- sq3_isoforms$isoform[sq3_isoforms$ORF_cat == "known-ORF"]
tx_novel <- sq3_isoforms$isoform[sq3_isoforms$ORF_cat == "novel-ORF"]
tx_not <- sq3_isoforms$isoform[sq3_isoforms$ORF_cat == "not-ORF"]

#----------------------------------------------------------------------------- #
## 2.3 Map transcript ID to ORF Length ----
isoform_length <- sq3_isoforms %>% dplyr::select(isoform, ORF_length) %>% tibble::deframe()


############################################################################## #
# ---- 3. Add information to the GTF ----

#----------------------------------------------------------------------------- #
## 3.1 Add ORF Category ----
mcols(gtf)$ORF_cat <- NA

mcols(gtf)[gtf$type == "transcript" & gtf$transcript_id %in% tx_known, "ORF_cat"] <- "known-ORF"
mcols(gtf)[gtf$type == "transcript" & gtf$transcript_id %in% tx_novel, "ORF_cat"] <- "novel-ORF"
mcols(gtf)[gtf$type == "transcript" & gtf$transcript_id %in% tx_not, "ORF_cat"] <- "not-ORF"


#----------------------------------------------------------------------------- #
## 3.2 Add ORF Length ----
mcols(gtf)$ORF_length <- NA
mcols(gtf)[gtf$ORF_cat %in% c("known-ORF", "novel-ORF"), "ORF_length"] <- isoform_length[mcols(gtf)[gtf$ORF_cat %in% c("known-ORF", "novel-ORF"), "transcript_id"]]


############################################################################## #
# ---- 4. Output the GTF ----
gtf %>% rtracklayer::export(output_gtf_path)


#
#
#
#
#
#
# gtf %>% tibble::as_tibble() %>%
#   dplyr::select(transcript_id, ORF_length) %>%
#   dplyr::filter(!is.na(ORF_length)) %>%
#   dplyr::left_join(sq3_isoforms, by = c("transcript_id" = "isoform")) %>%
#   dplyr::filter(ORF_length.x != ORF_length.y)
#
#
#
# table(lv$type)
# transcripts <- lv[lv$type == "transcript"]
# mcols(transcripts)$custom_attr <- "custom_value"
#
# lv$custom_attr <- NA
# mcols(lv)[lv$type == "transcript", "custom_attr"] <- "custom_value"
# rtracklayer::export(lv, "/home/grocamora/RytenLab-Research/38-Endome_generation/results/AGAT/sq3.annotated_exported.gtf")
#
# # Load necessary libraries
# library(dplyr)
# library(rtracklayer)
#
# # Input files
# gtf_file <- "input.gtf"
# tsv_file <- "annotations.tsv"
# output_file <- "output.gtf"
#
# # Read the TSV file
# tsv_data <- read.delim(tsv_file, header = TRUE, sep = "\t")
#
# # Define ORF categories
# tsv_data <- tsv_data %>%
#   mutate(ORF_category = case_when(
#     coding != "coding" ~ "not-ORF",
#     structural_category == "full-splice_match" & coding == "coding" ~ "known ORF",
#     TRUE ~ "novel ORF"
#   ))
#
# # Read the GTF file as a GRanges object
# gtf_data <- import.gff(gtf_file, format = "gtf")
#
# # Extract transcript features
# gtf_transcripts <- gtf_data[gtf_data$type == "transcript"]
#
# # Extract transcript_id
# gtf_transcripts$transcript_id <- sapply(gtf_transcripts$transcript_id, as.character)
#
# # Merge with annotation data
# gtf_transcripts <- merge(gtf_transcripts, tsv_data, by.x = "transcript_id", by.y = "isoform", all.x = TRUE)
#
# # Append ORF annotation
# gtf_transcripts$ORF_category <- ifelse(is.na(gtf_transcripts$ORF_category), "unknown", gtf_transcripts$ORF_category)
#
# gtf_data$ORF_category <- NA
# gtf_data$ORF_category[gtf_data$type == "transcript"] <- gtf_transcripts$ORF_category
#
# # Export the updated GTF
# export(gtf_data, output_file, format = "gtf")
#
# cat("Updated GTF file saved as", output_file, "\n")
