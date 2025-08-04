#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript classification by ORF status
##
## Aim: Complete the GTF with ORF information
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-03-26
##
## Copyright (c) Guillermo Rocamora Pérez, 2025
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
    input = list(gtf = "~/RytenLab-Research/snakefile-refactor/results_k15/04-ORF_Identification/ORFannotate/pigeon.annotated/ORFannotate_annotated_clean.gtf",
                 ref_annotation = "~/RytenLab-Research/Resources/GENCODE/gencode.v48.annotation.gtf",
                 protein_fa = "~/RytenLab-Research/snakefile-refactor/results_k15/04-ORF_Identification/ORFannotate/pigeon.annotated/protein.fa",
                 orf_summary = "~/RytenLab-Research/snakefile-refactor/results_k15/04-ORF_Identification/ORFannotate/pigeon.annotated/ORFannotate_summary.tsv",
                 pigeon_summary = "~/RytenLab-Research/snakefile-refactor/results_k15/03-Pigeon/Classify_Filter/pigeon.annotated_classification.filtered_lite_classification.txt"),
    output = list(gtf = "~/RytenLab-Research/snakefile-refactor/results_k15/04-ORF_Identification/ORF_Category/pigeon.annotated_orf.gtf",
                  isoform_summary = "~/RytenLab-Research/snakefile-refactor/results_k15/04-ORF_Identification/ORF_Category/pigeon.annotated_isoform_summary.tsv"),
    threads = 1
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
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
input_reference_gtf_path <- snakemake@input$ref_annotation
input_protein_seqs <- snakemake@input$protein_fa
input_orf_summary <- snakemake@input$orf_summary
input_pigeon_class <- snakemake@input$pigeon_summary

### Output Paths
output_gtf_path <- snakemake@output$gtf
output_tsv_path <- snakemake@output$isoform_summary

#### Create output directory
dir.create(dirname(output_gtf_path), showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.3 Script Parameters ----
#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----
#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----

StringSet_to_tibble <- function(string_set){
  string_set_df <- data.frame(
    header = names(string_set),
    sequence = as.character(string_set)
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      transcript_id = str_extract(header, "^[^ ]+"),
      gene_id = str_extract(header, "gene_id=([^;]+)") %>% str_remove("gene_id="),
      coding_prob = str_extract(header, "coding_prob=([0-9.]+)") %>% str_remove("coding_prob=") %>% as.numeric()
    ) %>%
    dplyr::select(transcript_id, gene_id, coding_prob, sequence)

  return(string_set_df)
}

############################################################################## #
# ---- 1. Load transcriptome and SQ3 Classification ----
gtf_df <- rtracklayer::import.gff(input_gtf_path) %>% tibble::as_tibble()
annotation_gtf <- rtracklayer::readGFF(input_reference_gtf_path) %>% tibble::as_tibble()
pigeon_class <- readr::read_delim(input_pigeon_class)

protein_fa <- Biostrings::readAAStringSet(input_protein_seqs)
orf_summary <- readr::read_tsv(input_orf_summary)

############################################################################## #
# ---- 2. Generate the Isoform information table ----

## GR: From the main GTF file, extract the transcripts and merge with
## ORFannotate summary table
isoform_summary <- gtf_df %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, stringtie_gene_id = gene_id, gene_name) %>%
  dplyr::left_join(orf_summary %>%
                     dplyr::select(transcript_id, has_orf, orf_len = orf_nt_len, utr5_len = utr5_nt_len, utr3_len = utr3_nt_len, coding_class, junction_count, NMD_sensitive),
                   by = "transcript_id")

## GR: Append information from the reference annotation.
isoform_summary <- isoform_summary %>%
  dplyr::left_join(annotation_gtf %>%
                     dplyr::filter(type == "transcript") %>%
                     dplyr::select(transcript_id, gene_id, gene_type, transcript_type),
                   by = "transcript_id") %>%
  dplyr::relocate(transcript_id, gene_id, gene_name, stringtie_gene_id, has_orf, orf_len, utr5_len, utr3_len, orfannotate_type = coding_class, ref_transcript_type = transcript_type, ref_gene_type = gene_type) %>%
  dplyr::mutate(in_ref = !is.na(gene_id),
                gene_id = ifelse(is.na(gene_id), stringtie_gene_id, gene_id),
                gene_name = ifelse(is.na(gene_name), gene_id, gene_name))

## GR: Append the ORF sequences
protein_df <- StringSet_to_tibble(protein_fa)
isoform_summary <- isoform_summary %>%
  dplyr::left_join(protein_df, by = c("transcript_id", "stringtie_gene_id" = "gene_id")) %>%
  dplyr::rename(orf_seq = sequence)

## GR: Append the Pigeon classification
isoform_summary <- isoform_summary %>%
  dplyr::left_join(pigeon_class %>% dplyr::select(isoform, structural_category, subcategory, ref_length, ref_exons),
                   by = c("transcript_id" = "isoform")) %>%
  dplyr::relocate(transcript_id, gene_id, gene_name, stringtie_gene_id,
                  structural_category, subcategory,
                  in_ref, has_orf, ref_length, orf_len, utr5_len, utr3_len, junction_count, ref_exons,
                  orfannotate_type, ref_transcript_type, ref_gene_type)


############################################################################## #
# ---- 3. Add relevant information to the gtf ----
column_to_gtf <- c("transcript_id", "gene_id", "gene_name", "structural_category", "in_ref", "orfannotate_type", "ref_transcript_type", "ref_gene_type")
NA_cols <- column_to_gtf[column_to_gtf != "transcript_id"]

gtf_gr <- gtf_df %>%
  dplyr::select(-gene_name, -gene_id) %>%
  dplyr::left_join(isoform_summary %>% dplyr::select(all_of(column_to_gtf)), by = "transcript_id") %>%
  # dplyr::mutate(across(all_of(NA_cols), ~if_else(type == "transcript", ., NA))) %>%
  dplyr::rename(ref_tx_type = ref_transcript_type) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)


############################################################################## #
# ---- 4. Output ----
gtf_gr %>% rtracklayer::export(output_gtf_path, format="GTF")
isoform_summary %>% readr::write_tsv(output_tsv_path)
