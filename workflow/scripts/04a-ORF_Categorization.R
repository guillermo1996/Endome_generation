#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript classification by ORF status
##
## Author: Guillermo Rocamora Pérez
## Date Created: 2025-03-26
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.1 (2026-07-16)
## _________________________________________________
##
## - Notes:
##   + Based on ORFannotate v1.0.0 output.
##   + Merges the ORFannotate summary, reference annotation, ORF protein
##   sequences and Pigeon classification into one per-transcript isoform
##   summary table.
##   + ref_transcript_id: the reference transcript of each transcript. The
##   transcript itself when it is a reference transcript, otherwise pigeon's
##   associated_transcript for full-splice (FSM) and incomplete-splice (ISM)
##   matches. NA for every other category.
##   + ref_isoform: TRUE for pigeon FSM transcripts, i.e. same intron chain as
##   a reference transcript (the ends may differ). Replaces in_ref.
##   + ref_transcript_type: biotype of the reference transcript, FSM only.
##   + ref_gene_type: biotype of the reference gene.
##   + A subset of the isoform summary (transcript_id, gene_id, gene_name,
##   structural_category, subcategory, ref_isoform, orfannotate_type,
##   ref_transcript_type, ref_gene_type) is joined back onto the GTF as
##   attributes.
##
## - Changelog:
##   + v1.2 (2026-09-29): reference transcript, gene and biotypes from pigeon's
##   FSM/ISM matches; in_ref replaced by ref_isoform; added ref_transcript_id,
##   gene_source and subcategory (GTF); gene_name is NA for merge loci.
##   + v1.1 (2026-07-16): adapted to work with different transcript assembly
##   merging software.
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
  ## dataset.group.merge_method folder created by the `test_data` rule. Outputs
  ## use the "interactive" prefix so they never overwrite the linked results.
  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  snakemake <- Snakemake(
    input = list(gtf = file.path(test_dir, "ORFannotate_annotated.gtf"),
                 ref_annotation = "/home/MinaRyten/Guillermo/Resources/GENCODE/gencode.v48.annotation.gtf",
                 protein_fa = file.path(test_dir, "protein.fa"),
                 orf_summary = file.path(test_dir, "ORFannotate_summary.tsv"),
                 pigeon_summary = file.path(test_dir, "test.pigeon_classification.filtered_lite_classification.txt")),
    output = list(gtf = file.path(test_dir, "interactive.orf_annotated.gtf"),
                  isoform_summary = file.path(test_dir, "interactive.isoform_summary.tsv")),
    threads = 1
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(Biostrings)
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

## FSM transcripts will contain reference genes and biotypes, while ISM will only take the reference gene
ref_isoform_categories <- c("full-splice_match")
ref_gene_categories <- c("full-splice_match", "incomplete-splice_match")

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

## Generate the information that will be appended to the isoform summary
orf_summary_filter <- orf_summary %>%
  dplyr::select(transcript_id, has_orf, 
    orf_len = orf_nt_len, utr5_len = utr5_nt_len, utr3_len = utr3_nt_len,
    coding_class, total_junctions, NMD_sensitive)

annotation_gtf_filter <- annotation_gtf %>% 
  dplyr::filter(type == "transcript") %>% 
  dplyr::select(transcript_id, gene_id, gene_type, gene_name, transcript_type)

pigeon_class_filter <- pigeon_class %>% 
  dplyr::select(transcript_id = isoform, structural_category, subcategory, associated_transcript, ref_length, ref_exons)

protein_df <- StringSet_to_tibble(protein_fa)

## From the main GTF file, extract the transcripts and merge with ORFannotate
## summary table, the ORF sequences and the Pigeon Classification. The protein
## sequences are joined by the merge locus, as written by ORFannotate.
isoform_summary <- gtf_df %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, merge_gene_id = gene_id) %>% 
  dplyr::left_join(orf_summary_filter, by = "transcript_id") %>% 
  dplyr::left_join(protein_df, by = c("transcript_id", "merge_gene_id" = "gene_id")) %>% 
  dplyr::left_join(pigeon_class_filter, by = "transcript_id")

## Reference transcript: the transcript itself when it is in the reference,
## otherwise pigeon's match for FSM/ISM. Its gene and biotypes come from the
## reference annotation.
isoform_summary <- isoform_summary %>%
  dplyr::mutate(ref_transcript_id = dplyr::case_when(
    transcript_id %in% annotation_gtf_filter$transcript_id ~ transcript_id,
    structural_category %in% ref_gene_categories ~ associated_transcript,
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(-associated_transcript) %>%
  dplyr::left_join(annotation_gtf_filter, by = c("ref_transcript_id" = "transcript_id"))

# Clean and sort the data
isoform_summary <- isoform_summary %>% 
  dplyr::rename(orfannotate_type = coding_class, orf_seq = sequence) %>% 
  dplyr::mutate(
    ref_isoform = structural_category %in% ref_isoform_categories,
    gene_source = dplyr::if_else(is.na(gene_id), "merge", "reference"),
    gene_id = dplyr::coalesce(gene_id, merge_gene_id),
    ref_transcript_type = dplyr::if_else(ref_isoform, transcript_type, NA_character_),
    ref_gene_type = gene_type
  ) %>% 
  dplyr::select(-transcript_type, -gene_type) %>% 
  dplyr::relocate(
    transcript_id, gene_id, gene_name, gene_source, merge_gene_id,
    structural_category, subcategory, ref_isoform, ref_transcript_id,
    coding_prob, orfannotate_type, ref_transcript_type, ref_gene_type,
    has_orf, ref_length, orf_len, utr5_len, utr3_len, total_junctions, ref_exons
  )

############################################################################## #
# ---- 3. Add relevant information to the gtf ----
column_to_gtf <- c("transcript_id", "gene_id", "gene_name", "ref_transcript_id", "structural_category", "subcategory", "ref_isoform", "orfannotate_type", "ref_transcript_type", "ref_gene_type")
columns_to_remove <- c("gene_name", "gene_id", "xloc", "ref_gene_id", "cmp_ref", "class_code", "tss_id", "exon_number", "contained_in", "cmp_ref_gene", "ref_gene_name", "ref_transcript_id")

gtf_gr <- gtf_df %>%
  dplyr::select(-any_of(columns_to_remove)) %>%
  dplyr::left_join(isoform_summary %>% dplyr::select(all_of(column_to_gtf)), by = "transcript_id") %>%
  dplyr::rename(ref_tx_type = ref_transcript_type) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

############################################################################## #
# ---- 4. Output ----
gtf_gr %>% rtracklayer::export(output_gtf_path, format="GTF")
isoform_summary %>% readr::write_tsv(output_tsv_path)
