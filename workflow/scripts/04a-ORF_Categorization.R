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
##   + in_ref: TRUE when the transcript's gene_id is found in the reference
##   annotation.
##   + A subset of the isoform summary (transcript_id, gene_id, gene_name,
##   structural_category, in_ref, orfannotate_type, ref_transcript_type,
##   ref_gene_type) is joined back onto the GTF as attributes.
##
## - Changelog:
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
  snakemake <- Snakemake(
    input = list(gtf = "data/test_data/ORFannotate_annotated.gtf",
                 ref_annotation = "~/RytenLab-Research/Resources/Gencode/gencode.v48.annotation.gtf",
                 protein_fa = "data/test_data/protein.fa",
                 orf_summary = "data/test_data/ORFannotate_summary.tsv",
                 pigeon_summary = "data/test_data/test.pigeon_classification.filtered_lite_classification.txt"),
    output = list(gtf = "data/test_data/test.orf_annotated.gtf",
                  isoform_summary = "data/test_data/Ebbert.control.isoform_summary.tsv"),
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

## SUGGESTED CHANGE (isomatch merge, pending review) ---------------------------
## With `transcript_merge_method: isomatch`, transcript_id is isomatch's own id
## (ISOMT_*), so the GENCODE join below (by transcript_id) finds no match: every
## transcript gets in_ref = FALSE and NA biotypes, and the `pc` filter in 04b
## removes all of them. 02a-isomatch_fix.R stores the matched reference
## transcript in the `ref_transcript_id` GTF attribute (NA for novel), which
## survives pigeon and ORFannotate. Several ISOMT_* can share one reference id
## (the reference isoform and its alternative-end versions). Suggested:
##
##   gtf_tx <- gtf_df %>% dplyr::filter(type == "transcript")
##   if (!"ref_transcript_id" %in% colnames(gtf_tx)) gtf_tx$ref_transcript_id <- gtf_tx$transcript_id
##
##   isoform_summary <- gtf_tx %>%
##     dplyr::select(transcript_id, ref_transcript_id, merge_gene_id = gene_id) %>%
##     dplyr::left_join(orf_summary_filter, by = "transcript_id") %>%
##     dplyr::left_join(annotation_gtf_filter, by = c("ref_transcript_id" = "transcript_id")) %>%
##     ...  # protein and pigeon joins unchanged (by transcript_id)
##
## and keep `ref_transcript_id` in the isoform summary / relocate() call. For
## StringTie it falls back to transcript_id, so the current behaviour is kept.
## Alternative for any merge method: pigeon's `associated_transcript` column
## (FSM and ISM matches) instead of the GTF attribute.
##
## Optional: the `isom_sample_cnt` / `isom_ref_source` attributes (transcript
## features only) could be carried into the isoform summary to filter
## reference-only transcripts without long-read support
## (isom_sample_cnt == 0) in 04b.
## -----------------------------------------------------------------------------

## Generate the information that will be appended to the isoform summary
orf_summary_filter <- orf_summary %>%
  dplyr::select(transcript_id, has_orf, 
    orf_len = orf_nt_len, utr5_len = utr5_nt_len, utr3_len = utr3_nt_len,
    coding_class, total_junctions, NMD_sensitive)

annotation_gtf_filter <- annotation_gtf %>% 
  dplyr::filter(type == "transcript") %>% 
  dplyr::select(transcript_id, gene_id, gene_type, gene_name, transcript_type)

pigeon_class_filter <- pigeon_class %>% 
  dplyr::select(transcript_id = isoform, structural_category, subcategory, ref_length, ref_exons)

protein_df <- StringSet_to_tibble(protein_fa)

## From the main GTF file, extract the transcripts and merge with ORFannotate
## summary table, the reference annotation, the ORF sequences and the Pigeon
## Classification
isoform_summary <- gtf_df %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, merge_gene_id = gene_id) %>% 
  dplyr::left_join(orf_summary_filter, by = "transcript_id") %>% 
  dplyr::left_join(annotation_gtf_filter, by = "transcript_id") %>% 
  dplyr::left_join(protein_df, by = c("transcript_id", "merge_gene_id" = "gene_id")) %>% 
  dplyr::left_join(pigeon_class_filter, by = c("transcript_id"))

# Clean and sort the data
isoform_summary <- isoform_summary %>% 
  dplyr::rename(orfannotate_type = coding_class, ref_transcript_type = transcript_type, ref_gene_type = gene_type) %>% 
  dplyr::mutate(
    ## SUGGESTED CHANGE (isomatch, pending review): with the join by
    ## ref_transcript_id above, in_ref means "matches a GENCODE transcript".
    ## Decide whether ISM matches should count (pigeon structural_category) or
    ## only full-splice matches, including alternative 5'/3' ends.
    in_ref = !is.na(gene_id),
    gene_id = ifelse(is.na(gene_id), merge_gene_id, gene_id),
    gene_name = ifelse(is.na(gene_id), gene_id, gene_name)
  ) %>% 
  dplyr::rename(orf_seq = sequence) %>% 
  dplyr::relocate(
    transcript_id, gene_id, gene_name, merge_gene_id,
    structural_category, subcategory,
    in_ref, coding_prob, orfannotate_type, ref_transcript_type, ref_gene_type,
    has_orf, ref_length, orf_len, utr5_len, utr3_len, total_junctions, ref_exons
  )

############################################################################## #
# ---- 3. Add relevant information to the gtf ----
column_to_gtf <- c("transcript_id", "gene_id", "gene_name", "structural_category", "in_ref", "orfannotate_type", "ref_transcript_type", "ref_gene_type")
columns_to_remove <- c("gene_name", "gene_id", "xloc", "ref_gene_id", "cmp_ref", "class_code", "tss_id", "exon_number", "contained_in", "cmp_ref_gene", "ref_gene_name")

gtf_gr <- gtf_df %>%
  dplyr::select(-any_of(columns_to_remove)) %>%
  dplyr::left_join(isoform_summary %>% dplyr::select(all_of(column_to_gtf)), by = "transcript_id") %>%
  dplyr::rename(ref_tx_type = ref_transcript_type) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

############################################################################## #
# ---- 4. Output ----
gtf_gr %>% rtracklayer::export(output_gtf_path, format="GTF")
isoform_summary %>% readr::write_tsv(output_tsv_path)
