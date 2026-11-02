#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript categorization by ORF status
##
## Aim: Merge the ORFannotate summary, the ORF protein sequences, the pigeon
## classification and the reference annotation into one isoform summary table,
## adding a subset of it to the GTF as attributes
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S04a-ORF_Categorization
## Date Created: 2025-03-26
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.5 (2026-10-29)
##
## _________________________________________________
##
## Notes:
##    - Based on the ORFannotate v1.0.0 output.
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.5 (2026-10-29): Added gene matching for NIC/NNC transcripts.
##    - v1.4 (2026-10-28): Standardization of the script.
##    - v1.3 (2026-10-20): Accepts Ensembl GTFs.
##    - v1.2 (2026-09-29): Reference transcript, gene and biotypes from the
##      pigeon FSM/ISM matches. in_ref replaced by ref_isoform. Added
##      ref_transcript_id, gene_source and subcategory.
##    - v1.1 (2026-07-16): Support for other transcriptome merging methods.
##    - v1.0 (2025-03-26): Initial version.
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
  prefix <- "Wood.control.st_ref"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(rs, "04-ORF_Identification/ORFannotate", prefix, "ORFannotate_annotated.gtf"),
      ref_annotation = "/home/MinaRyten/Guillermo/Resources/GENCODE/gencode.v48.annotation.gtf",
      protein_fa = file.path(rs, "04-ORF_Identification/ORFannotate", prefix, "protein.fa"),
      orf_summary = file.path(rs, "04-ORF_Identification/ORFannotate", prefix, "ORFannotate_summary.tsv"),
      pigeon_summary = file.path(rs, "03-Artifact_Removal/Classify_Filter", paste0(prefix, ".pigeon_classification.filtered_lite_classification.txt"))
    ),
    output = list(
      gtf = file.path(rs, "temp", paste0(prefix, ".orf_annotated.gtf")),
      isoform_summary = file.path(rs, "temp", paste0(prefix, ".isoform_summary.tsv"))
    ),
    params = list(),
    wildcards = list(
      prefix = prefix
    ),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(GenomicRanges)
  library(rtracklayer)
  library(Biostrings)
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

### FSM transcripts take the reference transcript, gene and biotypes.
ref_isoform_categories <- c("full-splice_match")

### ISM transcripts take the reference gene
ref_gene_categories <- c("full-splice_match", "incomplete-splice_match")

### NIC/NNC transcripts take pigeon's associated gene
ref_novel_categories <- c("novel_in_catalog", "novel_not_in_catalog")

### Isoform summary columns added to the GTF
column_to_gtf <- c("transcript_id", "gene_id", "gene_name", "ref_transcript_id", "structural_category", "subcategory",
                   "ref_isoform", "orfannotate_type", "ref_transcript_type", "ref_gene_type")

### Merge attributes removed from the GTF
columns_to_remove <- c("gene_name", "gene_id", "xloc", "ref_gene_id", "cmp_ref", "class_code", "tss_id", "exon_number",
                       "contained_in", "cmp_ref_gene", "ref_gene_name", "ref_transcript_id")

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' ORFannotate protein FASTA to a tibble
#'
#' @param string_set AAStringSet with ORFannotate headers
#'   ("<transcript_id> gene_id=...;coding_prob=...").
#'
#' @returns Tibble with `transcript_id`, `gene_id`, `coding_prob` and `sequence`.
orfFastaToTibble <- function(string_set){
  data.frame(
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
}

############################################################################## #
# ---- 1. Load the inputs ----
message("Loading the transcriptome, the reference annotation and the ORFannotate/pigeon outputs...")
gtf_df <- rtracklayer::import.gff(smk_inputs$gtf) %>% tibble::as_tibble()
annotation_gtf <- rtracklayer::readGFF(smk_inputs$ref_annotation) %>% tibble::as_tibble()
pigeon_class <- readr::read_delim(smk_inputs$pigeon_summary)

protein_fa <- Biostrings::readAAStringSet(smk_inputs$protein_fa)
orf_summary <- readr::read_tsv(smk_inputs$orf_summary)

############################################################################## #
# ---- 2. Isoform summary ----

#----------------------------------------------------------------------------- #
## 2.1 Select the columns of each input ----
orf_summary_filter <- orf_summary %>%
  dplyr::select(transcript_id, has_orf,
    orf_len = orf_nt_len, utr5_len = utr5_nt_len, utr3_len = utr3_nt_len,
    coding_prob, coding_class, total_junctions, NMD_sensitive)

### Ensembl GTFs use *_biotype instead of *_type
annotation_gtf_filter <- annotation_gtf %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::rename(dplyr::any_of(c(gene_type = "gene_biotype", transcript_type = "transcript_biotype"))) %>%
  dplyr::select(transcript_id, gene_id, gene_type, gene_name, transcript_type)

pigeon_class_filter <- pigeon_class %>%
  dplyr::select(transcript_id = isoform, structural_category, subcategory, associated_gene, associated_transcript, ref_length, ref_exons)

protein_df <- orfFastaToTibble(protein_fa)

#----------------------------------------------------------------------------- #
## 2.2 Join the inputs ----

### ORFannotate writes the merge locus as gene_id of the protein sequences
isoform_summary <- gtf_df %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, merge_gene_id = gene_id) %>%
  dplyr::left_join(orf_summary_filter, by = "transcript_id") %>%
  dplyr::left_join(protein_df %>% dplyr::select(-coding_prob), by = c("transcript_id", "merge_gene_id" = "gene_id")) %>%
  dplyr::left_join(pigeon_class_filter, by = "transcript_id")

### Add reference transcript
isoform_summary <- isoform_summary %>%
  dplyr::mutate(ref_transcript_id = dplyr::case_when(
    transcript_id %in% annotation_gtf_filter$transcript_id ~ transcript_id,
    structural_category %in% ref_gene_categories ~ associated_transcript,
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(-associated_transcript) %>%
  dplyr::left_join(annotation_gtf_filter, by = c("ref_transcript_id" = "transcript_id"))

### Associated NIC/NNC transcripts based on gene_id
annotation_gene_filter <- annotation_gtf_filter %>% dplyr::distinct(gene_id, gene_name, gene_type)

isoform_summary <- isoform_summary %>%
  dplyr::mutate(novel_gene_id = dplyr::if_else(
    is.na(gene_id) & structural_category %in% ref_novel_categories & associated_gene %in% annotation_gene_filter$gene_id,
    associated_gene, NA_character_)) %>%
  dplyr::left_join(annotation_gene_filter, by = c("novel_gene_id" = "gene_id"), suffix = c("", "_novel")) %>%
  dplyr::mutate(
    gene_id = dplyr::coalesce(gene_id, novel_gene_id),
    gene_name = dplyr::coalesce(gene_name, gene_name_novel),
    gene_type = dplyr::coalesce(gene_type, gene_type_novel)) %>%
  dplyr::select(-novel_gene_id, -gene_name_novel, -gene_type_novel, -associated_gene)

#----------------------------------------------------------------------------- #
## 2.3 Clean and sort ----
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
# ---- 3. Add the isoform summary to the GTF ----
gtf_gr <- gtf_df %>%
  dplyr::select(-dplyr::any_of(columns_to_remove)) %>%
  dplyr::left_join(isoform_summary %>% dplyr::select(dplyr::all_of(column_to_gtf)), by = "transcript_id") %>%
  dplyr::rename(ref_tx_type = ref_transcript_type) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

############################################################################## #
# ---- 4. Write the outputs ----
rtracklayer::export(gtf_gr, smk_outputs$gtf, format = "GTF")
readr::write_tsv(isoform_summary, smk_outputs$isoform_summary)
