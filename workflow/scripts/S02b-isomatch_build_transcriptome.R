#!/usr/bin/env Rscript

## _________________________________________________
##
## Isomatch build transcriptome
##
## Aim:  Build the isomatch transcriptome from the merge of the sample
## assemblies
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S02b-isomatch_fix
## Date Created: 2026-07-15
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v2.0 (2026-11-10)
##
## _________________________________________________
##
## Notes:
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v2.0 (2026-11-10): The merge no longer contains GENCODE. GENCODE is added 
##      here (use_reference). 
##    - v1.3 (2026-10-28): Standardization of the script.
##    - v1.2 (2026-09-29): Merged transcripts take the id of the reference
##      transcript they contain. Added merge_transcript_id and min_sample_cnt.
##    - v1.1 (2026-09-25): Added ref_transcript_id, isom_sample_cnt and
##      isom_ref_source.
##    - v1.0 (2026-07-15): Initial version.
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

  rs <- "results/Wood"
  prefix <- "Wood.control.iso_ref"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(rs, "02-Transcriptome_Assembly/isomatch_classify", paste0(prefix, ".annotated.gtf.gz")),
      classification = file.path(rs, "02-Transcriptome_Assembly/isomatch_classify", paste0(prefix, ".classification.txt.gz")),
      present_absent = file.path(rs, "02-Transcriptome_Assembly/isomatch_merge", paste0(prefix, ".present_absent.tsv.gz")),
      ref_gtf = "/home/MinaRyten/Guillermo/Resources/GENCODE/gencode.v48.annotation.gtf"
    ),
    output = list(
      gtf = file.path(rs, "temp", paste0(prefix, ".fixed.gtf"))
    ),
    params = list(
      use_reference = TRUE,
      min_sample_cnt = 1,
      end_match_bp = 50
    ),
    wildcards = list(
      dataset = "Wood",
      group = "control",
      merge_method = "iso_ref"
    ),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(rtracklayer)
  library(plyranges)
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
params <- snakemake@params

use_reference <- as.logical(params$use_reference)
min_sample_cnt <- as.integer(params$min_sample_cnt)
end_match_bp <- as.integer(params$end_match_bp)

### isomatch classify evidence kept as transcript attributes
evidence_cols <- c("diff_to_tss", "diff_to_tes", "within_cage_peak", "dist_to_cage_peak", "within_poly_a_peak", "dist_to_poly_a_site")

### Attributes set on the ### transcript features
tx_cols <- c("gene_name", "ref_transcript_id", "ref_gene_id", "isom_origin", "isom_sample_cnt", "ISOM_CATEGORY", "ISOM_SUBCATEGORY", evidence_cols, "merge_transcript_id")

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Set the "novel" and "NA" values of the ISOM_REF_* attributes to NA
isomRefNA <- function(x) dplyr::if_else(x %in% c("novel", "NA"), NA_character_, x)

#' Count the transcripts of each origin
countOrigin <- function(df){
  df %>%
    dplyr::count(isom_origin) %>%
    dplyr::mutate(label = sprintf("%s %d", isom_origin, n)) %>%
    dplyr::pull(label) %>%
    paste(collapse = ", ")
}

#' Keep the GTF columns shared by the GENCODE and the isomatch records
keepIdCols <- function(gr){
  mcols(gr) <- S4Vectors::DataFrame(
    source = as.character(gr$source),
    type = as.character(gr$type),
    gene_id = gr$gene_id,
    transcript_id = gr$transcript_id
  )
  gr
}

############################################################################## #
# ---- 1. Load the isomatch outputs ----
merge_gr <- rtracklayer::import(smk_inputs$gtf)

### One row per merged transcript, compared with GENCODE
classification <- readr::read_tsv(smk_inputs$classification, col_types = readr::cols(.default = "c"))

### One row per merged transcript, one 0/1 column per merge input file
present_absent <- readr::read_tsv(smk_inputs$present_absent)

### GENCODE transcript and exon features
ref_gr <- rtracklayer::import(smk_inputs$ref_gtf, format = "gtf", feature.type = c("transcript", "exon"),
                              colnames = c("source", "type", "gene_id", "transcript_id", "gene_name"))

############################################################################## #
# ---- 2. Origin of the merged transcripts ----

#----------------------------------------------------------------------------- #
## 2.1 Sample support ----
id_cols <- c("merged_tx_id", "merged_gene_id", "merged_exon_num", "src_tx_count_in_merged_group")
sample_cols <- setdiff(colnames(present_absent), id_cols)

present_absent <- present_absent %>%
  dplyr::mutate(isom_sample_cnt = rowSums(dplyr::across(dplyr::all_of(sample_cols)) > 0))

#----------------------------------------------------------------------------- #
## 2.2 Set origin ----
merged_tx <- classification %>%
  dplyr::transmute(
    merge_transcript_id = isoform_id,
    ISOM_CATEGORY = structural_category,
    ISOM_SUBCATEGORY = subcategory,
    ref_transcript_id = isomRefNA(ref_tx_id),
    ref_gene_id = isomRefNA(ref_gene_id),
    gene_name = isomRefNA(ref_gene_name),
    dplyr::across(dplyr::all_of(evidence_cols), isomRefNA)
  ) %>%
  dplyr::left_join(
    present_absent %>% dplyr::select(merge_transcript_id = merged_tx_id, merge_gene_id = merged_gene_id, isom_sample_cnt),
    by = "merge_transcript_id"
  ) %>%
  dplyr::mutate(
    fsm = ISOM_CATEGORY == "full-splice_match" & !is.na(ref_transcript_id),
    same_ends = abs(as.numeric(diff_to_tss)) <= end_match_bp & abs(as.numeric(diff_to_tes)) <= end_match_bp,
    isom_origin = dplyr::case_when(
      fsm & same_ends %in% TRUE ~ "reference",
      fsm ~ "alt_end",
      TRUE ~ "novel"
    )
  ) %>%
  dplyr::select(-fsm, -same_ends)

message(sprintf("Merged transcripts: %d (%s)", nrow(merged_tx), countOrigin(merged_tx)))

#----------------------------------------------------------------------------- #
## 2.3 Sample filter ----
merged_tx <- merged_tx %>% dplyr::filter(isom_origin == "reference" | isom_sample_cnt >= min_sample_cnt)

message(sprintf("Kept after the sample filter (reference or >= %d samples): %d (%s)", min_sample_cnt, nrow(merged_tx), countOrigin(merged_tx)))

############################################################################## #
# ---- 3. Build the transcriptome ----

#----------------------------------------------------------------------------- #
## 3.1 GENCODE transcripts ----
ref_support <- merged_tx %>%
  dplyr::filter(isom_origin == "reference") %>%
  dplyr::select(merge_transcript_id, transcript_id = ref_transcript_id) %>%
  dplyr::inner_join(present_absent, by = c("merge_transcript_id" = "merged_tx_id")) %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise(
    merge_transcript_id = paste(merge_transcript_id, collapse = ","),
    dplyr::across(dplyr::all_of(sample_cols), ~ any(.x > 0)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(isom_sample_cnt = rowSums(dplyr::across(dplyr::all_of(sample_cols)))) %>%
  dplyr::select(transcript_id, isom_sample_cnt, merge_transcript_id)

### With use_reference, every GENCODE transcript; otherwise only the observed ones
ref_meta <- mcols(ref_gr) %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript", use_reference | transcript_id %in% ref_support$transcript_id) %>%
  dplyr::transmute(transcript_id, gene_name, ref_transcript_id = transcript_id, ref_gene_id = gene_id, isom_origin = "reference") %>%
  dplyr::left_join(ref_support, by = "transcript_id") %>%
  dplyr::mutate(isom_sample_cnt = dplyr::coalesce(isom_sample_cnt, 0))

message(sprintf("GENCODE transcripts: %d (%d observed in the samples)", nrow(ref_meta), sum(ref_meta$isom_sample_cnt > 0)))

#----------------------------------------------------------------------------- #
## 3.2 Alternative-end versions ----

### Numbered within each GENCODE transcript by sample support
alt_meta <- merged_tx %>%
  dplyr::filter(isom_origin == "alt_end") %>%
  dplyr::arrange(ref_transcript_id, dplyr::desc(isom_sample_cnt), merge_transcript_id) %>%
  dplyr::group_by(ref_transcript_id) %>%
  dplyr::mutate(transcript_id = paste0(ref_transcript_id, "_alt", dplyr::row_number())) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene_id = ref_gene_id) %>% 
  dplyr::relocate(transcript_id)

message(sprintf("Alternative-end versions: %d of %d GENCODE transcripts", nrow(alt_meta), dplyr::n_distinct(alt_meta$ref_transcript_id)))

#----------------------------------------------------------------------------- #
## 3.3 Novel transcripts ----

### Keep the isomatch ids
novel_meta <- merged_tx %>%
  dplyr::filter(isom_origin == "novel") %>%
  dplyr::mutate(transcript_id = merge_transcript_id, gene_id = merge_gene_id)

message(sprintf("Novel transcripts: %d", nrow(novel_meta)))

#----------------------------------------------------------------------------- #
## 3.4 Combine ----

### Merged records of the alt_end and novel transcripts, with their new ids
new_ids <- dplyr::bind_rows(alt_meta, novel_meta) %>% dplyr::select(merge_transcript_id, transcript_id, gene_id)

merge_part <- merge_gr[merge_gr$transcript_id %in% new_ids$merge_transcript_id]
id_idx <- match(merge_part$transcript_id, new_ids$merge_transcript_id)
merge_part$gene_id <- new_ids$gene_id[id_idx]
merge_part$transcript_id <- new_ids$transcript_id[id_idx]

ref_part <- ref_gr[ref_gr$transcript_id %in% ref_meta$transcript_id]

gr <- c(keepIdCols(ref_part), keepIdCols(merge_part))

### Transcript attributes
tx_meta <- dplyr::bind_rows(ref_meta, alt_meta, novel_meta) %>% dplyr::select(transcript_id, dplyr::all_of(tx_cols))

stopifnot(!anyDuplicated(tx_meta$transcript_id))

metadata <- mcols(gr) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(tx_meta, by = "transcript_id") %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(tx_cols), ~ dplyr::if_else(type == "transcript", as.character(.x), NA_character_))) %>%
  S4Vectors::DataFrame()
mcols(gr) <- metadata

message(sprintf("Transcriptome: %d transcripts", nrow(tx_meta)))

############################################################################## #
# ---- 4. Output the GTF ----
rtracklayer::export(gr, smk_outputs$gtf, format = "gtf")