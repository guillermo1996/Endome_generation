#!/usr/bin/env Rscript

## _________________________________________________
##
## PacBio build transcriptome
##
## Aim: Build the PacBio transcriptome from the pooled isoseq collapse of the
## samples, classified against the reference annotation with pigeon
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S01b-PacBio_Build_Transcriptome
## Date Created: 2026-11-11
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.1 (2026-11-12)
##
## _________________________________________________
##
## Notes:
##    - Only for PacBio HiFi samples (pb* merge presets).
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.1 (2026-11-12): Per-sample read counts from the collapse read_stat.
##    - v1.0 (2026-11-11): Initial version.
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

  rs <- "results/Shimada/01-Alignment"
  prefix <- "Shimada.temporal_cortex.pb_ref"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(rs, "PacBio", "collapse", paste0(prefix, ".collapsed.gtf")),
      classification = file.path(rs, "PacBio", "pigeon", paste0(prefix, ".pigeon_classification.txt")),
      read_stat = file.path(rs, "PacBio", "collapse", paste0(prefix, ".collapsed.read_stat.txt")),
      read_ids = as.list(file.path(rs, "PacBio", "pbmm2", paste0(c("N1_T", "N28_T", "R2_T", "R8_T"), ".read_ids.tsv.gz"))),
      ref_gtf = "/home/MinaRyten/Guillermo/Resources/GENCODE/gencode.v48.annotation.gtf"
    ),
    output = list(
      gtf = file.path(rs, "PacBio", "temp", paste0(prefix, ".pacbio.gtf"))
    ),
    params = list(
      use_reference = TRUE,
      min_tpm = 1.0,
      min_sample_reads = 2,
      min_sample_cnt = 2,
      end_match_bp = 50
    ),
    wildcards = list(
      dataset = "Shimada",
      group = "temporal_cortex",
      merge_method = "pb_ref"
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
min_tpm <- as.numeric(params$min_tpm)
min_sample_reads <- as.integer(params$min_sample_reads)
min_sample_cnt <- as.integer(params$min_sample_cnt)
end_match_bp <- as.numeric(params$end_match_bp)

tx_cols <- c("gene_name", "ref_transcript_id", "ref_gene_id", "pb_sample_cnt", "pb_read_cnt",
             "PIGEON_CATEGORY", "PIGEON_SUBCATEGORY", "diff_to_tss", "diff_to_tes", "merge_transcript_id")
line_end_col <- "pb_origin"

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Count the transcripts of each origin
countOrigin <- function(df){
  df %>%
    dplyr::count(pb_origin) %>%
    dplyr::mutate(label = sprintf("%s %d", pb_origin, n)) %>%
    dplyr::pull(label) %>%
    paste(collapse = ", ")
}

#' Keep the GTF columns shared by the reference and the collapsed records
keepIdCols <- function(gr){
  mcols(gr) <- S4Vectors::DataFrame(
    source = as.character(gr$source),
    type = as.character(gr$type),
    gene_id = gr$gene_id,
    transcript_id = gr$transcript_id
  )
  gr
}

#' Count the samples that support each row
countSupport <- function(df, lib_size){
  counts <- as.matrix(df[names(lib_size)])
  fl_tpm <- sweep(counts, 2, lib_size, "/") * 1e6
  as.integer(rowSums(fl_tpm >= min_tpm & counts >= min_sample_reads))
}


############################################################################## #
# ---- 1. Load the inputs ----
message("Loading the collapsed transcripts, the pigeon classification, the read counts and the reference...")

### Pooled isoseq collapse (PB.<locus>.<n>)
collapse_gr <- rtracklayer::import(smk_inputs$gtf, format = "gtf")

### pigeon classification of the collapsed transcripts against the reference
classification <- readr::read_tsv(smk_inputs$classification, col_types = readr::cols(.default = "c"))

### Sample of each read
read_samples <- purrr::map_dfr(unlist(smk_inputs$read_ids), ~ readr::read_tsv(.x, col_names = c("read_id", "sample"), col_types = "cc"))
if (anyDuplicated(read_samples$read_id) > 0) {
  stop("Read names are not unique across the samples of the group, so reads cannot be assigned to samples")
}

### Collapsed transcript of each read
read_stat <- readr::read_tsv(smk_inputs$read_stat, col_select = c(id, pbid), col_types = readr::cols(.default = "c"))
unassigned <- setdiff(read_stat$id, read_samples$read_id)
if (length(unassigned) > 0) {
  stop(sprintf("%d collapsed reads are not in the alignments of any sample of the group", length(unassigned)))
}

### Full-length reads per collapsed transcript and sample
flnc_count <- read_stat %>%
  dplyr::inner_join(read_samples, by = c("id" = "read_id")) %>%
  dplyr::count(merge_transcript_id = pbid, sample) %>%
  tidyr::pivot_wider(names_from = sample, values_from = n, values_fill = 0)
samples <- setdiff(colnames(flnc_count), "merge_transcript_id")

### Library size of each sample: its full-length reads over every collapsed transcript
lib_size <- colSums(flnc_count[samples])
message(sprintf("Full-length reads per sample: %s", paste(sprintf("%s %d", samples, lib_size), collapse = ", ")))

### Reference transcript and exon features
ref_gr <- rtracklayer::import(smk_inputs$ref_gtf, format = "gtf", feature.type = c("transcript", "exon"),
                              colnames = c("source", "type", "gene_id", "transcript_id", "gene_name"))

ref_tx <- ref_gr[ref_gr$type == "transcript"] %>%
  mcols() %>%
  tibble::as_tibble() %>%
  dplyr::select(ref_transcript_id = transcript_id, ref_gene_id = gene_id, gene_name)

############################################################################## #
# ---- 2. Origin of the collapsed transcripts ----

#----------------------------------------------------------------------------- #
## 2.1 Sample support ----
support <- flnc_count %>%
  dplyr::mutate(
    pb_sample_cnt = countSupport(flnc_count, lib_size),
    pb_read_cnt = rowSums(dplyr::across(dplyr::all_of(samples)))
  )

### Locus of each collapsed transcript
collapse_genes <- collapse_gr[collapse_gr$type == "transcript"] %>%
  mcols() %>%
  tibble::as_tibble() %>%
  dplyr::distinct(merge_transcript_id = transcript_id, merge_gene_id = gene_id)

#----------------------------------------------------------------------------- #
## 2.2 Set origin ----
### reference: same intron chain as a reference transcript, both ends within end_match_bp
### alt_end: same intron chain, at least one end further away
### novel: every other pigeon category
collapsed_tx <- classification %>%
  dplyr::transmute(
    merge_transcript_id = isoform,
    PIGEON_CATEGORY = structural_category,
    PIGEON_SUBCATEGORY = subcategory,
    ref_transcript_id = dplyr::if_else(structural_category == "full-splice_match", associated_transcript, NA_character_),
    diff_to_tss = -as.numeric(diff_to_TSS),
    diff_to_tes = -as.numeric(diff_to_TTS)
  ) %>%
  dplyr::left_join(ref_tx, by = "ref_transcript_id") %>%
  dplyr::left_join(collapse_genes, by = "merge_transcript_id") %>%
  dplyr::left_join(support %>% dplyr::select(merge_transcript_id, pb_sample_cnt, pb_read_cnt), by = "merge_transcript_id") %>%
  dplyr::mutate(
    fsm = !is.na(ref_gene_id),
    same_ends = abs(diff_to_tss) <= end_match_bp & abs(diff_to_tes) <= end_match_bp,
    pb_origin = dplyr::case_when(
      fsm & same_ends %in% TRUE ~ "reference",
      fsm ~ "alt_end",
      TRUE ~ "novel"
    ),
    ref_transcript_id = dplyr::if_else(fsm, ref_transcript_id, NA_character_)
  ) %>%
  dplyr::select(-fsm, -same_ends)

message(sprintf("Collapsed transcripts: %d (%s)", nrow(collapsed_tx), countOrigin(collapsed_tx)))

#----------------------------------------------------------------------------- #
## 2.3 Sample filter ----
collapsed_tx <- collapsed_tx %>%
  dplyr::filter(pb_origin == "reference" | pb_sample_cnt >= min_sample_cnt)

message(sprintf("Kept after the sample filter (reference, or FL-TPM >= %s with >= %d reads in >= %d samples): %d (%s)",
                min_tpm, min_sample_reads, min_sample_cnt, nrow(collapsed_tx), countOrigin(collapsed_tx)))

############################################################################## #
# ---- 3. Build the transcriptome ----

#----------------------------------------------------------------------------- #
## 3.1 Reference transcripts ----
### Samples and reads over every collapsed transcript that matches the reference transcript
ref_support <- collapsed_tx %>%
  dplyr::filter(pb_origin == "reference") %>%
  dplyr::select(merge_transcript_id, transcript_id = ref_transcript_id) %>%
  dplyr::inner_join(flnc_count, by = "merge_transcript_id") %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise(
    merge_transcript_id = paste(merge_transcript_id, collapse = ","),
    dplyr::across(dplyr::all_of(samples), sum),
    .groups = "drop"
  )

ref_support <- ref_support %>%
  dplyr::mutate(
    pb_sample_cnt = countSupport(ref_support, lib_size),
    pb_read_cnt = rowSums(dplyr::across(dplyr::all_of(samples)))
  ) %>%
  dplyr::select(transcript_id, pb_sample_cnt, pb_read_cnt, merge_transcript_id)

### With use_reference, every reference transcript; otherwise only the observed ones
ref_meta <- ref_tx %>%
  dplyr::filter(use_reference | ref_transcript_id %in% ref_support$transcript_id) %>%
  dplyr::transmute(transcript_id = ref_transcript_id, gene_name, ref_transcript_id, ref_gene_id, pb_origin = "reference") %>%
  dplyr::left_join(ref_support, by = "transcript_id") %>%
  dplyr::mutate(
    pb_sample_cnt = dplyr::coalesce(pb_sample_cnt, 0),
    pb_read_cnt = dplyr::coalesce(pb_read_cnt, 0)
  )

message(sprintf("Reference transcripts: %d (%d observed in the samples)", nrow(ref_meta), sum(ref_meta$pb_sample_cnt > 0)))

#----------------------------------------------------------------------------- #
## 3.2 Alternative-end versions ----
alt_meta <- collapsed_tx %>%
  dplyr::filter(pb_origin == "alt_end") %>%
  dplyr::arrange(ref_transcript_id, dplyr::desc(pb_sample_cnt), dplyr::desc(pb_read_cnt), merge_transcript_id) %>%
  dplyr::group_by(ref_transcript_id) %>%
  dplyr::mutate(transcript_id = paste0(ref_transcript_id, "_alt", dplyr::row_number())) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene_id = ref_gene_id) %>%
  dplyr::relocate(transcript_id)

message(sprintf("Alternative-end versions: %d of %d reference transcripts", nrow(alt_meta), dplyr::n_distinct(alt_meta$ref_transcript_id)))

#----------------------------------------------------------------------------- #
## 3.3 Novel transcripts ----
### Keep the collapse ids (PB.<locus>.<n>)
novel_meta <- collapsed_tx %>%
  dplyr::filter(pb_origin == "novel") %>%
  dplyr::mutate(transcript_id = merge_transcript_id, gene_id = merge_gene_id)

message(sprintf("Novel transcripts: %d", nrow(novel_meta)))

#----------------------------------------------------------------------------- #
## 3.4 Combine ----

### Collapsed records of the alt_end and novel transcripts, with their new ids
new_ids <- dplyr::bind_rows(alt_meta, novel_meta) %>% dplyr::select(merge_transcript_id, transcript_id, gene_id)

collapse_part <- collapse_gr[collapse_gr$transcript_id %in% new_ids$merge_transcript_id]
id_idx <- match(collapse_part$transcript_id, new_ids$merge_transcript_id)
collapse_part$gene_id <- new_ids$gene_id[id_idx]
collapse_part$transcript_id <- new_ids$transcript_id[id_idx]
collapse_part$source <- "PacBio"

ref_part <- ref_gr[ref_gr$transcript_id %in% ref_meta$transcript_id]

gr <- c(keepIdCols(ref_part), keepIdCols(collapse_part))

### Transcript attributes
tx_meta <- dplyr::bind_rows(ref_meta, alt_meta, novel_meta) %>%
  dplyr::select(transcript_id, dplyr::all_of(c(tx_cols, line_end_col)))

stopifnot(!anyDuplicated(tx_meta$transcript_id))

### tx_cols on transcript lines only
metadata <- mcols(gr) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(tx_meta, by = "transcript_id") %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(tx_cols), ~ dplyr::if_else(type == "transcript", as.character(.x), NA_character_))) %>%
  dplyr::relocate(dplyr::all_of(line_end_col), .after = dplyr::last_col()) %>%
  S4Vectors::DataFrame()
mcols(gr) <- metadata

message(sprintf("Transcriptome: %d transcripts", nrow(tx_meta)))

############################################################################## #
# ---- 4. Output the GTF ----
dir.create(dirname(smk_outputs$gtf), recursive = TRUE, showWarnings = FALSE)
rtracklayer::export(gr, smk_outputs$gtf, format = "gtf")
