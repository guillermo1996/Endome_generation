#!/usr/bin/env Rscript

## _________________________________________________
##
## Export the ORF protein set for MMseqs2
##
## Aim: Write the deduplicated protein FASTA used by the all-vs-all MMseqs2
## search, and the map of its protein ids to transcripts
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S05b-Export_Proteins
## Date Created: 2026-08-19
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.1 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - One FASTA per (prefix, orf_filter): every truncation of the same ORF
##      filter shares the MMseqs2 search.
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.1 (2026-10-28): Standardization of the script.
##    - v1.0 (2026-08-19): Initial version.
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
  orf_filter <- "pc_novel"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      protein_fa = file.path(rs, "04-ORF_Identification/ORFannotate", prefix, "protein.fa"),
      gtf_filter = file.path(rs, "04-ORF_Identification/ORF_Filtration", sprintf("%s.%s.orf_filter.gtf", prefix, orf_filter))
    ),
    output = list(
      faa = file.path(rs, "temp", sprintf("%s.%s.proteins.faa", prefix, orf_filter)),
      protein_map = file.path(rs, "temp", sprintf("%s.%s.protein_map.tsv", prefix, orf_filter))
    ),
    params = list(),
    wildcards = list(
      prefix = prefix,
      orf_filter = orf_filter
    ),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(Biostrings)
  library(rtracklayer)
  library(digest)
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
wc <- snakemake@wildcards

run_id <- paste(wc$prefix, wc$orf_filter, sep = ".")

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----
lib_dir <- file.path(snakemake@scriptdir, "../lib")
source(file.path(lib_dir, "Phf-DuckDB.R"))
source(file.path(lib_dir, "Phf-Bin_Information_Score.R"))

############################################################################## #
# ---- 1. Load the protein set ----
message("Reading the ORFannotate protein FASTA...")
protein_fa_df <- StringSetToTibble(Biostrings::readAAStringSet(smk_inputs$protein_fa))

message("Reading the ORF-filtered GTF...")
kept_transcripts <- rtracklayer::import(smk_inputs$gtf_filter) %>%
  tibble::as_tibble() %>%
  dplyr::distinct(transcript_id) %>%
  dplyr::filter(!is.na(transcript_id)) %>%
  dplyr::pull(transcript_id)

message(sprintf("  %s proteins | %s transcripts retained by the '%s' filter in '%s'",
  format(nrow(protein_fa_df), big.mark = ","),
  format(length(kept_transcripts), big.mark = ","),
  wc$orf_filter, wc$prefix))

############################################################################## #
# ---- 2. Assign the protein ids ----

### Remove the trailing stop character
proteins_df <- protein_fa_df %>%
  dplyr::filter(transcript_id %in% kept_transcripts) %>%
  dplyr::filter(!is.na(sequence), nzchar(sequence)) %>%
  dplyr::mutate(
    run_id = run_id,
    aa_seq = sub("\\*+$", "", sequence),
    aa_len = nchar(aa_seq)
  ) %>%
  dplyr::filter(aa_len > 0) %>%
  dplyr::mutate(aa_id = assignId(., "aa_seq", prefix = "aa")) %>%
  dplyr::relocate(aa_id, aa_len, transcript_id, aa_seq, run_id)

message(sprintf("  %s distinct protein sequences", format(dplyr::n_distinct(proteins_df$aa_id), big.mark = ",")))

############################################################################## #
# ---- 3. Write the outputs ----

#----------------------------------------------------------------------------- #
## 3.1 Protein FASTA ----
message("Writing ", smk_outputs$faa, " ...")
protein_fasta <- proteins_df %>% dplyr::distinct(aa_id, aa_seq)
seqs <- Biostrings::AAStringSet(setNames(protein_fasta$aa_seq, protein_fasta$aa_id))
Biostrings::writeXStringSet(seqs, smk_outputs$faa)

#----------------------------------------------------------------------------- #
## 3.2 Protein to transcript map ----
proteins_df %>%
  dplyr::distinct(run_id, aa_id, transcript_id) %>%
  readr::write_tsv(smk_outputs$protein_map)

message("Done.")
