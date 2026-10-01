## _________________________________________________
##
## Export the ORF protein set for MMseqs2
##
## Aim: Produce the deduplicated protein FASTA that the all-vs-all MMseqs2
## search runs on
##
## Project: ENDome generation - Bin Information Content
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-08-19
##
## Latest Version: v1.0 (2026-08-19)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##
## Changelog:
##    - v1.0 (2026-08-19): Initial version.
##
## Contact: guillermorocamora@gmail.com
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
  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  test_orf_filter <- "ref_pc"

  snakemake <- Snakemake(
    input = list(
      protein_fa = file.path(test_dir, "protein.fa"),
      gtf_filter = file.path(test_dir, paste0("test.", test_orf_filter, ".orf_filter.gtf"))
    ),
    output = list(
      faa = file.path(test_dir, paste0("interactive.", test_orf_filter, ".proteins.faa")),
      protein_map = file.path(test_dir, paste0("interactive.", test_orf_filter, ".protein_map.tsv"))
    ),
    wildcards = list(
      prefix = basename(test_dir),
      orf_filter = test_orf_filter
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

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

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
### Snakemake exposes the directory the script itself lives in, so `lib/`
### resolves whether these files sit here or in workflow/scripts/.
lib_dir <- file.path(snakemake@scriptdir, "../lib")
source(file.path(lib_dir, "Phf-DuckDB.R"))
source(file.path(lib_dir, "Phf-Bin_Information_Score.R"))

############################################################################## #
# ---- 1. Load the protein set ----
message("Reading the ORFannotate protein FASTA...")
protein_fa_df <- StringSetToTibble(Biostrings::readAAStringSet(smk_inputs$protein_fa))

message("Reading the pre-truncated GTF to find retained transcripts...")
pre_gtf <- rtracklayer::import(smk_inputs$gtf_filter)

kept_transcripts <- pre_gtf %>%
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
## 3.1 The FASTA handed to MMseqs2 ----
message("Writing ", smk_outputs$faa, " ...")

protein_fasta <- proteins_df %>% 
  dplyr::select(aa_id, aa_seq) %>% 
  dplyr::filter(!is.na(aa_seq), nzchar(aa_seq)) %>% 
  dplyr::distinct()
seqs <- Biostrings::AAStringSet(setNames(protein_fasta$aa_seq, protein_fasta$aa_id))

if(!dir.exists(dirname(smk_outputs$faa))) dir.create(dirname(smk_outputs$faa), showWarnings = FALSE)
Biostrings::writeXStringSet(seqs, smk_outputs$faa)

#----------------------------------------------------------------------------- #
## 3.2 Store aa_id -> transcript_id ----
proteins_df %>% 
  dplyr::distinct(run_id, aa_id, transcript_id) %>% 
  readr::write_tsv(smk_outputs$protein_map)

message("Done.")