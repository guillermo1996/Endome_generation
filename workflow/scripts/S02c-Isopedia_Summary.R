#!/usr/bin/env Rscript

## _________________________________________________
##
## Isopedia summary
##
## Aim: Summarise the per-sample read evidence reported by `isopedia isoform`
## for every transcript of a step 02 assembly
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S02c-Isopedia_Summary
## Date Created: 2026-10-22
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.2 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - Written for the 13-field FORMAT of Isopedia v1.7.3:
##      RC_FSM_JC:RC_FSM_JC_TSS:RC_FSM_JC_TES:RC_FSM_JC_TSS_TES:RC_FSM_JC_EXACT:
##      RC_FSM_JC_WOBBLE_ONLY:RC_EM_ISM:CPM_FSM_JC:CPM_FSM_JC_TSS_TES:CPM_EST:
##      FRAC_FSM_JC_TSS:FRAC_FSM_JC_TES:FRAC_FSM_JC_TSS_TES
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.2 (2026-10-28): Standardization of the script.
##    - v1.1 (2026-10-28): dplyr/dbplyr instead of SQL queries. Uncompressed
##      outputs.
##    - v1.0 (2026-10-22): Initial version.
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
      resources = 'list',
      log = 'list',
      threads = 'numeric',
      scriptdir = 'character')
  )

  rs <- "results/Wood"
  prefix <- "Wood.control.st_ref"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      tsv = file.path(rs, "02-Transcriptome_Assembly/isopedia", paste0(prefix, ".isopedia.tsv.gz"))
    ),
    output = list(
      samples = file.path(rs, "temp", paste0(prefix, ".isopedia.samples.tsv")),
      transcripts = file.path(rs, "temp", paste0(prefix, ".isopedia.transcripts.tsv")),
      groups = file.path(rs, "temp", paste0(prefix, ".isopedia.groups.tsv"))
    ),
    params = list(
      min_read = 1,
      group_by = "tissue"
    ),
    wildcards = list(
      dataset = "Wood",
      group = "control",
      merge_method = "st_ref"
    ),
    resources = list(mem_gb = 16),
    threads = 4,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(DBI)
  library(duckdb)
  library(dbplyr)
  library(conflicted)
  library(tidyverse)
})

## Package options
options(dplyr.summarise.inform = FALSE)
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

isopedia_tsv <- smk_inputs$tsv
min_read <- as.numeric(params$min_read)
group_by <- as.character(params$group_by)
mem_gb <- snakemake@resources$mem_gb

### FORMAT fields used in the summary, named by their output column
format_fields <- c(
  fsm = "RC_FSM_JC",
  fsm_tss = "RC_FSM_JC_TSS",
  fsm_tes = "RC_FSM_JC_TES",
  fsm_tss_tes = "RC_FSM_JC_TSS_TES",
  em = "RC_EM_ISM",
  cpm_est = "CPM_EST"
)

############################################################################## #
# ---- 1. Parse the header ----

### The file starts with "##[SAMPLE]" lines (sample metadata), one "##[DBINFO]"
### line (total evidence per sample) and the "#" column header. The first
### transcript row is also read, as it holds the FORMAT field names.
message("Reading the header of ", isopedia_tsv, " ...")
con_txt <- gzfile(isopedia_tsv, "r")
header_lines <- c()
line <- readLines(con_txt, n = 1)

while (startsWith(line, "#")) {
  header_lines <- c(header_lines, line)
  line <- readLines(con_txt, n = 1)
}
first_row <- line
close(con_txt)

#----------------------------------------------------------------------------- #
## 1.1 Columns and FORMAT ----
col_line <- tail(grep("^#[^#]", header_lines, value = TRUE), 1)
columns <- strsplit(sub("^#", "", col_line), "\t", fixed = TRUE)[[1]]
format_idx <- match("FORMAT", columns)
if (is.na(format_idx)) stop("No FORMAT column in the Isopedia header", call. = FALSE)

fixed_cols <- columns[seq_len(format_idx - 1)]
sample_names <- columns[(format_idx + 1):length(columns)]
sample_keys <- sprintf("s_%05d", seq_along(sample_names))

first_fields <- strsplit(first_row, "\t", fixed = TRUE)[[1]]
if (length(first_fields) != length(columns)) {
  stop(sprintf("The first row has %i fields but the header has %i columns",
               length(first_fields), length(columns)), call. = FALSE)
}
format_names <- strsplit(first_fields[format_idx], ":", fixed = TRUE)[[1]]
missing_fields <- setdiff(format_fields, format_names)
if (length(missing_fields) > 0) {
  stop("FORMAT field(s) not found: ", paste(missing_fields, collapse = ", "),
       ". This script expects the output of Isopedia v1.7.3.", call. = FALSE)
}
field_pos <- match(format_fields, format_names) %>% setNames(names(format_fields))
zero_cell <- paste(rep("0", length(format_names)), collapse = ":")

fixed_clean <- fixed_cols %>%
  dplyr::recode("n_pos_samples/sample_size" = "isopedia_n_pos") %>%
  gsub("[^A-Za-z0-9_]+", "_", .)

#----------------------------------------------------------------------------- #
## 1.2 Sample metadata ----
sample_lines <- sub("^##\\[SAMPLE\\]", "", grep("^##\\[SAMPLE\\]", header_lines, value = TRUE))
if (length(sample_lines) > 1) {
  sample_meta <- readr::read_tsv(I(paste(sample_lines, collapse = "\n")), quote = "",
                                 col_types = readr::cols(.default = "c"))
  ### The metadata "sample" column is the BioSample accession. The first column
  ### (the index sample name) takes its place.
  colnames(sample_meta)[colnames(sample_meta) == "sample"] <- "biosample"
  colnames(sample_meta)[1] <- "sample"
} else {
  sample_meta <- tibble::tibble(sample = sample_names)
}

### Total evidence per sample: "name:count,name:count,..."
dbinfo_line <- sub("^##\\[DBINFO\\]", "", grep("^##\\[DBINFO\\]", header_lines, value = TRUE))
dbinfo <- tibble::tibble(entry = unlist(strsplit(dbinfo_line, ",", fixed = TRUE))) %>%
  tidyr::extract(entry, into = c("sample", "total_evidence"), regex = "^(.*):(\\d+)$", convert = TRUE)

### Sample group: case-insensitive match of `group_by` to a metadata column
group_col <- colnames(sample_meta)[tolower(colnames(sample_meta)) == tolower(group_by)][1]
if (is.na(group_col)) {
  warning(sprintf("Column '%s' not in the sample metadata (%s): every sample goes to group 'all'", group_by, paste(colnames(sample_meta), collapse = ", ")), call. = FALSE)
}

samples <- tibble::tibble(sample_key = sample_keys, sample = sample_names) %>%
  dplyr::left_join(sample_meta, by = "sample") %>%
  dplyr::left_join(dbinfo, by = "sample") %>%
  dplyr::mutate(
    sample_group = if (is.na(group_col)) "all" else .data[[group_col]],
    sample_group = dplyr::coalesce(dplyr::na_if(stringr::str_squish(tolower(sample_group)), ""), "unknown"),
    .after = sample
  )

############################################################################## #
# ---- 2. Read the evidence with DuckDB ----

### In-memory database that spills to tmp_dir when the memory limit is reached
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
tmp_dir <- file.path(dirname(smk_outputs$transcripts), paste0(".duckdb_tmp_", Sys.getpid()))
DBI::dbExecute(con, sprintf("SET temp_directory = '%s'", tmp_dir))
DBI::dbExecute(con, sprintf("SET threads = %i", as.integer(snakemake@threads)))
if (!is.null(mem_gb)) DBI::dbExecute(con, sprintf("SET memory_limit = '%iGB'", floor(0.8 * as.numeric(mem_gb))))
DBI::dbExecute(con, "SET preserve_insertion_order = false")

samples_db <- dplyr::copy_to(con, dplyr::select(samples, sample_key, sample_group), "samples")

read_file <- sprintf(
  "read_csv('%s', delim = '\t', header = false, skip = %i, quote = '', escape = '',
            auto_detect = false, max_line_size = 100000000, columns = {%s})",
  isopedia_tsv, length(header_lines),
  paste(sprintf("'%s': 'VARCHAR'", c(fixed_clean, "FORMAT", sample_keys)), collapse = ", ")
)
raw <- dplyr::tbl(con, dplyr::sql(paste("SELECT * FROM", read_file)))

### One row per transcript x sample
cells <- dplyr::tbl(con, dplyr::sql(sprintf(
  "UNPIVOT (SELECT trans_id, COLUMNS('^s_[0-9]+$') FROM %s)
   ON COLUMNS('^s_[0-9]+$') INTO NAME sample_key VALUE val",
  read_file)))

split_fields <- lapply(field_pos, function(i) rlang::expr(as.double(split_part(val, ":", !!i))))

#----------------------------------------------------------------------------- #
## 2.1 Evidence per transcript x sample group ----

### Cells without evidence are dropped before splitting the FORMAT fields, so
### the long table is never materialised
message("Aggregating the evidence per transcript and sample group (first pass over the file)...")
tx_groups <- cells %>%
  dplyr::filter(val != !!zero_cell) %>%
  dplyr::mutate(!!!split_fields) %>%
  dplyr::inner_join(samples_db, by = "sample_key") %>%
  dplyr::group_by(trans_id, sample_group) %>%
  dplyr::summarise(
    n_pos_fsm = sum(dplyr::if_else(fsm >= !!min_read, 1L, 0L), na.rm = TRUE),
    n_pos_est = sum(dplyr::if_else(fsm + em >= !!min_read, 1L, 0L), na.rm = TRUE),
    dplyr::across(c(fsm, fsm_tss, fsm_tes, fsm_tss_tes, em), ~ sum(.x, na.rm = TRUE)),
    cpm_est_sum = sum(cpm_est, na.rm = TRUE),
    cpm_est_max = max(cpm_est, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::collect()

#----------------------------------------------------------------------------- #
## 2.2 Transcript information ----
message("Reading the transcript columns (second pass over the file)...")
tx_info <- raw %>%
  dplyr::transmute(
    chrom = dplyr::if_else(substr(chrom, 1L, 3L) == "chr", chrom, paste0("chr", chrom)),
    start = as.integer(start) + 1L,   # Isopedia starts are 0-based
    end = as.integer(end),
    tx_length = as.integer(length),
    exon_count = as.integer(exon_count),
    trans_id, gene_id, detected, isopedia_n_pos
  ) %>%
  dplyr::collect()

DBI::dbDisconnect(con, shutdown = TRUE)
unlink(tmp_dir, recursive = TRUE)

n_dup <- sum(duplicated(tx_info$trans_id))
if (n_dup > 0) stop(n_dup, " duplicated trans_id in the Isopedia output", call. = FALSE)

############################################################################## #
# ---- 3. Summaries ----
group_sizes <- dplyr::count(samples, sample_group, name = "n_samples_group")
n_samples <- nrow(samples)

#----------------------------------------------------------------------------- #
## 3.1 Per transcript x sample group ----
groups_out <- tx_groups %>%
  dplyr::inner_join(group_sizes, by = "sample_group") %>%
  dplyr::transmute(
    trans_id, sample_group, n_samples_group,
    n_pos_fsm, n_pos_est,
    fsm, fsm_tss, fsm_tes, fsm_tss_tes, em,
    cpm_est_mean = cpm_est_sum / n_samples_group,
    cpm_est_max
  ) %>%
  dplyr::arrange(trans_id, sample_group)

#----------------------------------------------------------------------------- #
## 3.2 Per transcript ----
tx_evidence <- tx_groups %>%
  dplyr::group_by(trans_id) %>%
  dplyr::summarise(
    n_groups_pos_fsm = sum(n_pos_fsm > 0),
    n_pos_fsm = sum(n_pos_fsm),
    n_pos_est = sum(n_pos_est),
    dplyr::across(c(fsm, fsm_tss, fsm_tes, fsm_tss_tes, em), sum),
    cpm_est_mean = sum(cpm_est_sum) / n_samples,
    cpm_est_max = max(cpm_est_max),
    .groups = "drop"
  )

### Transcripts without evidence in any sample get zeros
count_cols <- c("n_pos_fsm", "n_pos_est", "n_groups_pos_fsm", "fsm", "fsm_tss", "fsm_tes",
                "fsm_tss_tes", "em", "cpm_est_mean", "cpm_est_max")

transcripts_out <- tx_info %>%
  dplyr::left_join(tx_evidence, by = "trans_id") %>%
  dplyr::mutate(
    n_samples = n_samples,
    dplyr::across(dplyr::all_of(count_cols), ~ tidyr::replace_na(.x, 0)),
    frac_tss = dplyr::if_else(fsm > 0, fsm_tss / fsm, NA_real_),
    frac_tes = dplyr::if_else(fsm > 0, fsm_tes / fsm, NA_real_),
    frac_tss_tes = dplyr::if_else(fsm > 0, fsm_tss_tes / fsm, NA_real_)
  ) %>%
  dplyr::select(
    chrom:isopedia_n_pos, n_samples, n_pos_fsm, n_pos_est, n_groups_pos_fsm,
    fsm, fsm_tss, fsm_tes, fsm_tss_tes, em, frac_tss, frac_tes, frac_tss_tes,
    cpm_est_mean, cpm_est_max
  ) %>%
  dplyr::arrange(chrom, start, trans_id)

############################################################################## #
# ---- 4. Write the outputs ----
readr::write_tsv(transcripts_out, smk_outputs$transcripts)
readr::write_tsv(groups_out, smk_outputs$groups)
readr::write_tsv(dplyr::select(samples, -sample_key), smk_outputs$samples)

n_detected <- sum(transcripts_out$n_pos_fsm > 0)
message(sprintf("Done: %s transcripts, %s (%.1f%%) with full-splice-match reads in >= 1 sample.",
  format(nrow(transcripts_out), big.mark = ","), format(n_detected, big.mark = ","),
  100 * n_detected / nrow(transcripts_out)))
