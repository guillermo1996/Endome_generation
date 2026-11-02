#!/usr/bin/env Rscript

## _________________________________________________
##
## Bin Information Content scores
##
## Aim: Score every bin of the ENDome by how much information is lost when its
## member transcripts are collapsed, from the Gower distances of its
## transcript pairs over the ORF/protein features
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S05d-Compute_Scores
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
##    - Writes the `bin_information` table into the DuckDB built by S05c.
##    - `info_score` is 1 - Rao's quadratic entropy, which penalises large
##      heterogeneous bins. `info_score_norm` normalises by the mass of the
##      comparable pairs, which cancels the bin size under uniform weights.
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
  width <- "500"
  txEnd <- "3p"
  endome <- sprintf("%s.%s.w%s.%s", prefix, orf_filter, width, txEnd)

  ### This script writes into the DuckDB: work on a copy, never on the result
  duckdb_copy <- file.path(rs, "temp", sprintf("%s.duckdb", endome))
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)
  if (!file.exists(duckdb_copy)) file.copy(file.path(rs, "05-Truncation/DuckDB", sprintf("%s.duckdb", endome)), duckdb_copy)

  snakemake <- Snakemake(
    input = list(
      duckdb = duckdb_copy,
      pairs = file.path(rs, "05-Truncation/MMseqs2", sprintf("%s.%s.pairs.tsv", prefix, orf_filter)),
      protein_map = file.path(rs, "05-Truncation/MMseqs2", sprintf("%s.%s.protein_map.tsv", prefix, orf_filter))
    ),
    output = list(
      done = file.path(rs, "temp", sprintf("%s.scores.done", endome))
    ),
    params = list(
      protein_metric = "bsr_max",
      length_scaling = "ratio",
      member_weights = "uniform",
      component_weights = list(
        protein = 1,
        cds_len = 0.1,
        start_codon = 0.2,
        stop_codon = 0.2,
        nmd = 0.4)
    ),
    wildcards = list(
      prefix = prefix,
      orf_filter = orf_filter,
      width = width,
      txEnd = txEnd
    ),
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
wc <- snakemake@wildcards

run_id <- paste(wc$prefix, wc$orf_filter, paste0("w", wc$width), wc$txEnd, sep = ".")

component_weights <- unlist(params$component_weights)
protein_metric <- params$protein_metric
length_scaling <- params$length_scaling

if (!identical(params$member_weights, "uniform")) {
  stop("Only `member_weights = \"uniform\"` is implemented.", call. = FALSE)
}

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----
lib_dir <- file.path(snakemake@scriptdir, "../lib")
source(file.path(lib_dir, "Phf-DuckDB.R"))
source(file.path(lib_dir, "Phf-Bin_Information_Score.R"))

############################################################################## #
# ---- 1. Load the database elements ----
message("Connecting to ", smk_inputs$duckdb, " ...")
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = smk_inputs$duckdb)
db <- getTbls(con)

#----------------------------------------------------------------------------- #
## 1.1 Bin membership ----
bin_members <- db$bins %>%
  dplyr::left_join(createBinMembership(db), by = "bin_id") %>%
  dplyr::select(bin_id, transcript_id, n_members) %>%
  dplyr::mutate(w = 1/n_members) %>%
  dplyr::compute(name = "tmp_bin_members", temporary = TRUE)

bin_members_local <- bin_members %>% dplyr::collect()

message(sprintf("%s bins | %s members | %s singleton bins | largest bin: %s members",
  format(dplyr::n_distinct(bin_members_local$bin_id), big.mark = ","),
  format(nrow(bin_members_local), big.mark = ","),
  format(sum(bin_members_local$n_members == 1), big.mark = ","),
  max(bin_members_local$n_members)))

#----------------------------------------------------------------------------- #
## 1.2 Transcript features ----
message("Extracting the transcript features...")

tx_features <- db$transcripts %>%
  dplyr::left_join(
    dplyr::select(db$cds, cds_id, cds_len, start_codon_canonical, stop_codon_canonical,
                  start_codon_sequence, stop_codon_sequence),
    by = "cds_id") %>%
  dplyr::left_join(dplyr::select(db$utr5, utr5_id, utr5_len), by = "utr5_id") %>%
  dplyr::left_join(dplyr::select(db$utr3, utr3_id, utr3_len), by = "utr3_id") %>%
  dplyr::transmute(
    transcript_id, aa_id, cds_len, total_junctions,
    structural_category, subcategory, ref_gene_type, ref_transcript_type,
    start_codon_canonical, stop_codon_canonical,
    start_codon_sequence, stop_codon_sequence,
    nmd_sensitive = NMD_sensitive
  ) %>%
  dplyr::semi_join(bin_members, by = "transcript_id") %>%
  dplyr::compute(name = "tmp_tx_features", temporary = TRUE)

############################################################################## #
# ---- 2. Protein similarity ----
message("Loading the MMseqs2 hits from ", smk_inputs$pairs, " ...")
DBI::dbExecute(con, sprintf("
  CREATE OR REPLACE TEMPORARY TABLE mmseqs_raw AS
  SELECT query, target, fident, alnlen, qcov, tcov, raw, bits, evalue
  FROM read_csv('%s', delim = '\t', header = true,
    types = {
      'query':  'VARCHAR', 'target': 'VARCHAR',
      'fident': 'FLOAT',   'alnlen': 'INTEGER',
      'qcov':   'FLOAT',   'tcov':   'FLOAT',
      'raw':    'INTEGER', 'bits':   'INTEGER',
      'evalue': 'DOUBLE'
    })", smk_inputs$pairs))

db <- getTbls(con)
message(sprintf("  %s hits loaded", format(dplyr::pull(dplyr::count(db$mmseqs_raw), n), big.mark = ",")))

#----------------------------------------------------------------------------- #
## 2.1 Pairwise protein similarity ----

### Bit-score ratios against the self-hit of each protein
self_bits <- db$mmseqs_raw %>%
  dplyr::filter(query == target) %>%
  dplyr::group_by(aa_id = query) %>%
  dplyr::summarise(bits_self = max(bits, na.rm = TRUE), .groups = "drop")

protein_pairs <- db$mmseqs_raw %>%
  dplyr::filter(query != target) %>%
  dplyr::mutate(aa_id_1 = pmin(query, target), aa_id_2 = pmax(query, target)) %>%
  dplyr::group_by(aa_id_1, aa_id_2) %>%
  dplyr::slice_max(bits, n = 1, with_ties = FALSE, na_rm = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::select(aa_id_1, aa_id_2, fident, qcov, tcov, bits) %>%
  dplyr::left_join(dplyr::rename(self_bits, aa_id_1 = aa_id, bits_1 = bits_self), by = "aa_id_1") %>%
  dplyr::left_join(dplyr::rename(self_bits, aa_id_2 = aa_id, bits_2 = bits_self), by = "aa_id_2") %>%
  dplyr::mutate(
    run_id = run_id,
    bsr_min = pmin(1, bits / pmin(bits_1, bits_2)),
    bsr_max = pmin(1, bits / pmax(bits_1, bits_2)),
    bsr_mixed = pmin(1, bits / sqrt(as.numeric(bits_1) * as.numeric(bits_2)))
  )

DBI::dbExecute(con, "DROP TABLE IF EXISTS protein_pairs")
protein_pairs <- protein_pairs %>% dplyr::compute(name = "protein_pairs", temporary = TRUE)
db <- getTbls(con)

message(sprintf("  %s protein pairs", format(dplyr::pull(dplyr::count(db$protein_pairs), n), big.mark = ",")))

#----------------------------------------------------------------------------- #
## 2.2 Average PID similarity per bin ----
if(isTRUE(params$compute_sim_pid)){
  shhh({ library(BiocParallel); library(pwalign); library(Biostrings) })

  message("Computing the within-bin PID similarity (slow)...")
  bin_proteins <- createBinProteins(db) %>% dplyr::collect()
  bpparam_multi <- BiocParallel::SnowParam(workers = snakemake@threads, progressbar = FALSE, type = "SOCK")

  bin_list <- bin_proteins %>%
    dplyr::select(bin_id, aa_id, aa_seq) %>%
    S4Vectors::split(., .$bin_id)

  ### Bins are sent to the workers in clusters to reduce the dispatch overhead
  n_clusters <- max(1, BiocParallel::bpnworkers(bpparam_multi) * 2)
  bin_clusters <- split(bin_list, ceiling(seq_along(bin_list)/(length(bin_list)/n_clusters)))

  helper_path <- file.path(lib_dir, "Phf-Bin_Information_Score.R")
  bin_similarity_pid <- BiocParallel::bplapply(bin_clusters, function(bin_cluster, helper_path){
    source(helper_path)
    purrr::map(bin_cluster, measureSeqSimilarity, min_width = 30, max_width = 7500, type = "global") %>% dplyr::bind_rows()
  }, helper_path = helper_path, BPPARAM = bpparam_multi) %>% dplyr::bind_rows()

  writeDuckTable(
    con, "bin_pid_similarity",
    spec = list(
      df = bin_similarity_pid %>% dplyr::mutate(run_id = run_id, .before = 0),
      pk = "bin_id",
      fk = c(bin_id = "bins.bin_id")
    )
  )
}

############################################################################## #
# ---- 3. Pairwise distances ----
message("Enumerating the within-bin transcript pairs and their distances...")

### Proteins left out of the MMseqs2 search have no comparable similarity
searched_aa <- readr::read_tsv(smk_inputs$protein_map) %>% dplyr::distinct(aa_id)

bin_pairs <- bin_members %>%
  dplyr::select(bin_id, transcript_id, w) %>%
  dplyr::inner_join(x = ., y = ., by = "bin_id", suffix = c("_i", "_j"), relationship = "many-to-many") %>%
  dplyr::filter(transcript_id_i < transcript_id_j) %>%
  dplyr::left_join(dplyr::rename_with(tx_features, ~ paste0(.x, "_i"), -transcript_id), by = c("transcript_id_i" = "transcript_id")) %>%
  dplyr::left_join(dplyr::rename_with(tx_features, ~ paste0(.x, "_j"), -transcript_id), by = c("transcript_id_j" = "transcript_id")) %>%
  dplyr::mutate(aa_id_1 = pmin(aa_id_i, aa_id_j), aa_id_2 = pmax(aa_id_i, aa_id_j)) %>%
  dplyr::left_join(
    db$protein_pairs %>% dplyr::select(aa_id_1, aa_id_2, protein_sim = !!sym(protein_metric)),
    by = c("aa_id_1", "aa_id_2")) %>%
  dplyr::collect() %>%
  dplyr::left_join(searched_aa %>% dplyr::mutate(searched_i = TRUE), by = c("aa_id_i" = "aa_id")) %>%
  dplyr::left_join(searched_aa %>% dplyr::mutate(searched_j = TRUE), by = c("aa_id_j" = "aa_id")) %>%
  dplyr::mutate(protein_sim = dplyr::case_when(
    aa_id_i == aa_id_j ~ 1,
    is.na(searched_i) | is.na(searched_j) ~ NA_real_,
    .default = dplyr::coalesce(protein_sim, 0)
  ))

bin_pairs_complete <- bin_pairs %>%
  dplyr::transmute(
    bin_id, transcript_id_i, transcript_id_j,
    d_protein = 1 - protein_sim,
    d_cds_len = relativeDelta(cds_len_i, cds_len_j, length_scaling),
    d_start_codon = codonDissimilarity(start_codon_canonical_i, start_codon_canonical_j,
                                       start_codon_sequence_i, start_codon_sequence_j),
    d_stop_codon = codonDissimilarity(stop_codon_canonical_i, stop_codon_canonical_j,
                                      stop_codon_sequence_i, stop_codon_sequence_j),
    d_nmd = discordance(nmd_sensitive_i, nmd_sensitive_j),
    aa_id_i, aa_id_j, w_i, w_j,
    start_codon_canonical_i, start_codon_canonical_j,
    stop_codon_canonical_i, stop_codon_canonical_j,
    nmd_sensitive_i, nmd_sensitive_j
  )
message(sprintf("%s within-bin transcript pairs scored", format(nrow(bin_pairs_complete), big.mark = ",")))

#----------------------------------------------------------------------------- #
## 3.1 Gower distance of every pair ----
component_cols <- grep("^d_", colnames(bin_pairs_complete), value = TRUE)

### Asymmetric variant: shared canonical codons and shared NMD insensitivity
### are not counted as agreement
bin_pairs_distance <- bin_pairs_complete %>%
  dplyr::mutate(
    d_ij = combineGower(., component_weights),
    d_ij_asymmetric = combineGower(
      dplyr::mutate(.,
        d_start_codon = dplyr::if_else(start_codon_canonical_i & start_codon_canonical_j, NA_real_, d_start_codon),
        d_stop_codon = dplyr::if_else(stop_codon_canonical_i & stop_codon_canonical_j, NA_real_, d_stop_codon),
        d_nmd = dplyr::if_else(!nmd_sensitive_i & !nmd_sensitive_j, NA_real_, d_nmd)),
      component_weights),
    .after = transcript_id_j
  ) %>%
  dplyr::mutate(run_id = run_id, .before = 0) %>%
  dplyr::select(run_id, bin_id,
    dplyr::starts_with("transcript"), d_ij, d_ij_asymmetric,
    dplyr::all_of(component_cols), dplyr::starts_with("aa_id"), w_i, w_j
  )

#----------------------------------------------------------------------------- #
## 3.2 Comparability of each component ----
component_report <- bin_pairs_distance %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(component_cols),
      list(
        comparable = \(x) mean(!is.na(x)),
        mean = \(x) mean(x, na.rm = TRUE)
      )
    )
  ) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = c("variable", ".value"),
                      names_pattern = "d_(.*)_(comparable|mean)") %>%
  dplyr::mutate(weight = component_weights[variable], .after = variable)
print(component_report)

############################################################################## #
# ---- 4. Bin Information Content ----

#----------------------------------------------------------------------------- #
## 4.1 Aggregate the pairwise distances ----
bin_rao <- bin_pairs_distance %>%
  dplyr::group_by(bin_id) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    n_pairs_scored = sum(!is.na(d_ij)),
    q_rao = 2 * sum(w_i * w_j * d_ij, na.rm = TRUE),
    pair_mass = 2 * sum(w_i * w_j * !is.na(d_ij)),
    q_rao_asymmetric = 2 * sum(w_i * w_j * d_ij_asymmetric, na.rm = TRUE),
    pair_mass_asymmetric = 2 * sum(w_i * w_j * !is.na(d_ij_asymmetric)),
    dplyr::across(dplyr::all_of(component_cols),
      \(x) sum(w_i * w_j * x, na.rm = TRUE) / sum(w_i * w_j * !is.na(x)), .names = "mean_{.col}"),
    .groups = "drop"
  )

#----------------------------------------------------------------------------- #
## 4.2 Bin information scores ----

### Singleton bins have no pairs and score 1
bin_information <- bin_members_local %>%
  dplyr::distinct(bin_id, n_members) %>%
  dplyr::left_join(bin_rao, by = "bin_id") %>%
  dplyr::mutate(
    dplyr::across(c(n_pairs, n_pairs_scored), \(x) dplyr::coalesce(x, 0)),
    dplyr::across(c(q_rao, pair_mass, q_rao_asymmetric, pair_mass_asymmetric), \(x) dplyr::coalesce(x, 0)),
    dplyr::across(dplyr::starts_with("mean_d_"), \(x) dplyr::if_else(is.nan(x), NA_real_, x)),
    q_rao_norm = dplyr::if_else(pair_mass > 0, q_rao / pair_mass, 0),
    info_score = 1 - q_rao,
    info_score_norm = 1 - q_rao_norm,
    info_score_asymmetric = 1 - dplyr::if_else(pair_mass_asymmetric > 0, q_rao_asymmetric / pair_mass_asymmetric, 0),
    run_id = run_id
  ) %>%
  dplyr::select(-pair_mass_asymmetric) %>%
  dplyr::relocate(run_id, bin_id, n_members, info_score_norm, info_score, info_score_asymmetric)

message("info_score_norm summary:")
print(summary(bin_information$info_score_norm))
message("info_score summary:")
print(summary(bin_information$info_score))

############################################################################## #
# ---- 5. Write the scores into the DuckDB ----
writeDuckTable(
  con, "bin_information",
  spec = list(
    df = bin_information,
    pk = "bin_id",
    fk = c(bin_id = "bins.bin_id")
  )
)

DBI::dbDisconnect(con, shutdown = TRUE)
message("Done.")
