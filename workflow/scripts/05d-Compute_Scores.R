## _________________________________________________
##
## Bin Information Content scores
##
## Aim: Score every bin of the ENDome by how much information is lost when its
##      member transcripts are collapsed. Each within-bin transcript pair is
##      given a Gower distance over the ORF/protein features, and the bin's score
##      aggregates those pairwise distances.
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
##    - Mutates the DuckDB built in 05c-Build_ENDome_DB.R
##    - Two scores are emitted. `info_score` normalises by the comparable pair
##      mass, which under uniform member weights cancels the bin size exactly. 
##      `info_score_rao` is Rao's quadratic entropy proper, which does penalise large heterogeneous bins
##
## Changelog:
##    - v1.0 (2026-08-19): Initial version. Ported from T03-ComputeScores.R.
##
## Contact: guillermorocamora@gmail.com
## _________________________________________________

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
  ## dataset.group.merge_method folder created by the `test_data` rule (its
  ## name is the prefix wildcard), and the orf_filter preset and truncation
  ## (width, txEnd) to test. This script WRITES into the DuckDB, so it works on
  ## an "interactive" copy of the linked database (multi-GB, copied once), never
  ## on the link itself, which points to the pipeline result.
  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  test_orf_filter <- "ref_pc"
  test_width <- "500"
  test_txEnd <- "3p"
  test_db <- paste0(test_orf_filter, ".w", test_width, ".", test_txEnd, ".duckdb")
  interactive_db <- file.path(test_dir, paste0("interactive.", test_db))
  if (!file.exists(interactive_db)) file.copy(file.path(test_dir, paste0("test.", test_db)), interactive_db)
  snakemake <- Snakemake(
    input = list(
      duckdb = interactive_db,
      pairs = file.path(test_dir, paste0("test.", test_orf_filter, ".pairs.tsv")),
      protein_map = file.path(test_dir, paste0("test.", test_orf_filter, ".protein_map.tsv"))
    ),
    output = list(
      done = paste0(interactive_db, ".scores.done")
    ),
    params = list(
      include_superseded = TRUE,
      compute_sim_pid = FALSE,
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
      prefix = basename(test_dir),
      orf_filter = test_orf_filter,
      width = test_width,
      txEnd = test_txEnd
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

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

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
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
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
## 1.2 Per-transcript features ----
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
# ---- 2. Import the MMseqs2 hits ----
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

# ### Drop the proteins this ENDome has no transcript for
# DBI::dbExecute(con, "
#   DELETE FROM mmseqs_raw
#   WHERE query  NOT IN (SELECT aa_id FROM proteins)
#      OR target NOT IN (SELECT aa_id FROM proteins)")

db <- getTbls(con)
message(sprintf("  %s hits retained", format(dplyr::pull(dplyr::count(db$mmseqs_raw), n), big.mark = ",")))

### Drop comparisons that are never going to be made in this dataset
# aa_gene_map <- readr::read_tsv(smk_inputs$protein_map) %>% dplyr::distinct(aa_id, gene_id)
# DBI::dbWriteTable(con, "aa_gene_map", aa_gene_map, overwrite = TRUE)

# n_before <- DBI::dbGetQuery(con, "SELECT count(*) n FROM mmseqs_raw")$n
# DBI::dbExecute(con, "
#   DELETE FROM mmseqs_raw
#   WHERE NOT EXISTS (
#     SELECT 1
#     FROM aa_gene_map q JOIN aa_gene_map t USING (gene_id)
#     WHERE q.aa_id = mmseqs_raw.query AND t.aa_id = mmseqs_raw.target)")
# n_after <- DBI::dbGetQuery(con, "SELECT count(*) n FROM mmseqs_raw")$n
# DBI::dbExecute(con, "DROP TABLE IF EXISTS aa_gene_map")
# message(sprintf("  %s of %s hits dropped as cross-gene",
# format(n_before - n_after, big.mark = ","), format(n_before, big.mark = ",")))

#----------------------------------------------------------------------------- #
## 2.1 Pairwise protein similarity ----
selfBits <- db$mmseqs_raw %>%
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
  dplyr::left_join(dplyr::rename(selfBits, aa_id_1 = aa_id, bits_1 = bits_self), by = "aa_id_1") %>%
  dplyr::left_join(dplyr::rename(selfBits, aa_id_2 = aa_id, bits_2 = bits_self), by = "aa_id_2") %>%
  dplyr::mutate(
    run_id = run_id,
    bsr_min = pmin(1, bits / pmin(bits_1, bits_2)),
    bsr_max = pmin(1, bits / pmax(bits_1, bits_2)),
    bsr_mixed = pmin(1, bits / sqrt(as.numeric(bits_1) * as.numeric(bits_2)))
  )

DBI::dbExecute(con, "DROP TABLE IF EXISTS protein_pairs")
protein_pairs <- protein_pairs %>% dplyr::compute(name = "protein_pairs", temporary = TRUE)
db <- getTbls(con)

message(sprintf("protein_pairs: %s pairs", format(dplyr::pull(dplyr::count(db$protein_pairs), n), big.mark = ",")))

#----------------------------------------------------------------------------- #
## 2.2 Average PID similarity per bin ----
if(isTRUE(params$compute_sim_pid)){
  shhh({ library(BiocParallel); library(pwalign); library(Biostrings) })

  message("Computing the within-bin PID similarity (this is slow)...")
  bin_proteins <- createBinProteins(db) %>% dplyr::collect()
  bpparam_multi <- BiocParallel::SnowParam(workers = snakemake@threads, progressbar = FALSE, type = "SOCK")

  bin_list <- bin_proteins %>%
    dplyr::select(bin_id, aa_id, aa_seq) %>%
    dplyr::collect() %>%
    S4Vectors::split(., .$bin_id)

  # Group the bins into clusters to reduce the dispatching overhead
  n_clusters <- max(1, BiocParallel::bpnworkers(bpparam_multi) * 2)
  bin_clusters <- split(bin_list, ceiling(seq_along(bin_list)/(length(bin_list)/n_clusters)))

  # Measure the similarity of every bin of every cluster
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
# ---- 3. Per-variable distances ----
message("Enumerating the within-bin transcript pairs and their distances...")

searched_aa <- readr::read_tsv(smk_inputs$protein_map) %>% dplyr::distinct(aa_id)

bin_pairs <- bin_members %>%
  dplyr::select(bin_id, transcript_id, w) %>% 
  dplyr::inner_join(x=., y=., by = "bin_id", suffix = c("_i", "_j"), relationship = "many-to-many") %>%
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

### Alternative metrics not implemented
# d_junctions = relativeDelta(total_junctions_i, total_junctions_j)
# d_utr_free = relativeDelta(utr_free_len_i, utr_free_len_j)
# d_utr_binned = relativeDelta(utr_binned_len_i, utr_binned_len_j)
# d_structural_cat = discordance(structural_category_i, structural_category_j)
# d_subcategory = discordance(subcategory_i, subcategory_j)
# d_biotype = pmax(discordance(ref_gene_type_i, ref_gene_type_j), discordance(ref_transcript_type_i, ref_transcript_type_j), na.rm = TRUE)

#----------------------------------------------------------------------------- #
## 3.1 Gower distance of every pair ----
component_cols <- grep("^d_", colnames(bin_pairs_complete), value = TRUE)

### Include asymmetric variant for external purpose
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
    starts_with("transcripts"), d_ij, d_ij_asymmetric, 
    all_of(component_cols), starts_with("aa_id"), w_i, w_j
  )

#----------------------------------------------------------------------------- #
## 3.2 Report the comparability of each component ----
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
  tidyr::pivot_longer(everything(), names_to = c("variable", ".value"),
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


message("Normalized info_score summary:")
print(summary(bin_information$info_score_norm))
message("Default info_score_rao summary:")
print(summary(bin_information$info_score))

############################################################################## #
# ---- 5. Store in the database ----
writeDuckTable(
  con, "bin_information",
  spec = list(
    df = bin_information,
    pk = "bin_id",
    fk = c(bin_id = "bins.bin_id")
  )
)

message("Done.")
DBI::dbDisconnect(con, shutdown = TRUE)
