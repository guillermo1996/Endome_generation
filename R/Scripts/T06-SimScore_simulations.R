## _________________________________________________
##
## Simulations of Similarity Scores
##
## Aim:
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-05-20
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-05-20)
##
## _________________________________________________
##
## - Notes:
##
## - Changelog:
##
## - Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

############################################################################## #
# ---- 0. Setup ----

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(Biostrings))
shhh(library(BiocParallel))
shhh(library(pwalign))
shhh(library(patchwork))
shhh(library(tidyverse))

## Package options - these are some examples and do not modify the script usage.
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)


#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
utr_site <- "3p"
n_threads <- 64

#----------------------------------------------------------------------------- #
## 0.3 Script Paths ----
main_path <- "~/RytenLab-Research/38-Endome_generation" # here::here()
script_path <- file.path(main_path, "R")

### Input Paths
sq3_class_path<- file.path(main_path, "results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt")
table_mapping_path <- file.path(script_path, paste0("Results/Table_mapping_", utr_site))

tx_tbl_path <- file.path(table_mapping_path, "tx_tbl.rds")
tr_tbl_path <- file.path(table_mapping_path, "tr_tbl.rds")
bin_tbl_path <- file.path(table_mapping_path, "bin_tbl.rds")
orf_tbl_path <- file.path(table_mapping_path, "orf_tbl.rds")

map_tx_tr_path <- file.path(table_mapping_path, "map_tx_tr.rds")
map_tr_bin_path <- file.path(table_mapping_path, "map_tr_bin.rds")

### Output Paths
results_path <- file.path(script_path, paste0("Results/SimScore_Simulations_", utr_site))

simScore_path <- function(file_prefix, seed = ""){
  if(seed != "") seed <- paste0("_s", seed)
  return(file.path(results_path, paste0(file_prefix, "_SimScore", seed, ".rds")))
}

plot_default_simScore_path <- file.path(script_path, "figures/SimScore_Simulations", paste0("default_simScore_", utr_site, ".png"))
plot_rl_simScore_path <- file.path(script_path, "figures/SimScore_Simulations", paste0("randomLabel_simScore_", utr_site, ".png"))

#### Create output directory
dir.create(dirname(simScore_path("")), showWarnings = F, recursive = T)
dir.create(dirname(plot_default_simScore_path), showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----

#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----
source(file.path(script_path, "R/hf_graph_and_themes.R"))
source(file.path(script_path, "R/Phf-Endome/hf_sim_score.R"))

############################################################################## #
# ---- 1. Section 1 of the script ----

#----------------------------------------------------------------------------- #
## 1.1 Load previous tables ----
tx_tbl <- readRDS(tx_tbl_path)
orf_tbl <- readRDS(orf_tbl_path)
# tr_tbl <- readRDS(tr_tbl_path)
# bin_tbl <- readRDS(bin_tbl_path)º
# sq3_class <- vroom::vroom(sq3_class_path)

############################################################################## #
# ---- 2. Base results simulation ----

#----------------------------------------------------------------------------- #
## 2.1 Run the Similarity Scoring ----
default_simScore_path <- simScore_path(file_prefix = "notRandom")
if(!file_test("-f", default_simScore_path)){
  # Validate input table
  validate_df(tx_tbl, c("bin_id", "orf_id", "orf_seq"))

  # Split the transcript table into bin groups for parallel processing
  bin_list <- tx_tbl %>%
    dplyr::select(bin_id, orf_id, orf_seq) %>%
    split(., .$bin_id)

  # Further split the bins into clusters to optimized the execution.
  bin_cluster <- split(bin_list, ceiling(seq_along(bin_list)/(length(bin_list)/(n_threads*2))))

  default_simScore <- BiocParallel::bplapply(bin_cluster, function(bin_cluster_i){
    cluster_results <- lapply(bin_cluster_i, measure_seq_similarity) %>% dplyr::bind_rows()

    return(cluster_results)
  }, BPPARAM = BiocParallel::MulticoreParam(workers = n_threads, progressbar = T, stop.on.error = T)) %>%
    dplyr::bind_rows()

  default_simScore %>% saveRDS(default_simScore_path)
}else{
  default_simScore <- readRDS(default_simScore_path)
}

default_simScore <- default_simScore %>% dplyr::filter(n_comp != 0)

#----------------------------------------------------------------------------- #
## 2.2 Plot the results ----
fig_subtitle <- paste0(substring(utr_site, 1, 1), "' UTR Truncation")
plot_sim_results(default_simScore, subtitle = fig_subtitle, fig_output = plot_default_simScore_path)

############################################################################## #
# ---- 3. Random Label Simulation ----

#----------------------------------------------------------------------------- #
## 3.1 Randomize transcript table ----
seed <- 1
set.seed(seed)

# Randomize the mapping between ORF IDs and ORF Sequences
random_orf_tbl <- orf_tbl %>% dplyr::mutate(orf_seq = sample(orf_seq))

# Update the transcript table with the new random ORF Sequences associated to
# the pre-existing ORF ID. This mapping simulates a randomization of the ORF
# Labelling, but ensures that two transcripts that previously mapped to the
# same sequence still do.
random_tx_tbl <- tx_tbl %>%
  dplyr::select(-orf_seq) %>%
  dplyr::left_join(random_orf_tbl, by = c("orf_id", "gene_id"))

#----------------------------------------------------------------------------- #
## 3.2 Run the Similarity Scoring ----
if(!file_test("-f", simScore_path(file_prefix = "randomLabel", seed))){
  # Validate input table
  validate_df(random_tx_tbl, c("bin_id", "orf_id", "orf_seq"))

  # Split the transcript table into bin groups for parallel processing
  bin_list <- random_tx_tbl %>%
    dplyr::select(bin_id, orf_id, orf_seq) %>%
    split(., .$bin_id)

  # Further split the bins into clusters to optimized the execution.
  bin_cluster <- split(bin_list, ceiling(seq_along(bin_list)/(length(bin_list)/(n_threads*2))))

  rl_simScore <- BiocParallel::bplapply(bin_cluster, function(bin_cluster_i){
    cluster_results <- lapply(bin_cluster_i, measure_seq_similarity) %>% dplyr::bind_rows()

    return(cluster_results)
  }, BPPARAM = BiocParallel::MulticoreParam(workers = n_threads, progressbar = T, stop.on.error = T)) %>%
    dplyr::bind_rows()

  rl_simScore %>% saveRDS(simScore_path(file_prefix = "randomLabel", seed))
}else{
  rl_simScore <- readRDS(simScore_path(file_prefix = "randomLabel", seed))
}

rl_simScore <- rl_simScore %>% dplyr::filter(n_comp != 0)

#----------------------------------------------------------------------------- #
## 3.3 Plot the results ----
fig_subtitle <- paste0(substring(utr_site, 1, 1), "' UTR Truncation")
plot_sim_results(rl_simScore, subtitle = fig_subtitle, fig_output = plot_rl_simScore_path)


############################################################################## #
# ---- 4. Automate Randomization ----

#----------------------------------------------------------------------------- #
## 4.1 Randomization functions ----

#' Randomize ORF Sequences by labels
#'
#' This randomization ensures that Transcripts with the same ORF sequence prior
#' to randomization end up with the same random sequence.
random_label <- function(tx_tbl, orf_tbl, seed = 0){
  # Set a seed for reproducibility
  set.seed(seed)

  # Randomize the mapping between ORF IDs and ORF Sequences
  random_orf_tbl <- orf_tbl %>% dplyr::mutate(orf_seq = sample(orf_seq))

  # Update transcript table with new random ORF Sequences.
  random_tx_tbl <- tx_tbl %>%
    dplyr::select(-orf_seq) %>%
    dplyr::left_join(random_orf_tbl, by = c("orf_id", "gene_id"))

  return(random_tx_tbl)
}

#' Randomize ORF Sequences by labels and genes
#'
#' This randomization ensures that Transcripts with the same ORF sequence prior
#' to randomization end up with the same random sequence. Randomization only
#' happens within genes.
random_label_gene <- function(tx_tbl, orf_tbl, seed = 0){
  # Set a seed for reproducibility
  set.seed(seed)

  # Randomize the mapping between ORF IDs and ORF Sequences
  random_orf_tbl <- orf_tbl %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(orf_seq = sample(orf_seq)) %>%
    dplyr::ungroup()

  # Update transcript table with new random ORF Sequences
  random_tx_tbl <- tx_tbl %>%
    dplyr::select(-orf_seq) %>%
    dplyr::left_join(random_orf_tbl, by = c("orf_id", "gene_id"))

  return(random_tx_tbl)
}


#' Randomize ORF Sequences
random_full <- function(tx_tbl, orf_tbl, seed = 0){
  # Set a seed for reproducibility
  set.seed(seed)

  # Randomize the mapping between Transcripts and the ORF Sequences
  random_tx_tbl <- tx_tbl %>% dplyr::mutate(orf_seq = sample(orf_seq))

  return(random_tx_tbl)
}

#' Randomize ORF Sequences by genes
random_full_gene <- function(tx_tbl, orf_tbl, seed = 0){
  # Set a seed for reproducibility
  set.seed(seed)

  # Randomize the mapping between Transcripts and the ORF Sequences within the
  # same gene
  random_tx_tbl <- tx_tbl %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(orf_seq = sample(orf_seq)) %>%
    dplyr::ungroup()

  return(random_tx_tbl)
}

#' Return the input table without randomization
not_random <- function(tx_tbl, orf_tbl, seed = 0) return(tx_tbl)

#----------------------------------------------------------------------------- #
## 4.2 Pipeline functions ----
compute_bin_similarity <- function(transcript_table, n_threads = 64, output_path = ""){
  # If output already exists, return the results
  if(file_test("-f", output_path)) return(readRDS(output_path))

  # Validate input table
  validate_df(transcript_table, c("bin_id", "orf_id", "orf_seq"))

  # Split the transcript table into bin groups for parallel processing
  bin_list <- transcript_table %>%
    dplyr::select(bin_id, orf_id, orf_seq) %>%
    split(., .$bin_id)

  # Further split the bins into clusters to optimized the execution.
  bin_cluster <- split(bin_list, ceiling(seq_along(bin_list)/(length(bin_list)/(n_threads*2))))

  simScore <- BiocParallel::bplapply(bin_cluster, function(bin_cluster_i){
    cluster_results <- lapply(bin_cluster_i, measure_seq_similarity) %>% dplyr::bind_rows()

    return(cluster_results)
  }, BPPARAM = BiocParallel::MulticoreParam(workers = n_threads, progressbar = T, stop.on.error = T)) %>%
    dplyr::bind_rows()

  # Save the results
  if(output_path != "") simScore %>% saveRDS(output_path)
  return(simScore)
}

run_permutations <- function(tx_tbl, orf_tbl, rand_function, file_prefix, run_seeds = 1){
  # Loop through all the seeds provided
  perm_results <- lapply(run_seeds, function(seed){
    iter_path <- simScore_path(file_prefix, seed)
    if(file_test("-f", iter_path)){
      # If output already exists, avoid running the randomization step
      perm_score <- readRDS(iter_path)
    }else{
      cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]  Starting seed ", seed), sep = "\n")

      # Randomize and compute the similarity for this seed
      random_tx_tbl <- rand_function(tx_tbl, orf_tbl, seed)
      perm_score <- compute_bin_similarity(random_tx_tbl, output_path = iter_path)
    }

    # Return the similarity for these permutations
    return(perm_score)
  }) %>% dplyr::bind_rows(.id = "perm") %>% dplyr::mutate(perm = as.integer(perm))

  # Return the permutation results
  return(perm_results)
}


#----------------------------------------------------------------------------- #
## 4.3 Run simulations ----

def_simScore <- run_permutations(tx_tbl, orf_tbl, rand_function = not_random, file_prefix = "notRandom", run_seeds = 1)
rl_simScore <- run_permutations(tx_tbl, orf_tbl, rand_function = random_label, file_prefix = "randomLabel", run_seeds = 1:50)
rlg_simScore <- run_permutations(tx_tbl, orf_tbl, rand_function = random_label_gene, file_prefix = "randomLabelGene", run_seeds = 1:50)

rf_simScore <- run_permutations(tx_tbl, orf_tbl, rand_function = random_full, file_prefix = "randomFull", run_seeds = 1:50)
rfg_simScore <- run_permutations(tx_tbl, orf_tbl, rand_function = random_full_gene, file_prefix = "randomFullGene", run_seeds = 1:50)

#----------------------------------------------------------------------------- #
## 4.4 Extract results and compare ----
def_simScore <- def_simScore %>% dplyr::filter(n_comp != 0)
rl_simScore <- rl_simScore %>% dplyr::filter(n_comp != 0)
rlg_simScore <- rlg_simScore %>% dplyr::filter(n_comp != 0)
rf_simScore <- rf_simScore %>% dplyr::filter(n_comp != 0)
rfg_simScore <- rfg_simScore %>% dplyr::filter(n_comp != 0)


default_df <- def_simScore
random_df <- rl_simScore
compare_model_means <- function(default_df, random_df) {
  obs_df <- default_df %>% dplyr::select(bin_id, sim_pid, n_comp)

  # Merge all permutation results
  random_sim_list <- random_df %>%
    dplyr::select(perm, bin_id, sim_pid, n_comp) %>%
    dplyr::group_by(n_comp, perm) %>%
    dplyr::summarise(random_mean = mean(sim_pid)) %>%
    dplyr::summarise(random_sim_list = list(random_mean),
                     random_mean = mean(random_mean))

  # Combine the two data.frame and calculate the empirical p-value and the
  # binomial p-value.
  model_df <- obs_df %>%
    dplyr::group_by(n_comp) %>%
    dplyr::summarise(obs_mean = mean(sim_pid)) %>%
    dplyr::left_join(random_sim_list, by = "n_comp") %>%
    dplyr::group_by(n_comp) %>%
    dplyr::mutate(p_empirical = (1 + sum(random_sim_list[[1]] >= obs_mean)) / (length(random_sim_list[[1]]) + 1),
                  p_binomial = 1 - pnorm(obs_mean, mean = random_mean, sd = sd(random_sim_list[[1]]))) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(n_comp, obs_mean, random_mean, p_empirical, p_binomial)

  return(model_df)
}

rl_pval <- compare_model_means(def_simScore, rl_simScore)
rlg_pval <- compare_model_means(def_simScore, rlg_simScore)
rf_pval <- compare_model_means(def_simScore, rf_simScore)
rfg_pval <- compare_model_means(def_simScore, rfg_simScore)

############################################################################## #
# ---- X. Additional Analyses ----

#----------------------------------------------------------------------------- #
## X.1 Random comparison between ORF ----
r <- BiocParallel::bplapply(1:32, function(i){
  r <- sapply(1:100, function(x){
    seqs <- sample(orf_tbl$orf_seq, size = 2, replace = T)
    pwalign::pairwiseAlignment(seqs[1], seqs[2], gapOpening = 11, gapExtension = 1, substitutionMatrix = "BLOSUM62") %>% pwalign::pid()
  })
}, BPPARAM = BiocParallel::MulticoreParam(workers = 32, progressbar = T))

cat("Average mean of random ORF comparisons:", mean(unlist(r))/100)


#----------------------------------------------------------------------------- #
## X.2 Results of Random Label Genes by genes with many ORF (not implemented yet) ----

# valid_genes <- orf_tbl %>% dplyr::group_by(gene_id) %>% dplyr::count() %>% dplyr::filter(n > 5) %>% dplyr::pull(gene_id)
# rlg_simScore %>%
#   dplyr::left_join(bin_tbl, by = "bin_id") %>%
#   dplyr::filter(gene_id %in% valid_genes) %>%
#   plotSimResults()
# plotSimResults(rlg_simScore)
