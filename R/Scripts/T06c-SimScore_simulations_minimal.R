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
shhh(library(pwalign))
shhh(library(patchwork))
shhh(library(BiocParallel))
shhh(library(tidyverse))

## Package options - these are some examples and do not modify the script usage.
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)

#----------------------------------------------------------------------------- #
## 0.2 Script Paths ----
main_path <- "/home/grocamora/RytenLab-Research/38-Endome_generation" # here::here()
script_path <- file.path(main_path, "R")

### Input Paths
sq3_class_path<- file.path(main_path, "results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt")
table_mapping_path <- file.path(script_path, "Results/Table_mapping")

tx_tbl_path <- file.path(table_mapping_path, "tx_tbl.rds")
tr_tbl_path <- file.path(table_mapping_path, "tr_tbl.rds")
bin_tbl_path <- file.path(table_mapping_path, "bin_tbl.rds")
orf_tbl_path <- file.path(table_mapping_path, "orf_tbl.rds")

map_tx_tr_path <- file.path(table_mapping_path, "map_tx_tr.rds")
map_tr_bin_path <- file.path(table_mapping_path, "map_tr_bin.rds")

### Output Paths
def_simScore_path <- file.path(script_path, paste0("Results/SimScore_Simulations2/notRandom_SimScore_s1.rds"))
rl_simScore_path <- function(seed) file.path(script_path, paste0("Results/SimScore_Simulations2/randomizeLabels_SimScore_s", seed, ".rds"))
simScore_path <- function(file_prefix, seed) file.path(script_path, paste0("Results/SimScore_Simulations2/", file_prefix, "_SimScore_s", seed, ".rds"))

#### Create output directory
dir.create(dirname(def_simScore_path), showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.3 Script Parameters ----

#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----

#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----
source(file.path(script_path, "R/Project_Specific/hf_table_mapping.R"))

randomizeLabels <- function(tx_tbl, orf_tbl, seed = 0){
  # Set a seed for reproducibility
  set.seed(seed)

  # Randomize the mapping between ORF IDs and ORF Sequences
  random_orf_tbl <- orf_tbl %>%
    dplyr::mutate(orf_seq = sample(orf_seq))

  # Update the transcript table with the new random ORF Sequences associated to
  # the pre-existing ORF ID. This mapping simulates a randomization of the ORF
  # Labelling, but ensures that two transcripts that previously mapped to the
  # same sequence still do.
  random_tx_tbl <- tx_tbl %>%
    dplyr::select(-orf_seq) %>%
    dplyr::left_join(random_orf_tbl, by = c("orf_id", "gene_id"))

  return(random_tx_tbl)
}

randomizeLabelsGene <- function(tx_tbl, orf_tbl, seed = 0){
  # Set a seed for reproducibility
  set.seed(seed)

  # Randomize the mapping between ORF IDs and ORF Sequences
  random_orf_tbl <- orf_tbl %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(orf_seq = sample(orf_seq))

  # random_orf_tbl <- orf_tbl %>%
  #   dplyr::left_join(dplyr::select(tx_tbl, gene_id, orf_id), by = "orf_id") %>%
  #   dplyr::relocate(gene_id) %>%
  #   dplyr::distinct(gene_id, orf_id, .keep_all = T) %>%
  #   dplyr::arrange(gene_id) %>%
  #   dplyr::group_by(gene_id) %>%
  #   dplyr::mutate(orf_id = sample(orf_id)) %>%
  #   dplyr::ungroup()

  # Update the transcript table with the new random ORF Sequences associated to
  # the pre-existing ORF ID. This mapping simulates a randomization of the ORF
  # Labelling, but ensures that two transcripts that previously mapped to the
  # same sequence still do.
  random_tx_tbl <- tx_tbl %>%
    dplyr::select(-orf_seq) %>%
    dplyr::left_join(random_orf_tbl, by = c("orf_id", "gene_id"))

  return(random_tx_tbl)
}

notRandom <- function(tx_tbl, orf_tbl, seed = 0) return(tx_tbl)

#----------------------------------------------------------------------------- #
## 4.2 Pipeline functions ----
compute_bin_similarity <- function(transcripts, n_threads = 128*2, workers = 64, output_path = ""){
  if(file_test("-f", output_path)) return(readRDS(output_path))

  orf_groups <- transcripts %>%
    dplyr::select(bin_id, orf_id, orf_seq) %>%
    split(., .$bin_id)

  # Further split the groups into clusters to optimized the execution.
  n_threads = n_threads
  orf_core_groups <- split(orf_groups, ceiling(seq_along(orf_groups)/(length(orf_groups)/n_threads)))

  # Run the similarity scoring
  orf_sim_output <- BiocParallel::bplapply(orf_core_groups, function(orf_core_group){
    orf_iter_results <- lapply(orf_core_group, measureSeqSimilarity) %>% dplyr::bind_rows()

    return(orf_iter_results)
  }, BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = T)) %>% dplyr::bind_rows()

  # Save the results
  if(output_path != "") orf_sim_output %>% saveRDS(output_path)
  return(orf_sim_output)
}

runPermutations <- function(tx_tbl, orf_tbl, rand_function, file_prefix, run_seeds = 1){
  perm_results <- lapply(run_seeds, function(seed){
    random_tx_tbl <- rand_function(tx_tbl, orf_tbl, seed)
    perm_score <- compute_bin_similarity(random_tx_tbl, output_path = simScore_path(file_prefix, seed)) %>%
      dplyr::mutate(n_comp = n_orf - short_seqs - long_seqs) %>%
      dplyr::relocate(bin_id, n_orf, n_comp)
  }) %>% dplyr::bind_rows(.id = "perm") %>% dplyr::mutate(perm = as.integer(perm))

  return(perm_results)
}

############################################################################## #
# ---- 1. Section 1 of the script ----

#----------------------------------------------------------------------------- #
## 1.1 Load previous tables ----
tx_tbl <- readRDS(tx_tbl_path)
tr_tbl <- readRDS(tr_tbl_path)
bin_tbl <- readRDS(bin_tbl_path)
orf_tbl <- readRDS(orf_tbl_path)
# sq3_class <- vroom::vroom(sq3_class_path)

# map_tx_tr <- readRDS(map_tx_tr_path)
# map_tr_bin_path <- readRDS(map_tr_bin_path)


#----------------------------------------------------------------------------- #
## 1.2 Load previous tables ----

def_simScore <- runPermutations(tx_tbl, orf_tbl, rand_function = notRandom, file_prefix = "notRandom", run_seeds = 1)
rl_simScore <- runPermutations(tx_tbl, orf_tbl, rand_function = randomizeLabels, file_prefix = "randomizeLabels", run_seeds = 1:50)
rlg_simScore <- runPermutations(tx_tbl, orf_tbl, rand_function = randomizeLabelsGene, file_prefix = "randomizeLabelsGene", run_seeds = 1:50)


#
# compute_pvals(default_df, randomLabel_df)
#
# orf_sim_output <- randomLabel_df
#
# plotSimResults <- function(orf_sim_output){
#   library(ggtext)
#   counts_df <- orf_sim_output %>%
#     dplyr::mutate(nTx = ifelse(n_orf <= 15, n_orf, 16)) %>%
#     dplyr::mutate(nTx = factor(nTx)) %>%
#     dplyr::mutate(nTx = fct_recode(nTx, "\u200B>15" = "16")) %>%
#     dplyr::count(nTx, name = "n") %>%
#     dplyr::mutate(label = paste0(nTx, "<br><span style = 'font-size:7pt'>n: ", format(n, big.mark = ""), "</span>"))
#
#   axis_labels <- setNames(counts_df$label, counts_df$nTx)
#
#   orf_sim_output %>%
#     dplyr::mutate(nTx = ifelse(n_orf <= 15, n_orf, 16)) %>%
#     dplyr::mutate(nTx = factor(nTx)) %>%
#     dplyr::mutate(nTx = fct_recode(nTx, "\u200B>15" = "16")) %>%
#     ggplot(aes(x = nTx, y = sim_null)) +
#     geom_violin(fill = "#8ab1cf", color = "black", alpha = 0.7, width = 1.05, adjust = 0.65) +
#     stat_summary(aes(fill = "Mean"), fun = mean, geom = "point", shape = 21, size = 2.5, color = "black", show.legend = T) +
#     scale_fill_manual(name = "", values = c("Mean" = "#f8766d")) +
#     scale_y_continuous(breaks = seq(0, 1, 0.25)) +
#     scale_x_discrete(labels = axis_labels) +
#     labs(x = "Number of Transcripts\n(Elements in bin)", y = "Similarity Score",
#          title = "Distribution of Similarity Scores by number of Transcripts") +
#     theme_minimal(base_size = 10) +
#     theme(
#       panel.grid.major.x = element_blank(),
#       plot.title = element_text(size = 14, face = "bold"),
#       axis.text.x = element_markdown(size = 9),
#       axis.text.y = ggplot2::element_text(color = "black", size = 9),
#       axis.title.x = ggplot2::element_text(face = "bold", size = 11, margin = margin(5, 0, 0, 0)),
#       axis.title.y = ggplot2::element_text(face = "bold", size = 11, margin = margin(0, 10, 0, 0)),
#     )
# }
