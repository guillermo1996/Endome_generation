## _________________________________________________
##
## Comparisons of Similarity Scores
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
## Latest version: v1.0 (2025-05-29)
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
utr_site <- "5p"
n_threads <- 64

#----------------------------------------------------------------------------- #
## 0.3 Script Paths ----
main_path <- "~/RytenLab-Research/38-Endome_generation" # here::here()
script_path <- file.path(main_path, "R")

### Input Paths
simScore_3p_path <- file.path(script_path, "Results/SimScore_Simulations_3p/notRandom_SimScore.rds")
simScore_5p_path <- file.path(script_path, "Results/SimScore_Simulations_5p/notRandom_SimScore.rds")

### Output Paths
plot_default_simScore_path <- file.path(script_path, "figures/SimScore_Simulations", paste0("simScore_comparison.png"))
plot_binSize_3p_path <- file.path(script_path, "figures/SimScore_Simulations/binSize_3p.png")
plot_binSize_5p_path <- file.path(script_path, "figures/SimScore_Simulations/binSize_5p.png")
plot_binSize_comparison_path <- file.path(script_path, "figures/SimScore_Simulations/binSize_comparison.png")

#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----

#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----
source(file.path(script_path, "R/hf_graph_and_themes.R"))
source(file.path(script_path, "R/Phf-Endome/hf_sim_score.R"))

############################################################################## #
# ---- 1. Section 1 of the script ----

#----------------------------------------------------------------------------- #
## 1.1 Load previous results ----
simScore_3p <- readRDS(simScore_3p_path) %>% dplyr::mutate(UTR = factor("3'", levels = c("5'", "3'")))
simScore_5p <- readRDS(simScore_5p_path) %>% dplyr::mutate(UTR = factor("5'", levels = c("5'", "3'")))

############################################################################## #
# ---- 2. Similarity Distribution ----
plot_sim_comparisons(simScore_3p, simScore_5p, fig_output = plot_default_simScore_path)


############################################################################## #
# ---- 3. Number of ORFs Distribution ----
plot_bin_size(simScore_3p, subtitle = "3' UTR Trucation", fig_output = plot_binSize_3p_path)
plot_bin_size(simScore_5p, subtitle = "5' UTR Trunction", fig_output = plot_binSize_5p_path)

plot_bin_size_comparison(simScore_3p, simScore_5p, fig_output = plot_binSize_comparison_path)
