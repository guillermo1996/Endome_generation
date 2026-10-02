## _________________________________________________
##
## scUTRquant pseudoalignment rates
##
## Aim: Collect the kallisto `run_info.json` of every scUTRquant target and sample
##      for one {dataset}.{group}.{merge_method}, tabulate the pseudoalignment
##      rates and plot them per ENDome ({orf_filter}.w{width}.{txEnd}).
##
## Project: ENDome generation - UTR Quantification
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-01
##
## Latest Version: v1.0 (2026-10-01)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - The targets and samples are read from the resolved scUTRquant config
##      (`target`, `sample_file`), and the run_info.json files from scUTRquant's
##      own layout: <config dir>/data/kallisto/{target}/{sample_id}/run_info.json.
##      Only the targets of this {prefix} (= dataset.group.merge_method) are kept.
##    - One point per sample. Boxplots are drawn once there are enough samples to
##      summarise (params$min_samples_boxplot).
##    - Layout follows Deprecated/Scripts/T07-Kallisto_pseudoalignment_plots.R in
##      41-ENDome_generation_extra (x = width, fill = truncation site, facet = ORF
##      filter), with p_pseudoaligned and p_unique as rows.
##
## Changelog:
##    - v1.0 (2026-10-01): Initial version.
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
  test_dir <- "debug_results/Ebbert.control/06-UTR_Quantification-0ed3e"
  test_prefix <- "Ebbert.control.st_ref"
  snakemake <- Snakemake(
    input = list(
      config_file = file.path(test_dir, "scUTRquant/scUTRquant_config.yaml")
    ),
    output = list(
      tsv = file.path("claude_tmp", paste0(test_prefix, ".pseudoalignment.tsv")),
      png = file.path("claude_tmp", paste0(test_prefix, ".pseudoalignment.png"))
    ),
    params = list(
      min_samples_boxplot = 3
    ),
    wildcards = list(
      prefix = test_prefix
    ),
    log = list(),
    threads = 1,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages

shhh({
  library(conflicted)
  library(tidyverse)
  library(jsonlite)
  library(yaml)
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

if (length(snakemake@log) > 0) {
  log_con <- file(snakemake@log[[1]], open = "wt"); sink(log_con); sink(log_con, type = "message")
}

min_samples_boxplot <- params$min_samples_boxplot

#### Truncation-site colours (first two of ggsci's Simpsons palette, as in T07)
txEnd_colours <- c("3p" = "#FED439", "5p" = "#709AE1")
txEnd_labels <- c("3p" = "3' ENDome", "5p" = "5' ENDome")
metric_labels <- c("p_pseudoaligned" = "Pseudoaligned [%]", "p_unique" = "Uniquely pseudoaligned [%]")

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Keep the targets of one prefix and split them into {orf_filter}, {width} and {txEnd}
#'
#' @param target scUTRquant target names, "{prefix}.{orf_filter}.w{width}.{txEnd}"
#' @param prefix The {dataset}.{group}.{merge_method} prefix
#' @return A tibble with target, orf_filter, width and txEnd
parse_targets <- function(target, prefix) {
  tibble::tibble(target = target) %>%
    dplyr::filter(startsWith(target, paste0(prefix, "."))) %>%
    dplyr::mutate(suffix = substring(target, nchar(prefix) + 2)) %>%
    tidyr::extract(suffix, into = c("orf_filter", "width", "txEnd"), regex = "^(.+)\\.w([^.]+)\\.(\\dp)$")
}

#' Order truncation widths numerically, with non-numeric widths (e.g. "UTR") last
width_levels <- function(width) {
  width <- unique(width)
  num <- suppressWarnings(as.numeric(width))
  c(width[!is.na(num)][order(num[!is.na(num)])], sort(width[is.na(num)]))
}

############################################################################## #
# ---- 1. Locate run_info.json files ----

#----------------------------------------------------------------------------- #
## 1.1 Targets and samples from the scUTRquant config ----
scutrquant_config <- yaml::read_yaml(smk_inputs$config_file)
scutrquant_dir <- dirname(smk_inputs$config_file)

targets <- parse_targets(unlist(scutrquant_config$target), wc$prefix) %>%
  dplyr::filter(!is.na(txEnd))

if (nrow(targets) == 0) stop("No scUTRquant targets found for prefix '", wc$prefix, "'.", call. = FALSE)

samples <- readr::read_csv(scutrquant_config$sample_file) %>% dplyr::pull(sample_id)

#----------------------------------------------------------------------------- #
## 1.2 Expected run_info.json paths ----
run_info_files <- tidyr::expand_grid(targets, sample_id = samples) %>%
  dplyr::mutate(run_info = file.path(scutrquant_dir, "data/kallisto", target, sample_id, "run_info.json"))

missing_files <- run_info_files %>% dplyr::filter(!file.exists(run_info))
if (nrow(missing_files) > 0) {
  stop("Missing run_info.json files:\n  ", paste(missing_files$run_info, collapse = "\n  "), call. = FALSE)
}

############################################################################## #
# ---- 2. Pseudoalignment table ----
pseudo_table <- run_info_files %>%
  dplyr::mutate(info = purrr::map(run_info, jsonlite::fromJSON)) %>%
  tidyr::hoist(info, "n_processed", "n_pseudoaligned", "n_unique", "p_pseudoaligned", "p_unique", "kallisto_version") %>%
  dplyr::select(-info, -run_info) %>%
  dplyr::mutate(prefix = wc$prefix, dataset_name = scutrquant_config$dataset_name, .before = 1)

readr::write_tsv(pseudo_table, smk_outputs$tsv)

############################################################################## #
# ---- 3. Pseudoalignment plot ----
n_samples <- length(samples)

plot_data <- pseudo_table %>%
  tidyr::pivot_longer(c(p_pseudoaligned, p_unique), names_to = "metric", values_to = "rate") %>%
  dplyr::mutate(width = factor(width, levels = width_levels(width)),
                metric = factor(metric, levels = names(metric_labels)))

median_labels <- plot_data %>%
  dplyr::group_by(metric, orf_filter, width, txEnd) %>%
  dplyr::summarise(rate = median(rate), .groups = "drop")

dodge <- position_dodge(width = 0.75)

pseudo_plot <- plot_data %>%
  ggplot(aes(x = width, y = rate, fill = txEnd)) +
  {if (n_samples >= min_samples_boxplot) geom_boxplot(width = 0.6, position = dodge, outlier.shape = NA)} +
  geom_point(aes(group = txEnd), shape = 21, size = 2.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75, seed = 17)) +
  geom_text(aes(label = paste0(round(rate, 1), "%"), group = txEnd), data = median_labels,
            position = dodge, vjust = -1.2, fontface = "bold", size = 3.5) +
  scale_fill_manual(name = "Truncation site", values = txEnd_colours, labels = txEnd_labels) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  scale_x_discrete(labels = function(x) ifelse(grepl("^\\d+$", x), paste0(x, " nt"), x)) +
  facet_grid(rows = vars(metric), cols = vars(orf_filter), scales = "free_y",
             labeller = labeller(metric = metric_labels)) +
  labs(x = "Truncation width", y = "Pseudoalignment rate [%]",
       title = "Kallisto Bus - Pseudoalignment rate",
       subtitle = paste0("ENDome: ", wc$prefix, "\n",
                         "Evaluation dataset: ", scutrquant_config$dataset_name, " (", n_samples, " sample", ifelse(n_samples == 1, "", "s"), ")",
                         ifelse(n_samples >= min_samples_boxplot, "", " - one point per sample, labels are medians"))) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.key.size = unit(1.5, "line"),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(face = "bold")
  )

ggsave(smk_outputs$png, plot = pseudo_plot, width = 11, height = 7, units = "in", dpi = 300)

if (length(snakemake@log) > 0) {
  sink(type = "message"); sink(); close(log_con)
}
