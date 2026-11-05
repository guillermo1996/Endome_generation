#!/usr/bin/env Rscript

## _________________________________________________
##
## scUTRquant pseudoalignment rates
##
## Aim: Gather the kallisto pseudoalignment rates of every ENDome of one
## {dataset}.{group}.{merge_method}, tabulate them and plot them per ENDome
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S06a-Pseudoalignment_Rates
## Date Created: 2026-10-01
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.4 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - {orf_filter}, {width} and {txEnd} are read from the end of the input
##      file names.
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.4 (2026-10-28): Standardization of the script.
##    - v1.3 (2026-10-16): Region column (params$regions) and one plot per
##      metric.
##    - v1.2 (2026-10-16): Reads the per-ENDome run_info tables.
##    - v1.1 (2026-10-16): One scUTRquant config per ENDome.
##    - v1.0 (2026-10-01): Initial version.
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

  rs <- "results/Wood"
  prefix <- "Wood.control.st_ref"
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      run_info = list.files(file.path(rs, "06-UTR_Quantification/pseudoalignment"),
                            pattern = paste0("^", prefix, "\\..*\\.run_info\\.tsv$"), full.names = TRUE)
    ),
    output = list(
      tsv = file.path(rs, "temp", paste0(prefix, ".pseudoalignment.tsv")),
      p_pseudoaligned = file.path(rs, "temp", paste0(prefix, ".pseudoalignment.png")),
      p_unique = file.path(rs, "temp", paste0(prefix, ".unique.png"))
    ),
    params = list(
      min_samples_boxplot = 3,
      dataset_name = "hardy_snRNAseq",
      regions = c("ACG", "IPL")
    ),
    wildcards = list(
      dataset = "Wood",
      group = "control",
      merge_method = "st_ref"
    ),
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

min_samples_boxplot <- params$min_samples_boxplot
prefix <- paste(wc$dataset, wc$group, wc$merge_method, sep = ".")
regions <- unlist(params$regions)
draw_regions <- length(regions) > 0

### Plot colours and labels
txEnd_colours <- c("3p" = "#FED439", "5p" = "#709AE1")
txEnd_labels <- c("3p" = "3' ENDome", "5p" = "5' ENDome")
metric_labels <- c("p_pseudoaligned" = "Pseudoaligned [%]", "p_unique" = "Uniquely pseudoaligned [%]")
orf_filter_labels <- c("ref_pc" = "Reference protein-coding", "pc_novel" = "Coding, incl. novel isoforms",
                       "pc_all" = "All coding", "all" = "Unfiltered")

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Order the truncation widths numerically, with non-numeric widths (e.g. "UTR") last
widthLevels <- function(width) {
  width <- unique(width)
  num <- suppressWarnings(as.numeric(width))
  c(width[!is.na(num)][order(num[!is.na(num)])], sort(width[is.na(num)]))
}

#' Box plot of one pseudoalignment metric per ENDome
#'
#' @param plot_data Pseudoalignment table, with width as an ordered factor.
#' @param metric Column to plot ("p_pseudoaligned" or "p_unique").
#' @param title Plot title.
#' @param n_samples Number of scUTRquant libraries.
#'
#' @returns A ggplot, with the regions as facet rows when `draw_regions` is TRUE.
plotMetric <- function(plot_data, metric, title, n_samples) {
  group_cols <- c("orf_filter", "width", "txEnd", if (draw_regions) "region")
  median_labels <- plot_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(rate = median(.data[[metric]]), label_y = max(.data[[metric]]), .groups = "drop")

  plot_data %>%
    ggplot(aes(x = width, y = .data[[metric]], fill = txEnd)) +
    {if (n_samples >= min_samples_boxplot) geom_boxplot(width = 0.5, position = position_dodge(width = 0.55), outlier.shape = NA)} +
    geom_text(aes(y = label_y, label = paste0(round(rate, 1), "%"), group = txEnd), data = median_labels,
      position = position_dodge(width = 0.55), vjust = -0.6, fontface = "bold", size = 3) +
    scale_fill_manual(name = "Truncation site", values = txEnd_colours, labels = txEnd_labels) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
    scale_x_discrete(labels = function(x) ifelse(grepl("^\\d+$", x), paste0(x, " nt"), x)) +
    facet_grid(rows = if (draw_regions) vars(region) else NULL, cols = vars(orf_filter), scales = "free_y",
      labeller = labeller(orf_filter = orf_filter_labels)) +
    labs(x = "Truncation width [bp]", y = metric_labels[[metric]],
         title = title,
         subtitle = paste0("ENDome: ", prefix, "\n",
                           "Evaluation dataset: ", params$dataset_name, " (", n_samples, " sample", ifelse(n_samples == 1, "", "s"), ")",
                           ifelse(n_samples >= min_samples_boxplot, "", " - one point per sample, labels are medians"))) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.key.size = unit(1.5, "line"),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      strip.text.y = element_text(angle = 0)
    )
}

############################################################################## #
# ---- 1. Pseudoalignment table ----
run_info_files <- unlist(smk_inputs$run_info)
names(run_info_files) <- sub("\\.run_info\\.tsv$", "", basename(run_info_files))

pseudo_table <- run_info_files %>%
  purrr::map(readr::read_tsv) %>%
  dplyr::bind_rows(.id = "endome") %>%
  tidyr::extract(endome, into = c("orf_filter", "width", "txEnd"), regex = "\\.([^.]+)\\.w([^.]+)\\.(\\dp)$", remove = FALSE) %>%
  dplyr::mutate(region = if (draw_regions) stringr::str_extract(sample_id, paste(regions, collapse = "|")) else NA_character_,
                region = tidyr::replace_na(region, "unknown"), .after = sample_id) %>%
  dplyr::mutate(prefix = prefix, dataset_name = params$dataset_name, .before = 1)

if (any(is.na(pseudo_table$txEnd))) {
  stop("Could not parse {orf_filter}.w{width}.{txEnd} from: ",
       paste(unique(pseudo_table$endome[is.na(pseudo_table$txEnd)]), collapse = ", "), call. = FALSE)
}

readr::write_tsv(pseudo_table, smk_outputs$tsv)

############################################################################## #
# ---- 2. Pseudoalignment plots ----
n_samples <- dplyr::n_distinct(pseudo_table$sample_id)

plot_data <- pseudo_table %>%
  dplyr::mutate(width = factor(width, levels = widthLevels(width)))

pseudoaligned_plot <- plotMetric(plot_data, "p_pseudoaligned", "Kallisto Bus - Pseudoalignment rate", n_samples)
unique_plot <- plotMetric(plot_data, "p_unique", "Kallisto Bus - Unique pseudoalignment rate", n_samples)

plot_height <- if (draw_regions) 3 + 2.5 * dplyr::n_distinct(plot_data$region) else 5.5
ggsave(smk_outputs$p_pseudoaligned, plot = pseudoaligned_plot, width = 11, height = plot_height, units = "in", dpi = 300)
ggsave(smk_outputs$p_unique, plot = unique_plot, width = 11, height = plot_height, units = "in", dpi = 300)
