## _________________________________________________
##
## scUTRquant pseudoalignment rates
##
## Aim: Gather the kallisto pseudoalignment rates of every ENDome of one
##      {dataset}.{group}.{merge_method}, tabulate them and plot them per ENDome
##      ({orf_filter}.w{width}.{txEnd}).
##
## Project: ENDome generation - UTR Quantification
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-01
##
## Latest Version: v1.3 (2026-10-16)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - Input: one {endome_name}.run_info.tsv per ENDome (rule collect_run_info),
##      one row per scUTRquant library. {orf_filter}, {width} and {txEnd} are
##      read from the end of the file name.
##    - params$regions (optional): the plots are split by the first of these
##      strings found in each library name ("unknown" if none).
##    - Boxplots are only drawn with >= params$min_samples_boxplot libraries.
##
## Changelog:
##    - v1.0 (2026-10-01): Initial version.
##    - v1.1 (2026-10-16): One config per ENDome (input$config_files).
##    - v1.2 (2026-10-16): Reads the per-ENDome run_info TSVs instead of the
##                         scUTRquant configs (test_data compatible).
##    - v1.3 (2026-10-16): region column (params$regions) and one plot per metric
##                         (output$pseudoaligned, output$unique).
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
  test_dir <- "data/test_data/Wood.control.st_ref"

  snakemake <- Snakemake(
    input = list(
      run_info = list.files(test_dir, pattern = "^test\\..*\\.run_info\\.tsv$", full.names = TRUE)
    ),
    output = list(
      tsv = file.path(test_dir, "interactive.pseudoalignment.tsv"),
      p_pseudoaligned = file.path(test_dir, "interactive.pseudoalignment_rate.png"),
      p_unique = file.path(test_dir, "interactive.pseudoalignment_unique.png")
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

min_samples_boxplot <- params$min_samples_boxplot
prefix <- paste(wc$dataset, wc$group, wc$merge_method, sep = ".")
regions <- unlist(params$regions)              # null in Snakemake -> NULL here
draw_regions <- length(regions) > 0

#### Truncation-site colours (first two of ggsci's Simpsons palette, as in T07)
txEnd_colours <- c("3p" = "#FED439", "5p" = "#709AE1")
txEnd_labels <- c("3p" = "3' ENDome", "5p" = "5' ENDome")
metric_labels <- c("p_pseudoaligned" = "Pseudoaligned [%]", "p_unique" = "Uniquely pseudoaligned [%]")
orf_filter_labels <- c("ref_pc" = "Reference protein-coding", "pc_novel" = "Coding, incl. novel isoforms",
                       "pc_all" = "All coding", "all" = "Unfiltered")

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Order truncation widths numerically, with non-numeric widths (e.g. "UTR") last
width_levels <- function(width) {
  width <- unique(width)
  num <- suppressWarnings(as.numeric(width))
  c(width[!is.na(num)][order(num[!is.na(num)])], sort(width[is.na(num)]))
}

#' Box plot of one pseudoalignment metric per ENDome
#'
#' @param plot_data The pseudoalignment table, with width as an ordered factor
#' @param metric Column to plot ("p_pseudoaligned" or "p_unique")
#' @param title Plot title
#' @return A ggplot (regions as facet rows when draw_regions is TRUE)
plot_metric <- function(plot_data, metric, title) {
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
# ---- 2. Pseudoalignment plot ----
n_samples <- dplyr::n_distinct(pseudo_table$sample_id)

plot_data <- pseudo_table %>%
  dplyr::mutate(width = factor(width, levels = width_levels(width)))

pseudoaligned_plot <- plot_metric(plot_data, "p_pseudoaligned", "Kallisto Bus - Pseudoalignment rate")
unique_plot <- plot_metric(plot_data, "p_unique", "Kallisto Bus - Unique pseudoalignment rate")

plot_height <- if (draw_regions) 3 + 2.5 * dplyr::n_distinct(plot_data$region) else 5.5
ggsave(smk_outputs$p_pseudoaligned, plot = pseudoaligned_plot, width = 11, height = plot_height, units = "in", dpi = 300)
ggsave(smk_outputs$p_unique, plot = unique_plot, width = 11, height = plot_height, units = "in", dpi = 300)