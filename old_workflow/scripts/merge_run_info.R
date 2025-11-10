#!/usr/bin/env Rscript

################################################################################
## Mock `snakemake` preamble
################################################################################

if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots=c(
      input='list',
      output='list',
      params='list',
      wildcards='list',
      log="list",
      threads='numeric'
    )
  )
  snakemake <- Snakemake(
    input=list(input_dir = "~/RytenLab-Research/38-Endome_generation/results_k15/06-Evaluation/kallisto/sq3.annotated.w500.3p"),
    output=list(png = "~/RytenLab-Research/38-Endome_generation/results_k15/figures/sq3.annotated.w500.3p.pseudoaligned_rate.png")
  )
}

################################################################################
## Libraries and Parameters
################################################################################

library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)

# snakemake@input$input_files <- c("results_k15/06-Evaluation/kallisto/sq3.annotated.w500.3p/borah_IPL_S3/run_info.json",
#                                  "results_k15/06-Evaluation/kallisto/sq3.annotated.w500.3p/bovon_ACG_S1/run_info.json")
################################################################################
## Load Data and plot
################################################################################

raw_data <- snakemake@input$input_files %>%
  lapply(function(x) {
    if(basename(x) == "run_info.json"){
      jsonlite::fromJSON(x) %>%
        as.data.frame() %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(sample_id = basename(dirname(x))) %>%
        dplyr::relocate(sample_id)
    }
  }) %>%
  dplyr::bind_rows()

plot_data <- raw_data %>%
  tidyr::pivot_longer(c(p_pseudoaligned, p_unique), names_to = "label", values_to = "perc") %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(mean_perc = median(perc, na.rm = T)) %>%
  dplyr::ungroup()

plot <- plot_data %>%
  ggplot(aes(x = reorder(sample_id, perc), y = perc, fill = label)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) +
  geom_hline(aes(yintercept = mean_perc, color = label)) +
  geom_text(aes(x = 1, y = mean_perc*1.05, label = paste0("median = ", round(mean_perc, 2)))) +
  scale_y_continuous(n.breaks = 10) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=4))


################################################################################
## Export results to disk
################################################################################

dir.create(dirname(snakemake@output$png), showWarnings = F, recursive = T)
ggsave(snakemake@output$png, plot, width = 12, height = 9, dpi = 300, units = "in", bg = "white")
