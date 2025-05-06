## _________________________________________________
##
## Helper Functions: Graphs and Themes
##
## Aim: Include in one single file the functions and variables needed to bring a
## consistent theme to the graphs used in the analysis.
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 06/06/2023
##
## Copyright (c) Guillermo Rocamora Pérez, 2023
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.1 (18/01/2024)
## _________________________________________________
##
## - Notes:
##
## Many themes will be further modified per graph as required by the situation.
## This script includes only a guideline on the themes to use and other
## variables required to properly run the graphs.
##
## - Changelog:
##
##
## - Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________
##

################################################################################
# ---- 1. Required libraries ----
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(ggplot2)))

################################################################################
# ---- 2. Custom ggplot themes ----

## Default theme employed across all analyses
custom_gg_theme <- theme(
  plot.title = element_text(size = 14, face = "bold"),
  panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.text.x = ggplot2::element_text(color = "black", size = 9, angle = 0, hjust = 0.5),
  axis.text.y = ggplot2::element_text(color = "black", size = 9),
  axis.title.x = ggplot2::element_text(face = "bold", size = 11, margin = margin(5, 0, 0, 0)),
  axis.title.y = ggplot2::element_text(face = "bold", size = 11, margin = margin(0, 10, 0, 0)),
  panel.grid.minor = element_line(color = "#aaaaaa", linewidth = 0.05, linetype = 2),
  panel.grid.major.y = element_line(color = "#444444", linewidth = 0.05, linetype = 2),
  panel.grid.major.x = element_line(color = "#444444", linewidth = 0.05),
  panel.background = element_rect(fill = "#FBFBFB"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "top",
  legend.key = element_rect(color = "black", linewidth = 0.01),
  legend.key.size = unit(1, "lines"),
  strip.text.x = element_text(color = "black", face = "bold", size = 9),
  strip.background = element_rect(color = "black", linewidth = 1, linetype = "solid"),
  strip.text.y = element_text(color = "black", face = "bold", size = 9),
  plot.margin = margin(0.5, 0.5, 0.2, 0.5, "cm"),
  legend.box.just = "center",
  legend.margin = margin(0, 0, 0, 0, unit = "lines"),
  legend.spacing.y = unit(0.1, "lines")
)

## Theme modified to better represent titles and subtitles. It sets the legend
## on top right for more compress plots
custom_gg_theme_subtitle <- custom_gg_theme +
  theme(
    plot.subtitle = element_text(size = 10, vjust = 1, color = "black", margin = margin(b = -15, t = -5)),
    plot.title = element_text(size = 14, face = "bold", margin = margin(b = 9)),
    legend.justification = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.box.just = "center",
    legend.margin = margin(0, 0, -0.5, -0.5, unit = "lines")
  )

patchwork_annotation_theme_guide <- theme(
  plot.title = element_text(size = 18, face = "bold", margin = margin(b = 0.5,  unit = "lines")),
  plot.subtitle = element_text(size = 11, vjust = 1, color = "black", margin = margin(b = -1, t = 0,  unit = "lines")),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.justification = "right",
  legend.box.just = "center",
  legend.margin = margin(0, 0, 1, 0, unit = "lines"),
  plot.tag.position = c(0, 1),
  plot.tag = element_text(hjust = 0, vjust = -0.75)
)



################################################################################
# ---- 3. Project Specific variables & parameters ----

## 3.0 Mixed projects ----
var_labels <- c("prop_count" = "Proportion of Junction Counts",
                "prop_junc" = "Proportion of Unique Junctions")

## 3.1 PPMI Splicing ----
### Based on "jco" colors
PPMI_group_colors <- c(
  "All samples" = "#868686E5",
  "Control" = "#868686E5",
  "Case" = "#0073C2E5",
  "Healthy Control" = "#868686E5",
  "PD" = "#0073C2E5",
  "Genetic Cohort Unaffected" = "#CD534CE5",
  "Genetic Cohort PD" = "#EFC000E5",
  "Negative" = "#A73030FF",
  "Positive" = "#4A6990FF"
)

### Names for legends
PPMI_group_names <- c(
  "All samples" = "All samples",
  "case_control_other_latest" = "PD Status",
  "study_arm" = "Study arm"
)

### Labels for facets/axis
PPMI_group_labels <- c(
  "All samples" = "All samples",
  "Control" = "Control",
  "Case" = "Case",
  "Healthy Control" = "Healthy Control",
  "Genetic Cohort Unaffected" = "Genetic Cohort Unaffected",
  "Genetic Cohort PD" = "Genetic Cohort PD",
  "PD" = "PD"
)

## 3.2 ASAP Splicing ---
### Wood Splicing ----
Wood_group_colors <- c(
  "All samples" = "#868686E5",
  "Control" = "#868686E5",
  "PD" = "#0073C2E5"
)

Wood_group_names <- c(
  "All samples" = "All samples",
  # "case_control_other_latest" = "PD Status",
  "Group" = "PD Status"
)

Wood_group_labels <- c(
  "All samples" = "All samples",
  "Control" = "Control",
  "PD" = "Case"
)

Wood_tissue_labels <- c("C_CTX" = "Cingulate Cortex", "CAU" = "Caudate Nucleus",
                        "F_CTX" = "Frontal Cortex", "P_CTX" = "Parietal Cortex",
                        "PARA" = "Parahippocampal Gyrus", "PUT" = "Putamen",
                        "SN" = "Substantia Nigra", "T_CTX" = "Temporal Cortex")

### Hardy Splicing ----
Hardy_group_colors <- c(
  "All samples" = "#868686E5",
  "Control" = "#868686E5",
  "PD and PDD" = "#0073C2E5",
  "PD" = "#EFC000E5",
  "PDD" = "#CD534CE5",
  "Long" = "#EFC000E5",
  "Short" = "#CD534CE5",
  "Long Dementia" = "#7ba6da",
  "Short Dementia" = "#CD534CE5",
  "Long No Dementia" = "#0073C2E5",
  "Short No Dementia" = "#EFC000E5"
)

Hardy_group_names <- c(
  "All samples" = "All samples",
  "pd_status" = "PD Status",
  "GroupB" = "PD Status"
)

Hardy_group_labels <- c(
  "All samples" = "All samples",
  "Control" = "Control",
  "PD" = "PD",
  "PDD" = "PDD"
)

Hardy_tissue_labels <- c("ACG" = "Anterior Cingulate Gyrus",
                         "MFG" = "Middle Frontal Gyrus",
                         "IPL" = "Inferior Parietal Lobe",
                         "MTG" = "Middle Temporal Gyrus")

### Combination
ASAP_group_colors <- c(Wood_group_colors, Hardy_group_colors) %>% .[!duplicated(names(.))]
ASAP_group_names <- c(Wood_group_names, Hardy_group_names) %>% .[!duplicated(names(.))]
ASAP_group_labels <- c(Wood_group_labels, Hardy_group_labels) %>% .[!duplicated(names(.))]

ASAP_tissue_labels <- c(Wood_tissue_labels, Hardy_tissue_labels)
ASAP_tissue_labels2 <- lapply(seq_along(ASAP_tissue_labels), function(i){
  paste0(ASAP_tissue_labels[i], "\n(", names(ASAP_tissue_labels)[i], ")")
}) %>% unlist %>% setNames(names(ASAP_tissue_labels))
ASAP_tissue_labels3 <- lapply(seq_along(ASAP_tissue_labels), function(i){
  paste0(ASAP_tissue_labels[i], " (", names(ASAP_tissue_labels)[i], ")")
}) %>% unlist %>% setNames(names(ASAP_tissue_labels))

################################################################################
# ---- 4. Project Specific Functions ----

## 4.1 Splicing ----

#' Plot the junction categories
#'
#' @param annotated_prop_df dataframe, includes the proportions for each sample
#'   and category.
#' @param grouping_var string, field in the dataframe to group and compare the
#'   different samples.
#' @param prop_var string, variable to plot. Defaults to "prop_junc".
#' @param y_discontinuity numeric vector of two elements. The plot is condensed
#'   to remove the values between the two specified there. Recommended to be
#'   decimal numbers (e.g c(0.3, 0.6)), and it might not work as intended.
#'   Defaults to NULL.
#' @param custom_limits numeric vector of two elements, the custom limits to
#'   apply to the plot. Might not work as expected.
#' @param add_pval boolean, whether to add a Wilcoxon test between the different
#'   groups. Defaults to FALSE.
#' @param pval_y numeric, position to write the p-value of the previous test.
#'   Defaults to NULL.
#' @param ref_group string, name of the reference group to compare against in
#'   the Wilcoxon tests. Defaults to NULL.
#'
#' @return ggplot2 object.
#' @export
# plotJunctionCategories <- function(annotated_prop_df,
#                                    grouping_var,
#                                    x_var = "type",
#                                    prop_var = "prop_junc",
#                                    y_discontinuity = NULL,
#                                    custom_limits = c(0, 1),
#                                    add_pval = F,
#                                    pval_y = NULL,
#                                    ref_group = NULL) {
#   junc_prop_df <- annotated_prop_df %>% dplyr::rename(prop_var = !!sym(prop_var))
#
#   if (!is.null(ref_group)) junc_prop_df[[grouping_var]] <- junc_prop_df[[grouping_var]] %>% forcats::fct_relevel(., ref_group)
#
#   # Fix for null "grouping_var"
#   if (is.null(grouping_var)) {
#     junc_prop_df$`All samples` <- "All samples"
#     grouping_var <- "All samples"
#     add_pval <- F
#   }
#
#   # If a discontinuity in the Y-axis is introduced, we need to modify the values
#   # being plotted.
#   if (!is.null(y_discontinuity)) {
#     skip_len <- y_discontinuity[2] - y_discontinuity[1]
#     junc_prop_df <- junc_prop_df %>%
#       dplyr::mutate(prop_var = ifelse(prop_var < y_discontinuity[2],
#         prop_var + skip_len,
#         prop_var
#       ))
#   }
#
#   # Main plot
#   p <- junc_prop_df %>%
#     ggplot(aes(x = .data[[x_var]], y = prop_var, fill = .data[[grouping_var]])) +
#     geom_boxplot() +
#     scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)), limits = custom_limits, breaks = seq(0, 1, 0.2)) +
#     labs(x = "", y = "Proportion of junctions") +
#     scale_fill_manual(
#       values = grouping_colors,
#       guide = guide_legend(ncol = 2), name = groping_names[grouping_var],
#       labels = grouping_labels
#     )
#
#   # Modify the plot to include the Y-axis discontinuity.
#   # TODO: Fix custom limits to be compatible with the Y-discontinuity
#   if (!is.null(y_discontinuity)) {
#     p <- p +
#       scale_y_continuous(
#         limits = c(skip_len, custom_limits[2]),
#         breaks = seq(skip_len, custom_limits[2], 0.1),
#         labels = c(
#           seq(0, y_discontinuity[1], 0.1) %>% .[-length(.)],
#           paste0(y_discontinuity[2], "\n", y_discontinuity[1]),
#           seq(y_discontinuity[2], custom_limits[2], 0.1)[-1]
#         )
#       ) +
#       annotate(
#         geom = "segment", x = -Inf, xend = Inf, y = y_discontinuity[2] - 0.01, yend = y_discontinuity[2] - 0.01,
#         lty = "dashed", color = "#656565", lwd = 0.4
#       ) +
#       annotate(
#         geom = "segment", x = -Inf, xend = Inf, y = y_discontinuity[2] + 0.01, yend = y_discontinuity[2] + 0.01,
#         lty = "dashed", color = "#656565", lwd = 0.4
#       )
#   }
#
#   # Add P-values
#   if (add_pval) {
#     p <- p +
#       ggpubr::geom_pwc(aes(group = .data[[grouping_var]]),
#         label = "p.adj.signif",
#         p.adjust.method = "BH", p.adjust.by = "panel",
#         vjust = 0.6, label.size = 6, hide.ns = T, ref.group = ref_group, y.position = pval_y,
#         symnum.args = list(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns"))
#       )
#   }
#
#   return(p)
# }


#' Plot the junction categories
#'
#' @param prop_df dataframe, includes the proportions for each sample and
#'   category.
#' @param x_var string, name of the column to represent in the X-axis. Defaults
#'   to "type".
#' @param y_var string, name of the column to represent in the Y-axis. Defaults
#'   to "prop_junc".
#' @param grouping_var string, name of the column in the dataframe to group and
#'   compare the different samples. Defaults to NULL.
#' @param custom_limits numeric vector of two elements, list of lower and upper
#'   limits to the plots. It might not be applied if drawn boxes/significances
#'   are not within these limits. Defaults to NULL.
#' @param y_break numeric vector of two elements. The plot is condensed to
#'   remove the range specified in this argument. The input is expected to be as
#'   "c(0.3, 0.5)". Behaviour unstable, might cause issues when facetting.
#'   Defaults to NULL.
#' @param add_pval boolean, whether to add a Wilcoxon test between the different
#'   groups. Defaults to FALSE.
#' @param ref_group string, name of the reference group to compare against in
#'   the Wilcoxon tests. Defaults to NULL.
#' @param experiment string, global name of the experiment. Employed to color
#'   code the Groups and rename the legend. Defaults to "PPMI".
#'
#' @return ggplot2 object.
#' @export
plotJunctionCategories <- function(prop_df,
                                   x_var = "type",
                                   y_var = "prop_junc",
                                   grouping_var = NULL,
                                   custom_limits = NULL,
                                   y_break = NULL,
                                   add_pval = F,
                                   wilcox_groups = c("type"),
                                   stat_test = NULL,
                                   ref_group = NULL,
                                   experiment = "PPMI") {
  # Remove unused levels from "x_var" if is factor
  if (is.factor(prop_df[[x_var]])) prop_df <- prop_df %>% dplyr::mutate(!!x_var := droplevels(.data[[x_var]]))

  # Add a column name "split" to determine if the axis will be broken
  y_break_median <- ifelse(is.null(y_break), 0, median(y_break))
  prop_df <- prop_df %>% dplyr::mutate(split = case_when(!is.null(y_break) ~ .data[[y_var]] < y_break_median, .default = T))

  # If no grouping var is provided, a temporary column named "All samples" is
  # created to group all samples into it
  if (is.null(grouping_var)) {
    grouping_var <- "All samples"
    prop_df$`All samples` <- "All samples"
  }

  # Draw the main plot
  prop_plot <- prop_df %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_boxplot(aes(fill = .data[[grouping_var]])) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)), breaks = seq(0, 1, 0.1)) +
    custom_gg_theme

  # Modify the scale fill based on the experiment
  prop_plot <- modifyScaleColor(prop_plot, grouping_var, experiment)

  # If custom limits are applied, add a "scale_y_continuous" with modified
  # expand
  if (!is.null(custom_limits)) {
    prop_plot <- prop_plot +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), limits = custom_limits, breaks = seq(0, 1, 0.1))
  }

  # If a Y-axis break is added, we add a facet and draw invisible points to
  # force the plot to have the desired custom limits
  if (!is.null(y_break)) {
    # Modify the Y-break values to match with the later "expand" in
    # "scale_y_continuous"
    y_break <- y_break / c(1.07, 1)

    # Create a dataframe with points that will be drawn to adjust the graph
    # limits
    custom_limits_df <- data.frame(
      x = prop_df[[x_var]][1],
      y = c(y_break, custom_limits)
    ) %>%
      dplyr::mutate(split = y < median(y_break)) %>%
      dplyr::rename(!!y_var := y, !!x_var := x)

    # Add the invisible dots and use a facet to split the graph, breaking the
    # axis. It is neccesary to remove the "strip.text.y"
    prop_plot <- prop_plot +
      geom_point(data = custom_limits_df, shape = NA) +
      # facet_grid(.data[[y_var]] < median(y_break) ~ ., scales = "free_y", space = "free_y") +
      facet_grid(split ~ ., scales = "free_y", space = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = NULL, breaks = seq(0, 1, 0.1)) +
      theme(strip.text.y = element_blank())
  }

  # Add p-values of the Wilcoxon comparison between the groping_var
  if (add_pval) {
    if(is.null(stat_test)){
      # We do a manual test first that later will be added to the plot. It is
      # important to add X and Y positions of the brackets, and to correct the
      # Y-position to ignore non-significant elements. There is probably an option
      # to ignore them in the default functions, but I have not found an automatic
      # way to do so.
      stat_test <- prop_df %>%
        dplyr::group_by(split, across(all_of(wilcox_groups))) %>%
        rstatix::wilcox_test(as.formula(paste0(y_var, "~", grouping_var)), ref.group = ref_group) %>%
        rstatix::adjust_pvalue(method = "BH") %>%
        rstatix::add_significance(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns"))
    }else{
      stat_test <- stat_test %>%
        dplyr::filter(.y. == y_var)
    }

    # Fix Y-position
    stat_test <- stat_test %>%
      rstatix::add_xy_position(x = x_var, step.increase = 0.06, fun = "max") %>%
      fix_pval_y_position(x_var)

    # Add the manual stat test
    prop_plot <- prop_plot +
      ggpubr::stat_pvalue_manual(stat_test, label = "p.adj.signif", hide.ns = T, label.size = 6, vjust = 0.6, tip.length = 0.02)

    if(!is.null(custom_limits)){
      if(max(stat_test$y.position) > custom_limits[2]){
        prop_plot <- prop_plot +
          scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)), breaks = seq(0, 1, 0.1))
      }
    }
  }

  # Return the plot
  return(prop_plot)
}


#' Plot the junction categories with facet
#'
#' @param prop_df dataframe, includes the proportions for each sample and
#'   category.
#' @param x_var string, name of the column to represent in the X-axis. Defaults
#'   to "type".
#' @param y_var string, name of the column to represent in the Y-axis. Defaults
#'   to "prop_junc".
#' @param grouping_var string, name of the column in the dataframe to group and
#'   compare the different samples. Defaults to NULL.
#' @param split_var string, name of the column in the dataframe to split the
#'   graph in facets. Defaults to NULL.
#' @param add_pval boolean, whether to add a Wilcoxon test between the different
#'   groups. Defaults to FALSE.
#' @param ref_group string, name of the reference group to compare against in
#'   the Wilcoxon tests. Defaults to NULL.
#' @param n_row numeric, number of rows in the facet.
#' @param experiment string, global name of the experiment. Employed to color
#'   code the Groups and rename the legend. Defaults to "PPMI".
#'
#' @return ggplot2 object.
#' @export
plotJunctionCategoriesSplit <- function(prop_df,
                                        x_var = "type",
                                        y_var = "prop_junc",
                                        grouping_var = NULL,
                                        split_var = NULL,
                                        add_pval = F,
                                        ref_group = NULL,
                                        n_row = NULL,
                                        step.increase = 0.18,
                                        experiment = "PPMI") {
  # If no grouping var is provided, a temporary column named "All samples" is
  # created to group all samples into it
  if (is.null(grouping_var)) {
    grouping_var <- "All samples"
    prop_df$`All samples` <- "All samples"
  }

  # Draw the main plot with a facet split
  prop_plot <- prop_df %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_boxplot(aes(fill = .data[[grouping_var]])) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    facet_wrap(reformulate(split_var), scales = "free", nrow = n_row) +
    custom_gg_theme

  # Modify the scale fill based on the experiment
  prop_plot <- modifyScaleColor(prop_plot, grouping_var, experiment)

  # Add p-values of the Wilcoxon comparison between the groping_var
  if (add_pval) {
    # We do a manual test first that later will be added to the plot. It is
    # important to add X and Y positions of the brackets, and to correct the
    # Y-position to ignore non-significant elements. There is probably an option
    # to ignore them in the default functions, but I have not found an automatic
    # way to do so.
    stat_test <- prop_df %>%
      dplyr::group_by(type) %>%
      rstatix::wilcox_test(as.formula(paste0(y_var, "~", grouping_var)), ref.group = ref_group) %>%
      rstatix::adjust_pvalue(method = "BH") %>%
      rstatix::add_significance(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns")) %>%
      add_pval_y_position(x_var, y_var) %>%
      rstatix::add_x_position(x = x_var) %>%
      dplyr::mutate(xmin = xmin - x + 1, xmax = xmax - x + 1) %>%
      dplyr::filter(p.adj.signif != "ns")

    if(nrow(stat_test) > 0){
      stat_test <- stat_test %>%
        dplyr::group_by(across(all_of(x_var))) %>%
        dplyr::mutate(group_id = seq(n())) %>%
        dplyr::mutate(y.position = max_y + (max_y-min_y)*step.increase*group_id) %>%
        dplyr::ungroup()

      # Add the manual stat test
      prop_plot <- prop_plot +
        ggpubr::stat_pvalue_manual(stat_test, label = "p.adj.signif", hide.ns = T, label.size = 6, vjust = 0.6, tip.length = 0.02)
    }
  }

  # Return the plot
  return(prop_plot)
}


#' Fix Y-position in rstatix::add_xy_position
#'
#' Modifies the position of significant differences by removing the space of
#' non-significant comparisons.
#'
#' @param stat_test dataframe, as obtained by the test from rstatix.
#' @param x_var string, name of the column to represent in the X-axis.
#'
#' @return dataframe with fixed Y-positions for significant comparisons.
#' @export
fix_pval_y_position <- function(stat_test, x_var) {
  if (all(stat_test$p.adj.signif == "ns")) {
    return(stat_test)
  }

  stat_test <- stat_test %>%
    dplyr::group_by(across(all_of(x_var))) %>%
    dplyr::mutate(y_position_list = list(y.position)) %>%
    dplyr::filter(p.adj.signif != "ns") %>%
    dplyr::mutate(group_id = seq(n())) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y.position = y_position_list[[group_id]]) %>%
    dplyr::select(-y_position_list, -group_id) %>%
    dplyr::ungroup()

  return(stat_test)
}

add_pval_y_position <- function(stat_test, x_var, y_var){
  minmax_y_df <- attributes(stat_test)$args$data %>%
    dplyr::group_by(.data[[x_var]]) %>%
    dplyr::summarise(max_y = max(.data[[y_var]]),
                     min_y = min(.data[[y_var]]))

  stat_test <- stat_test %>%
    dplyr::left_join(minmax_y_df, by = "type")

  return(stat_test)
}

plotJunctionCategoriesMerge <- function(figA, figB, title = "", subtitle = "PPMI dataset"){
  plot_merged <- figA / figB +
    patchwork::plot_annotation(tag_levels = "A",
                               tag_prefix = "Fig.",
                               tag_suffix = ":",
                               title = title,
                               subtitle = subtitle) +
    patchwork::plot_layout(guides = "collect") &
    custom_gg_theme +
    theme(
      strip.text.y = element_blank(),
      plot.tag.position = c(0, 1),
      plot.tag = element_text(hjust = 0, vjust = 0),
    )

  return(plot_merged)
}

modifyScaleColor <- function(plot, grouping_var, experiment = "PPMI"){
  # Modify the scale fill based on the experiment
  if(experiment == "PPMI"){
    plot <- plot +
      scale_fill_manual(
        name = PPMI_group_names[grouping_var],
        values = PPMI_group_colors,
        labels = PPMI_group_labels,
        guide = guide_legend(ncol = 2)
      )
  }else if(experiment == "Wood"){
    plot <- plot +
      scale_fill_manual(
        name = Wood_group_names[grouping_var],
        values = Wood_group_colors,
        labels = Wood_group_labels,
        guide = guide_legend(ncol = 2)
      )
  }else if(experiment == "Hardy"){
    plot <- plot +
      scale_fill_manual(
        name = Hardy_group_names[grouping_var],
        values = Hardy_group_colors,
        labels = Hardy_group_labels,
        guide = guide_legend(ncol = 2)
      )
  }

  return(plot)
}

# Legacy functions needed for some old projects to work

#' Plot the junction categories
#'
#' @param prop_df dataframe, includes the proportions for each sample and
#'   category.
#' @param x_var string, name of the column to represent in the X-axis. Defaults
#'   to "type".
#' @param y_var string, name of the column to represent in the Y-axis. Defaults
#'   to "prop_junc".
#' @param grouping_var string, name of the column in the dataframe to group and
#'   compare the different samples. Defaults to NULL.
#' @param custom_limits numeric vector of two elements, list of lower and upper
#'   limits to the plots. It might not be applied if drawn boxes/significances
#'   are not within these limits. Defaults to NULL.
#' @param y_break numeric vector of two elements. The plot is condensed to
#'   remove the range specified in this argument. The input is expected to be as
#'   "c(0.3, 0.5)". Behaviour unstable, might cause issues when facetting.
#'   Defaults to NULL.
#' @param add_pval boolean, whether to add a Wilcoxon test between the different
#'   groups. Defaults to FALSE.
#' @param ref_group string, name of the reference group to compare against in
#'   the Wilcoxon tests. Defaults to NULL.
#' @param experiment string, global name of the experiment. Employed to color
#'   code the Groups and rename the legend. Defaults to "PPMI".
#'
#' @return ggplot2 object.
#' @export
plotJunctionCategories_legacy <- function(prop_df,
                                   x_var = "type",
                                   y_var = "prop_junc",
                                   grouping_var = NULL,
                                   custom_limits = NULL,
                                   y_break = NULL,
                                   add_pval = F,
                                   ref_group = NULL,
                                   experiment = "PPMI") {
  # Remove unused levels from "x_var" if is factor
  if (is.factor(prop_df[[x_var]])) prop_df <- prop_df %>% dplyr::mutate(!!x_var := droplevels(.data[[x_var]]))

  # Add a column name "split" to determine if the axis will be broken
  y_break_median <- ifelse(is.null(y_break), 0, median(y_break))
  prop_df <- prop_df %>% dplyr::mutate(split = case_when(!is.null(y_break) ~ .data[[y_var]] < y_break_median, .default = T))

  # If no grouping var is provided, a temporary column named "All samples" is
  # created to group all samples into it
  if (is.null(grouping_var)) {
    grouping_var <- "All samples"
    prop_df$`All samples` <- "All samples"
  }

  # Draw the main plot
  prop_plot <- prop_df %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_boxplot(aes(fill = .data[[grouping_var]])) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)), breaks = seq(0, 1, 0.1)) +
    custom_gg_theme

  # Modify the scale fill based on the experiment
  if(experiment == "PPMI"){
    prop_plot <- prop_plot +
      scale_fill_manual(
        name = PPMI_group_names[grouping_var],
        values = PPMI_group_colors,
        labels = PPMI_group_labels,
        guide = guide_legend(ncol = 2)
      )
  }else if(experiment == "Wood"){
    prop_plot <- prop_plot +
      scale_fill_manual(
        name = Wood_group_names[grouping_var],
        values = Wood_group_colors,
        labels = Wood_group_labels,
        guide = guide_legend(ncol = 2)
      )
  }

  # If custom limits are applied, add a "scale_y_continuous" with modified
  # expand
  if (!is.null(custom_limits)) {
    prop_plot <- prop_plot +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), limits = custom_limits, breaks = seq(0, 1, 0.1))
  }

  # If a Y-axis break is added, we add a facet and draw invisible points to
  # force the plot to have the desired custom limits
  if (!is.null(y_break)) {
    # Modify the Y-break values to match with the later "expand" in
    # "scale_y_continuous"
    y_break <- y_break / c(1.07, 1)

    # Create a dataframe with points that will be drawn to adjust the graph
    # limits
    custom_limits_df <- data.frame(
      x = prop_df[[x_var]][1],
      y = c(y_break, custom_limits)
    ) %>%
      dplyr::mutate(split = y < median(y_break)) %>%
      dplyr::rename(!!y_var := y, !!x_var := x)

    # Add the invisible dots and use a facet to split the graph, breaking the
    # axis. It is neccesary to remove the "strip.text.y"
    prop_plot <- prop_plot +
      geom_point(data = custom_limits_df, shape = NA) +
      # facet_grid(.data[[y_var]] < median(y_break) ~ ., scales = "free_y", space = "free_y") +
      facet_grid(split ~ ., scales = "free_y", space = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = NULL, breaks = seq(0, 1, 0.1)) +
      theme(strip.text.y = element_blank())
  }

  # Add p-values of the Wilcoxon comparison between the groping_var
  if (add_pval) {
    # We do a manual test first that later will be added to the plot. It is
    # important to add X and Y positions of the brackets, and to correct the
    # Y-position to ignore non-significant elements. There is probably an option
    # to ignore them in the default functions, but I have not found an automatic
    # way to do so.
    stat_test <- prop_df %>%
      dplyr::group_by(type, split) %>%
      rstatix::wilcox_test(as.formula(paste0(y_var, "~", grouping_var)), ref.group = ref_group) %>%
      rstatix::adjust_pvalue(method = "BH") %>%
      rstatix::add_significance(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns")) %>%
      # add_y_position(step.increase = 0.06)
      rstatix::add_xy_position(x = x_var, step.increase = 0.06) %>%
      fix_pval_y_position(x_var)

    # Add the manual stat test
    prop_plot <- prop_plot +
      ggpubr::stat_pvalue_manual(stat_test, label = "p.adj.signif", hide.ns = T, label.size = 6, vjust = 0.6, tip.length = 0.02)

    if(!is.null(custom_limits)){
      if(max(stat_test$y.position) > custom_limits[2]){
        prop_plot <- prop_plot +
          scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)), breaks = seq(0, 1, 0.1))
      }
    }
  }

  # Return the plot
  return(prop_plot)
}


#' Plot the junction categories with facet
#'
#' @param prop_df dataframe, includes the proportions for each sample and
#'   category.
#' @param x_var string, name of the column to represent in the X-axis. Defaults
#'   to "type".
#' @param y_var string, name of the column to represent in the Y-axis. Defaults
#'   to "prop_junc".
#' @param grouping_var string, name of the column in the dataframe to group and
#'   compare the different samples. Defaults to NULL.
#' @param split_var string, name of the column in the dataframe to split the
#'   graph in facets. Defaults to NULL.
#' @param add_pval boolean, whether to add a Wilcoxon test between the different
#'   groups. Defaults to FALSE.
#' @param ref_group string, name of the reference group to compare against in
#'   the Wilcoxon tests. Defaults to NULL.
#' @param n_row numeric, number of rows in the facet.
#' @param experiment string, global name of the experiment. Employed to color
#'   code the Groups and rename the legend. Defaults to "PPMI".
#'
#' @return ggplot2 object.
#' @export
plotJunctionCategoriesSplit_legacy <- function(prop_df,
                                        x_var = "type",
                                        y_var = "prop_junc",
                                        grouping_var = NULL,
                                        split_var = NULL,
                                        add_pval = F,
                                        ref_group = NULL,
                                        n_row = NULL,
                                        step.increase = 0.18,
                                        experiment = "PPMI") {
  # If no grouping var is provided, a temporary column named "All samples" is
  # created to group all samples into it
  if (is.null(grouping_var)) {
    grouping_var <- "All samples"
    prop_df$`All samples` <- "All samples"
  }

  # Draw the main plot with a facet split
  prop_plot <- prop_df %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_boxplot(aes(fill = .data[[grouping_var]])) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    facet_wrap(reformulate(split_var), scales = "free", nrow = n_row) +
    custom_gg_theme

  # Modify the scale fill based on the experiment
  if(experiment == "PPMI"){
    prop_plot <- prop_plot +
      scale_fill_manual(
        name = PPMI_group_names[grouping_var],
        values = PPMI_group_colors,
        labels = PPMI_group_labels,
        guide = guide_legend(ncol = 2)
      )
  }else if(experiment == "Wood"){
    prop_plot <- prop_plot +
      scale_fill_manual(
        name = Wood_group_names[grouping_var],
        values = Wood_group_colors,
        labels = Wood_group_labels,
        guide = guide_legend(ncol = 2)
      )
  }

  # Add p-values of the Wilcoxon comparison between the groping_var
  if (add_pval) {
    # We do a manual test first that later will be added to the plot. It is
    # important to add X and Y positions of the brackets, and to correct the
    # Y-position to ignore non-significant elements. There is probably an option
    # to ignore them in the default functions, but I have not found an automatic
    # way to do so.
    stat_test <- prop_df %>%
      dplyr::group_by(type) %>%
      rstatix::wilcox_test(as.formula(paste0(y_var, "~", grouping_var)), ref.group = ref_group) %>%
      rstatix::adjust_pvalue(method = "BH") %>%
      rstatix::add_significance(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns")) %>%
      rstatix::add_y_position(fun = "max", step.increase = 0, scales = "free_y") %>% dplyr::rename(max_y = y.position) %>%
      rstatix::add_y_position(fun = "min", step.increase = 0, scales = "free_y") %>% dplyr::rename(min_y = y.position) %>%
      rstatix::add_x_position(x = x_var) %>%
      dplyr::mutate(xmin = xmin - x + 1, xmax = xmax - x + 1) %>%
      dplyr::filter(p.adj.signif != "ns")

    if(nrow(stat_test) > 0){
      stat_test <- stat_test %>%
        dplyr::group_by(across(all_of(x_var))) %>%
        dplyr::mutate(group_id = seq(n())) %>%
        dplyr::mutate(y.position = max_y + (max_y-min_y)*step.increase*group_id) %>%
        dplyr::ungroup()

      # Add the manual stat test
      prop_plot <- prop_plot +
        ggpubr::stat_pvalue_manual(stat_test, label = "p.adj.signif", hide.ns = T, label.size = 6, vjust = 0.6, tip.length = 0.02)
    }
  }

  # Return the plot
  return(prop_plot)
}