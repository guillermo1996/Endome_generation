#' Covariate Correct the logcounts
#'
#' Apply the `limma::removeBatchEffect()` function to remove the contribution
#' from covariates to the log2(fpm + 1) feature counts. Genes that do not pass
#' the fpm > 1 in 50% if the samples threshold are removed.
#'
#' @param raw_data dataframe, feature counts.
#' @param meta dataframe, metadata of the samples.
#' @param covariates character vector, list of covariates to correct the
#'   logcounts.
#' @param blacklsited_genes character vector, list of genes from the ENCODE
#'   blacklisted regions. If not provided, no genes are removed. By default,
#'   NULL.
#' @param overwrite boolean, whether to overwrite previously generated results
#'   from the function. If set to FALSE, the function looks for the files in
#'   disk and loads them if possible. By default, FALSE.
#' @param logcounts_path string, path to store the uncorrected log counts.
#'   Defaults to "".
#' @param logcounts_adj_path string, path to store the corrected log counts.
#'   Defaults to "".
#'
#' @return list with the corrected and uncorrected logcounts dataframes.
#' @export
logcountsCovariateCorrection <- function(raw_data,
                                         meta,
                                         experimental_condition,
                                         covariates,
                                         blacklisted_genes = NULL,
                                         overwrite = F,
                                         logcounts_path = "",
                                         logcounts_adj_path = ""){
  # Check if output already exists and the "overwrite" parameter. If the file is
  # found in disk and "overwrite" is set to F, then we read the previously
  # generated output of the function.
  if(file_test("-f", logcounts_adj_path) & !overwrite){
    if(file_test("-f", logcounts_path)){
      logcounts <- readRDS(logcounts_path)
    }else{
      logcounts <- NULL
    }

    logcounts_adj <- readRDS(logcounts_adj_path)
    return(list(logcounts = logcounts, logcounts_adj = logcounts_adj))
  }

  # Ensure matching row names and sample IDs (not needed)
  if(!identical(colnames(raw_data), meta$sample_id)){
    raw_data <- raw_data[, intersect(meta$sample_id, colnames(raw_data))]
    meta <- meta[match(colnames(raw_data), meta$sample_id), ]
  }
  stopifnot(identical(colnames(raw_data), meta$sample_id))

  # Prepare the design matrix. The "0" term before the experimental condition
  # facilitates the later creation of contrasts, while not affecting the final
  # results. As explained "If experimental_condition is a factor, then the
  # models with and without the intercept term are equivalent" (source: Section
  # 4.3 of
  # https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)
  meta_covs <- meta %>%
    dplyr::select(sample_id, all_of(experimental_condition), all_of(covariates)) %>%
    dplyr::mutate(across(where(is.numeric), scale),
                  across(where(is.character), as.factor)) %>%
    tibble::column_to_rownames("sample_id")

  design = model.matrix(reformulate(c(0, experimental_condition, covariates)), data = meta_covs)

  ## Remove columns that introduce linear combinations
  linear_combos <- caret::findLinearCombos(design)
  if(!is.null(linear_combos$remove)) design <- design[, -linear_combos$remove]

  # Prepare the feature count matrix
  dds <- DESeq2::DESeqDataSetFromMatrix(raw_data, colData = meta_covs, design = ~ 1)

  ## Remove genes from the blacklisted regions
  dds <- dds[!rownames(dds) %in% blacklisted_genes, ]

  ## Remove low expressed genes and apply log2(fpm + 1)
  isexpr <- rowSums(DESeq2::fpm(dds) > 1) >= 0.5 * ncol(dds)
  logcounts = log2(DESeq2::fpm(dds)[isexpr, ] + 1)

  ## Correct the feature counts
  design_treatment <- design[, grepl(paste(c(experimental_condition, "\\(Intercept\\)"), collapse = "|"), colnames(design)), drop = F]
  design_batch <- design[, !colnames(design) %in% colnames(design_treatment), drop = F]
  logcounts_adj <- limma::removeBatchEffect(logcounts, covariates = design_batch, design = design_treatment)

  # Store the logcounts in disk
  if(logcounts_path != "") logcounts %>% saveRDS(logcounts_path)
  if(logcounts_adj_path != "") logcounts_adj %>% saveRDS(logcounts_adj_path)

  # Return both the corrected and uncorrected logcounts
  return(list(logcounts = logcounts, logcounts_adj = logcounts_adj))
}

getOutliersPCA <- function(pca,
                           var_explained = 85,
                           max_pc = NULL,
                           set_pc = NULL,
                           z_score_threshold = 3){
  highest_pc <- which(summary(pca)$importance[3, ] > var_explained/100)[1]
  if(!is.null(max_pc)) highest_pc <- min(highest_pc, max_pc)
  if(!is.null(set_pc)) highest_pc <- set_pc

  pc_outliers <- pca$x %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    .[, seq(highest_pc + 1)] %>%
    dplyr::mutate(across(paste0("PC", seq(highest_pc)), list(zscore = ~abs((. - mean(.))/sd(.))))) %>%
    dplyr::mutate(max_outlier_pc = highest_pc,
                  z_score_pca_outlier = if_any(ends_with("zscore"), ~ . > z_score_threshold))

  return(pc_outliers)
}

# Function extracted from "autoelbow" python package
# https://pypi.org/project/autoelbow/
elbow_finder <- function(x_values, y_values) {
  calc_distance <- function(x1, y1, a, b, c){
    d = abs((a*x1 + b*y1 + c))/(sqrt(a*a + b*b))
    return(d)
  }

  a <- head(y_values, 1) - tail(y_values, 1)
  b <- tail(x_values, 1) - head(x_values, 1)
  c1 <- head(x_values, 1) - tail(y_values, 1)
  c2 <- tail(x_values, 1) - head(y_values, 1)
  c <- c1 - c2

  distance <- c()
  for(i in seq_along(x_values)){
    distance <- c(distance, calc_distance(x_values[i], y_values[i], a, b, c))
  }

  return(which.max(distance))
}

getOutliersWGCNA <- function(quantlog, z_score_threshold = 2){
  mat <- quantlog

  # Calculate biweight midcorrelation matrix between samples
  bicor.mat <- WGCNA::bicor(mat)

  # Calculate signed adjacency matrix between samples
  normadj.mat <- (0.5 + 0.5*bicor.mat)^2

  # Calculate connectivity
  connectivity.mat <- WGCNA::fundamentalNetworkConcepts(normadj.mat)
  ku <- connectivity.mat$Connectivity
  z.ku <- (ku - mean(ku))/sqrt(var(ku))

  # Extract results
  sample_connectivity_res <- z.ku %>%
    stack() %>%
    tibble::as_tibble() %>%
    dplyr::rename(sample_id = ind,
                  wgcna_connectivity_z_score = values) %>%
    dplyr::mutate(z_score_wgcna_outlier = abs(wgcna_connectivity_z_score) > z_score_threshold) %>%
    dplyr::relocate(sample_id, wgcna_connectivity_z_score)

  return(sample_connectivity_res)
}

plotPCgroup <- function(pca,
                        meta,
                        covariate,
                        number_of_PCs = 8,
                        overlap_pairs = F){
  # Sub-function to plot based on the PCA and two selected PCs
  plotPCpair <- function(PC_pair = c(pc_a, pc_b), pca_group, covariate = NULL){
    pca_plot <- pca_group %>%
      ggplot(aes(x = .data[[PC_pair[1]]], y = .data[[PC_pair[2]]])) +
      geom_point(aes(fill = .data[[covariate]]), alpha = 0.85, shape = 21, size = 2.5) +
      custom_gg_theme +
      theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

    if(is.numeric(pca_group[[covariate]])){
      pca_plot <- pca_plot +
        ggsci::scale_fill_gsea(name = covariate) +
        guides(fill = guide_colourbar(frame.colour = "black", ticks.colour = "black",
                                      ticks.linewidth = 0.5, frame.linewidth = 0.5,
                                      barwidth = 10, title.vjust = 0.82))
    }else{
      pca_plot <- pca_plot +
        ggsci::scale_fill_simpsons(name = covariate) +
        guides(fill = guide_legend(override.aes = list(size = 4)))
    }

    return(pca_plot)
  }

  # Combine the PCA with the metadata to have access to the covariate
  pca_group <- pca$x %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    dplyr::select(sample_id, all_of(paste0("PC", seq(number_of_PCs)))) %>%
    dplyr::left_join(meta, by = "sample_id")

  # List of available PCs
  PC_list <- grep("^PC\\d+", colnames(pca_group), value = T)

  # Let's say we have the first four PCs: 1, 2, 3, 4. We have two approaches to
  # group them for plotting: i) overlapping pairs would make 3 plots (1, 2), (2,
  # 3), (3, 4) and ii) non-overlapping paris would make 2 plots (1, 2) and (3,
  # 4). Which approach to take is decided by the 'overlap_pairs' parameter.
  if(overlap_pairs){
    PC_pairs <- mapply(function(x, y) list(c(x, y)), PC_list[-length(PC_list)], PC_list[-1])
  }else{
    PC_pairs <- split(PC_list, sort(rep(seq(length(PC_list)/2), 2))) %>% setNames(., lapply(., function(x) x[[1]]))
  }

  # Use patchwork to join all figures
  pca_plot <- lapply(PC_pairs, FUN = plotPCpair, pca_group, covariate) %>%
    patchwork::wrap_plots(guides = "collect") &
    theme(legend.position = "top",
          axis.title.y = element_text(margin = margin(0, 3, 0, 13)),
          axis.title.x = element_text(margin = margin(3, 0, 13, 0)))

  return(pca_plot)
}