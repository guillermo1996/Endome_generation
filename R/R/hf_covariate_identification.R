#' Extract correlation between numerical covariates and the PCs
#'
#' Adapted from the covariate identification script by Aine.
#'
#' @param pca_out dataframe, PCs of the studied samples.
#' @param meta_clean dataframe, metadata of the studied samples.
#' @param highest_pc numeric, number of PCs to employ.
#'
#' @return dataframe with the correlation between numerical covariates and the
#'   PCs.
#' @export
extractNumCor <- function(pca_out, meta_clean, highest_pc = 10) {
  ## A: Spearman Rho values:
  cors.r.numericVars <- pca_out %>%
    dplyr::select(sample_id, all_of(paste0("PC", seq(highest_pc)))) %>%
    dplyr::left_join(meta_clean %>% dplyr::select(sample_id, where(is.numeric)), by = "sample_id") %>%
    dplyr::select(-sample_id) %>%
    as.matrix() %>%
    Hmisc::rcorr(type = "spearman") %>%
    .[["r"]] %>%
    as.data.frame() %>%
    dplyr::select(matches("^PC\\d+")) %>%
    tibble::rownames_to_column("meta_var") %>%
    dplyr::filter(!grepl("^PC\\d+", meta_var)) %>%
    tidyr::pivot_longer(cols = -meta_var, names_to = "name", values_to = "value") %>%
    dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name)))) %>%
    dplyr::group_by(meta_var) %>%
    dplyr::mutate(max.PC.cor = value[which.max(abs(value))], mean.PC.cor = mean(value)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(PC = name)

  ## A: P-values
  cors.p.numericVars <- pca_out %>%
    dplyr::select(sample_id, all_of(paste0("PC", seq(highest_pc)))) %>%
    dplyr::left_join(meta_clean %>% dplyr::select(sample_id, where(is.numeric)), by = "sample_id") %>%
    dplyr::select(-sample_id) %>%
    as.matrix() %>%
    Hmisc::rcorr(type = "spearman") %>%
    .[["P"]] %>%
    as.data.frame() %>%
    dplyr::select(matches("^PC\\d+")) %>%
    tibble::rownames_to_column("meta_var") %>%
    dplyr::filter(!grepl("^PC\\d+", meta_var)) %>%
    tidyr::pivot_longer(cols = -meta_var, names_to = "name", values_to = "value") %>%
    dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name)))) %>%
    dplyr::rename(PC = name)

  ## A: bind p and stat
  cors.num <- dplyr::left_join(cors.p.numericVars %>% dplyr::rename(p = value),
    cors.r.numericVars %>% dplyr::rename(stat = value),
    by = c("meta_var", "PC")
  ) %>%
    dplyr::mutate(var_type = "continuous") %>%
    dplyr::group_by(meta_var) %>%
    dplyr::mutate(min.pvalue = min(p)) %>%
    dplyr::relocate(var_type, .after = last_col())

  return(cors.num)
}


#' Extract correlation between categorical covariates and the PCs
#'
#' Adapted from the covariate identification script by Aine.
#'
#' @param pca_out dataframe, PCs of the studied samples.
#' @param meta_clean dataframe, metadata of the studied samples.
#' @param highest_pc numeric, number of PCs to employ.
#'
#' @return dataframe with the correlation between categorical covariates and the
#'   PCs.
#' @export
extractCatCor <- function(pca_out, meta_clean, highest_pc) {
  ## A: get PC/categorical correlations using K-W Chi-squared
  ## G: added some regular expression to avoid issues with columns that are not
  ## PCs.
  cors.cat <- pca_out %>%
    dplyr::select(sample_id, all_of(paste0("PC", seq(highest_pc)))) %>%
    dplyr::left_join(meta_clean %>% dplyr::select(sample_id, !where(is.numeric)), by = "sample_id") %>%
    dplyr::select(-sample_id) %>%
    run_stat_test_for_all_cols_get_res(in_df = ., stat_test = "kruskal") %>%
    dplyr::filter(
      grepl("^PC\\d+", var1),
      !grepl("^PC\\d+", var2)
    ) %>%
    dplyr::select(PC = var1, meta_var = var2, p, stat) %>%
    dplyr::group_by(meta_var) %>%
    dplyr::mutate(
      max.PC.cor = stat[which.max(abs(stat))],
      mean.PC.cor = mean(stat),
      min.pvalue = min(p)
    ) %>%
    dplyr::mutate(PC = as.factor(as.numeric(gsub("PC", "", PC)))) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(meta_var, PC) %>%
    dplyr::mutate(var_type = "categorical")

  return(cors.cat)
}


#' Merge categorical and numerical correlation with PCs
#'
#' @param pca prcomp object.
#' @param cors.num dataframe, correlation between numerical covariates and PCs.
#' @param cors.cat dataframe, correlation between categorical covariates and
#'   PCs.
#'
#' @return dataframe with the correlation between covariates and the PCs.
#' @export
mergeNumCatCor <- function(pca, cors.num, cors.cat) {
  ## A: generate dataframe containing covariate, PC axis, rho & P-value
  ## corresponding to PC-covariate, variance explained by each PC first, binding
  ## cont (spearman) and cat (K-S) data to get df containing stat assessment for
  ## all putative covariates
  cors.all <- dplyr::bind_rows(cors.num, cors.cat) %>%
    dplyr::mutate(PC = as.factor(PC)) %>%
    dplyr::mutate(
      padj = p.adjust(p, method = "fdr"),
      is.sig = case_when(padj < 0.05 ~ stat,
        .default = NA_integer_
      )
    ) %>%
    # A: add variance explained
    dplyr::left_join(
      summary(pca)$importance %>%
        t() %>%
        tibble::as_tibble(rownames = "PC") %>%
        janitor::clean_names() %>%
        dplyr::rename(PC = pc) %>%
        dplyr::mutate(PC = as.factor(gsub("PC", "", PC))),
      by = "PC"
    ) %>%
    # A: calculate weighted correlation - weights PC-covariate correlation
    # coefficient by the % variance explained by that PC
    dplyr::mutate(variance_weighted_correlation = proportion_of_variance * (stat**2)) %>%
    # A: add adjusted P-value - adjust for number of putative covariates
    dplyr::group_by(PC) %>%
    dplyr::mutate(padj = p.adjust(p, "fdr")) %>%
    dplyr::ungroup()

  return(cors.all)
}


#' Run stat test between dataframe columns
#'
#' Initially developed by Aine Fairbrother-Browne. This function calculates the
#' collinearity between categorical ~ categorical and numerical ~ categorical
#' covariates.
#'
#' @param in_df dataframe.
#' @param stat_test string, either "kruskal" for numerical ~ categorical or
#'   "chisq" for categorical ~ categorical.
#' @param colnames1 string vector, names of the categorical/numerical columns.
#'   Can be left as NULL.
#' @param colnames2 string vector, names of the categorical columns. Can be left
#'   as NULL.
#'
#' @return dataframe with the collinearity between the columns.
#' @export
run_stat_test_for_all_cols_get_res <- function(in_df, stat_test, colnames1 = NULL, colnames2 = NULL) {
  if (is.null(colnames1) | is.null(colnames2)) {
    colnames1 <- colnames(in_df)
    colnames2 <- colnames(in_df)
  }

  collect_stat <- expand.grid(colnames1, colnames2) %>%
    as.data.frame() %>%
    dplyr::mutate(value = NA) %>%
    tidyr::pivot_wider(names_from = "Var2", values_from = "value") %>%
    tibble::column_to_rownames("Var1")

  collect_p <- expand.grid(colnames1, colnames2) %>%
    as.data.frame() %>%
    dplyr::mutate(value = NA) %>%
    tidyr::pivot_wider(names_from = "Var2", values_from = "value") %>%
    tibble::column_to_rownames("Var1")

  for (col_i in colnames1) {
    for (col_j in colnames2) {
      if (stat_test == "kruskal") {
        # A: run the test - numeric ~ categorical
        res <- kruskal.test(in_df[[col_i]] ~ in_df[[col_j]])

        # A: store result
        collect_stat[col_i, col_j] <- res$statistic %>% as.numeric()
        collect_p[col_i, col_j] <- res$p.value %>% as.numeric()
      }

      if (stat_test == "chisq") {
        # A: run the test - categorical ~ categorical
        res <- stats::chisq.test(in_df[[col_i]], in_df[[col_j]])

        # A: store result
        collect_stat[col_i, col_j] <- res$statistic %>% as.numeric()
        collect_p[col_i, col_j] <- res$p.value %>% as.numeric()
      }
    }
  }

  # A: bind the stat and p-value into a single dataframe
  out_df <- dplyr::left_join(
    collect_p %>%
      tibble::rownames_to_column("var1") %>%
      tidyr::pivot_longer(2:ncol(.), names_to = "var2", values_to = "p"),
    collect_stat %>%
      tibble::rownames_to_column("var1") %>%
      tidyr::pivot_longer(2:ncol(.), names_to = "var2", values_to = "stat"),
    by = c("var1", "var2")
  ) %>%
    dplyr::filter(!(var1 == var2)) %>%
    dplyr::mutate(test_run = stat_test)

  return(out_df)
}
