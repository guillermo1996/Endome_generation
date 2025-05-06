## _________________________________________________
##
## Helper Functions: Notebooks
##
## Aim: Include into a single file all the functions required to develop
## consistent notebooks across different analyses.
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 18/01/2024
##
## Copyright (c) Guillermo Rocamora Pérez, 2024
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (18/01/2024)
## _________________________________________________
##
## - Notes:
##
## Many of these functions might require modifications per specific project
## requirements. These are just baselines to follow but not definite versions of
## them.
##
##
## - Changelog:
##
## + v1.0: initial file.
##
##
## - Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

#' Print table for a dataframe
#'
#' This function is a wrapper that includes everything needed for consistent
#' tables across different reports and analyses. Many parameters can be employed
#' to precisely control how the table should look like.
#'
#' @param df dataframe to print.
#' @param style string, style to represent the table. Accepted values are "html"
#'   and "md", each one with a different backend ("DT" and "kableExtra"
#'   respectively). Defaults to "html".
#' @param limit numeric, maximum number of rows to show. Defaults to NULL, or no limit.
#' @param random boolean, whether to shuffle the rows before printing. Default to FALSE.
#' @param seed numeric, seed to use if the dataframe is shuffled. Defaults to NULL, or no seed.
#' @param order see \code{DT::datatable()}.
#' @param autoWidth see \code{DT::datatable()}.
#' @param pageLength see \code{DT::datatable()}.
#' @param rownames see \code{DT::datatable()}.
#' @param compact boolean, whether to include "compact" in the DT classes. See
#'   \code{DT::datatable()}.
#' @param rowGroup see \code{DT::datatable()}.
#' @param hideRowGroup boolean, when grouping rows, whether to include the
#'   grouping field. Defaults to FALSE.
#' @param dom see \code{DT::datatable()}. Defaults to "tpf".
#' @param full_width see \code{kableExtra::kable_classic}.
#'
#' @return prints the input dataframe as a table.
#' @export
printPrettyDf <- function(df, style = "html", 
                          seed = NULL, random = F, limit = NULL,
                          compact = T, autoWidth = F, pageLength = 10, order = list(), rownames = F, dom = "ftipr", rowGroup = NULL, hideRowGroup = F,
                          full_width = T){
  # Verify the style
  if(!style %in% c("html", "md")) stop("Only valid styles are 'html' and 'md'.")
  
  # Limit or randomize the dataframe based on seed/random/limit parameters
  if(!is.null(seed)) set.seed(seed)
  if(random) df <- df %>% dplyr::slice_sample(prop = 1)
  if(!is.null(limit)) df <- df %>% dplyr::slice_head(n = limit)
  
  # If HTML style is provided, DT::datatable is employed
  if(style == "html"){
    # Generate the class of the datatable
    dt_class <- "display cell-border nowrap"
    if(compact) dt_class <- stringr::str_c(dt_class, "compact", sep = " ")
    
    # Generate the options of the datatable
    DT_options <- list(scrollX = T,
                       bLengthChange = T,
                       autoWidth = autoWidth,
                       pageLength = pageLength,
                       order = order,
                       dom = dom)
    
    # Generate the datatable
    dt_object <- DT::datatable(df, class = dt_class, rownames = rownames, options = DT_options)
    
    # If rowGroup is provided, add the options and regenerate the datatable
    if(!is.null(rowGroup)){
      DT_options$rowGroup <- list(dataSrc = c(rowGroup))
      DT_options$columnDefs <- list(list(visible = F, targets = c(ifelse(hideRowGroup, rowGroup, -1))))
      dt_object <- DT::datatable(df, class = dt_class, rownames = rownames, options = DT_options, extension = "RowGroup")
    }
    
    if(requireNamespace("knitr", quietly = TRUE) && tryCatch(knitr::is_html_output(), error = function(e) FALSE)) return(htmltools::tagList(dt_object))
    return(dt_object)
  }
  
  ## If MD style is provided, kableExtra::kbl is employed
  if(style == "md"){
    kable_object <- kableExtra::kbl(df, booktabs = T, linesep = "") %>% 
      kableExtra::kable_classic(full_width = full_width, "hover", "striped", html_font = "Cambria", font_size = 14) %>% 
      kableExtra::row_spec(0, bold = T, font_size = 16)
    
    return(kable_object)
  }
  
  return(NULL)
}

# printPrettyDf <- function(df,
#                           style = "html",
#                           limit = NULL,
#                           random = F,
#                           seed = NULL,
#                           order = list(),
#                           autoWidth = F,
#                           pageLength = 10,
#                           rownames = F,
#                           compact = T,
#                           rowGroup = NULL,
#                           hideRowGroup = F,
#                           dom = "tpf",
#                           full_width = T) {
#   if (!style %in% c("html", "md")) {
#     stop("Only valid styles are 'html' and 'md'")
#   }
# 
#   if (random) {
#     if (!is.null(seed)) set.seed(seed)
#     df <- if (!is.null(limit)) df %>% dplyr::slice_sample(n = limit) else df %>% dplyr::slice_sample(prop = 1)
#   }
#   if (!is.null(limit)) df <- df %>% dplyr::slice_head(n = limit)
# 
#   if (style == "html") {
#     dt_class <- "display cell-border nowrap"
#     if (compact) dt_class <- stringr::str_c(dt_class, "compact", sep = " ")
# 
#     if (is.null(rowGroup)) {
#       return(DT::datatable(df,
#         options = list(
#           scrollX = T,
#           autoWidth = autoWidth,
#           bLengthChange = F,
#           pageLength = pageLength,
#           order = order,
#           dom = dom
#         ),
#         class = dt_class, rownames = rownames
#       ))
#     } else {
#       hideRow <- if (hideRowGroup) rowGroup else -1
# 
#       return(DT::datatable(df,
#         extensions = "RowGroup",
#         options = list(
#           rowGroup = list(dataSrc = c(rowGroup)),
#           bLengthChange = F,
#           columnDefs = list(list(visible = F, targets = c(hideRow))),
#           scrollX = T,
#           autoWidth = autoWidth,
#           pageLength = pageLength,
#           order = order,
#           dom = dom
#         ),
#         class = dt_class, rownames = rownames
#       ))
#     }
#   }
# 
#   if (style == "md") {
#     df %>%
#       kableExtra::kbl(booktabs = T, linesep = "") %>%
#       kableExtra::kable_classic(full_width = full_width, "hover", "striped", html_font = "Cambria", font_size = 14) %>%
#       kableExtra::row_spec(0, bold = T, font_size = 16) %>%
#       return()
#   }
# }
