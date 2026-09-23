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
## Latest version: v1.3 (26/08/2026)
## _________________________________________________
##
## - Notes:
#
## - Changelog:
## 
##    + v1.3 (26/08/2026): Improved write_functions to report the helper 
##    functions in the notebooks
##    + v1.2 (30/09/2025): Added write_functions function.
##    + v1.1 (28/05/2025): Updated to improve consistency with the DRI Server.
##    + v1.0: initial file.
##
## - Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________


################################################################################
# ---- 1. Required libraries ----
suppressWarnings(suppressMessages(library(kableExtra)))
suppressWarnings(suppressMessages(library(DT)))

################################################################################
# ---- 2. Print tables ----

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
renderTable <- function(df, style = "auto",
                        # --- shared settings
                        caption = NULL, digits = 3, big.mark = ",", align = NULL, escape = TRUE, exclude_cols = NULL, 
                        seed = NULL, random = FALSE, limit = NULL,
                        # --- DT settings
                        compact = TRUE, autoWidth = FALSE, pageLength = 10, order = list(), rownames = FALSE, dom = "ftipr",
                        rowGroup = NULL, hideRowGroup = FALSE, filter = "none", scrollX = TRUE, scrollY = NULL,
                        lengthMenu = c(10, 25, 50, 100),
                        # --- kableExtra settings
                        theme = "bootstrap", full_width = TRUE, font_size = 14, html_font = NULL, row_padding = "3x"){
  ## ---- 0. Input handling ---------------------------------------------------
  if(is.null(df)) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  style <- match.arg(tolower(as.character(style)), c("auto", "html", "md"))
  backend <- switch(style,
                    html = "dt",
                    md = "kable",
                    auto = if(nrow(df) > 20) "dt" else "kable")
  if(backend == "dt" && !requireNamespace("DT", quietly = TRUE)) {
    warning("Package 'DT' is not installed; falling back to kableExtra.", call. = FALSE)
    backend <- "kable"
  }
  
  ## ---- 1. Subset / shuffle -------------------------------------------------
  if(!is.null(seed)) set.seed(seed)
  if(isTRUE(random)) df <- df %>% dplyr::slice_sample(prop = 1)
  if(!is.null(limit) && limit < nrow(df)) df <- df %>% dplyr::slice_head(n = limit)
  
  num_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], exclude_cols)
  int_cols <- num_cols[vapply(df[num_cols], .isWholeCol, logical(1))]
  dec_cols <- setdiff(num_cols, int_cols)
  
  ## Columns opted out of formatting: pass them through untouched
  if(length(exclude_cols) > 0){
    keep <- intersect(exclude_cols, names(df))
    for(cl in keep) if (is.numeric(df[[cl]])) df[[cl]] <- as.character(df[[cl]])
    if(is.null(align) && backend == "kable") align <- ifelse(names(df) %in% keep | vapply(df, is.numeric, logical(1)), "r", "l")
  }
  
  ## ---- 2. DT backend -------------------------------------------------------
  if(backend == "dt"){
    dt_class <- "display cell-border nowrap"
    if(isTRUE(compact)) dt_class <- paste(dt_class, "compact")
    extensions <- character(0)
    
    dt_options <- list(
      scrollX = scrollX,
      bLengthChange = TRUE,
      autoWidth = autoWidth,
      pageLength = pageLength,
      lengthMenu = lengthMenu,
      order = order,
      dom = dom,
      deferRender = TRUE
    )
    if(!is.null(scrollY)) dt_options$scrollY <- scrollY
    
    ## Row grouping: accept a column name or a 0-based index
    if (!is.null(rowGroup)) {
      extensions <- c(extensions, "RowGroup")
      grp_idx <- .dtColIndex(rowGroup, df, rownames)
      dt_options$rowGroup <- list(dataSrc = grp_idx)
      # Rows must be sorted by the grouping column for RowGroup to make sense
      if(length(order) == 0) dt_options$order <- list(list(grp_idx, "asc"))
      if(isTRUE(hideRowGroup)) dt_options$columnDefs <- list(list(visible = FALSE, targets = grp_idx))
    }
    
    dt_object <- DT::datatable(
      df, 
      class = dt_class, rownames= rownames, caption = caption,
      filter = filter, escape = escape, extensions = extensions,
      options = dt_options)
    
    ## Format numbers on render, so sorting still uses the underlying values.
    if(length(int_cols) > 0){
      dt_object <- DT::formatRound(dt_object, columns = int_cols, digits = 0, interval = 3, mark = big.mark)
    }
    if(length(dec_cols) > 0) {
      dt_object <- if(is.null(digits)){
        DT::formatCurrency(dt_object, columns = dec_cols, currency = "", digits = NULL, interval = 3, mark = big.mark)
      }else{
        DT::formatRound(dt_object, columns = dec_cols, digits = digits, interval = 3, mark = big.mark)
      }
    }
    
    if (.inKnitr()) return(htmltools::tagList(dt_object))
    return(dt_object)
  }
  
  
  ## ---- 3. kableExtra backend ----------------------------------------------
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("Package 'kableExtra' is required for the kable backend.", call. = FALSE)
  }
  
  fmt <- "html"
  theme <- match.arg(tolower(theme), c("bootstrap", "classic"))
  
  kable_object <- kableExtra::kbl(
    df, format = fmt, caption = caption, align = align,
    digits = if (is.null(digits)) getOption("digits") else digits,
    format.args = list(big.mark = big.mark),
    escape = escape,
    booktabs = TRUE,
    linesep = ""
  )
  
  if(theme == "bootstrap"){
    kable_object <- kableExtra::kable_styling(
      kable_object,
      bootstrap_options = c("striped", "hover", "condensed", "bordered"),
      latex_options = c("striped", "hold_position"),
      full_width = full_width,
      position = "left",
      font_size = font_size
    )
  }else{
    kable_object <- kableExtra::kable_classic(
      kable_object,
      lightable_options = c("striped", "hover"),
      full_width = full_width,
      position = "left",
      html_font = if (is.null(html_font)) "Cambria" else html_font,
      font_size = font_size
    )
  }
  
  # Make the table more compact
  pad_css <- if(isTRUE(compact) && !is.null(row_padding)){
    sprintf("padding-top: %s; padding-bottom: %s; line-height: 1.2;", row_padding, row_padding)
  } else NULL
  
  kable_object <- kableExtra::row_spec(kable_object, 0, bold = TRUE, font_size = font_size + 2)
  if(!is.null(pad_css)) kable_object <- kableExtra::row_spec(kable_object, seq_len(nrow(df)), extra_css = pad_css)
  
  # Long static tables get a scroll box instead of running off the page
  if (fmt == "html" && nrow(df) > 25) {
    kable_object <- kableExtra::scroll_box(kable_object, width = "100%", height = "500px")
  }
  
  kable_object
}

# TRUE when knitting (any format)
.inKnitr <- function() {
  requireNamespace("knitr", quietly = TRUE) && !is.null(knitr::opts_knit$get("out.format"))
}

# TRUE for integer columns and for doubles that only hold whole values
.isWholeCol <- function(x) {
  if (is.integer(x)) return(TRUE)
  if (!is.numeric(x)) return(FALSE)
  x <- x[is.finite(x)]                     # ignore NA/Inf when deciding
  length(x) == 0 || all(x == round(x))
}


# DT column indices are 0-based and shift by one when rownames are shown
.dtColIndex <- function(col, df, rownames){
  offset <- if(isTRUE(rownames)) 1 else 0
  if(is.character(col)){
    idx <- match(col, names(df))
    if(any(is.na(idx))) stop("rowGroup column not found in the dataframe: ", paste(col[is.na(idx)], collapse = ", "), call. = FALSE)
    return(as.integer(idx - 1 + offset))
  }
  as.integer(col)
}


################################################################################
# ---- 3. Print functions ----

#' Register helper files for `write_functions()` and `write_sections()`
#'
#' Indexes each file's top-level functions (by name, roxygen block included) and
#' its banner sections (`## ---- label ----`), recording which functions sit in
#' which section. This is the `knitr::read_chunk()` of this setup: it only reads
#' the files, so keep the usual `source()` calls next to it. Re-registering a
#' file refreshes its entries.
#'
#' @param ... Paths to the helper files.
#'
#' @returns Invisibly, the names available to the writers.
#' @export
register_helper_files <- function(...) {
  for (f in unlist(list(...), use.names = FALSE)) {
    lines <- readLines(f, warn = FALSE)
    exprs <- parse(text = lines, keep.source = TRUE)
    refs  <- attr(exprs, "srcref")
    hits  <- grep(.hf_rx, lines)
    labs  <- sub(.hf_rx, "\\1", lines[hits])

    # Drop any previous registration of this file
    .hf$funs     <- Filter(function(e) !identical(e$file, f), .hf$funs)
    .hf$sections <- Filter(function(e) !identical(e$file, f), .hf$sections)

    # Sections: header (rule line + banner) and body up to the next banner
    for (j in seq_along(hits)) {
      top <- if (hits[j] > 1L && grepl("^\\s*#{5,}\\s*$", lines[hits[j] - 1L])) hits[j] - 1L else hits[j]
      to  <- if (j < length(hits)) hits[j + 1L] - 1L else length(lines)
      .hf$sections[[labs[j]]] <- list(file = f, funs = character(0),
                                      header = lines[seq(top, hits[j])],
                                      body   = .hf_trim(lines[.hf_span(hits[j] + 1L, to)]))
    }

    # Functions: roxygen block included, tagged with the section they sit in
    for (i in seq_along(exprs)) {
      e <- exprs[[i]]
      if (!is.call(e) || length(e) < 3L || !as.character(e[[1L]])[1L] %in% c("<-", "=")) next
      if (!is.call(e[[3L]]) || !identical(as.character(e[[3L]][[1L]])[1L], "function")) next
      nm   <- as.character(e[[2L]])
      r    <- as.integer(refs[[i]])
      from <- r[1L]
      while (from > 1L && grepl("^\\s*#'", lines[from - 1L])) from <- from - 1L
      .hf$funs[[nm]] <- list(file = f, code = lines[seq(from, r[3L])])
      k <- sum(hits < from)
      if (k > 0L) .hf$sections[[labs[k]]]$funs <- c(.hf$sections[[labs[k]]]$funs, nm)
    }

    .hf$files <- unique(c(.hf$files, f))
  }
  invisible(hf_registered())
}

#' Write helper functions into a folded, non-evaluated code chunk
#'
#' @param labels Character vector of function names, or of section labels — a
#'   section expands to every function it contains, without its header. Section
#'   numbering, case and punctuation are ignored, so "2. Table Helpers" and
#'   "table helpers" are equivalent. Defaults to every registered function.
#' @param file Optional path or basename; restricts the default to one file.
#'
#' @returns Called for its output. Requires `results = 'asis'`.
#' @export
write_functions <- function(labels = NULL, file = NULL) {
  if (is.null(labels)) labels <- names(.hf_from(.hf$funs, file))
  .hf_emit(unlist(lapply(labels, .hf_fun_code), use.names = FALSE))
}

#' Write whole helper sections, headers included, into a folded code chunk
#'
#' @param labels Character vector of section labels. A function name is also
#'   accepted and writes that function alone. Defaults to every registered
#'   section.
#' @param file Optional path or basename; restricts the default to one file.
#' @param skip_empty Whether to drop sections holding no functions (library
#'   preambles, and so on). Applies only to the default set, never to labels
#'   you asked for by name.
#'
#' @returns Called for its output. Requires `results = 'asis'`.
#' @export
write_sections <- function(labels = NULL, file = NULL, skip_empty = TRUE) {
  if (is.null(labels)) {
    secs <- .hf_from(.hf$sections, file)
    if (skip_empty) secs <- Filter(function(s) length(s$funs) > 0L, secs)
    labels <- names(secs)
  }
  .hf_emit(unlist(lapply(labels, .hf_sec_code), use.names = FALSE))
}

#' Report what is registered, per file and section
#'
#' @returns Invisibly, a data frame of file / section / function.
#' @export
hf_status <- function() {
  if (!length(.hf$files)) {
    message("hf: nothing registered in this session.")
    return(invisible(data.frame()))
  }
  out <- data.frame(file = character(0), section = character(0), fun = character(0))
  for (f in .hf$files) {
    secs <- .hf_from(.hf$sections, f)
    orph <- setdiff(names(.hf_from(.hf$funs, f)), unlist(lapply(secs, `[[`, "funs"), use.names = FALSE))
    cat(basename(f), "\n", sep = "")
    for (nm in c(names(secs), if (length(orph)) "(no section)")) {
      fs <- if (identical(nm, "(no section)")) orph else secs[[nm]]$funs
      cat(sprintf("  %-34s %s\n", nm, if (length(fs)) paste(fs, collapse = ", ") else "-"))
      out <- rbind(out, data.frame(file = basename(f), section = nm,
                                   fun = if (length(fs)) fs else NA_character_))
    }
  }
  invisible(out)
}

#' Names currently available to the writers
#' @export
hf_registered <- function() sort(c(names(.hf$funs), names(.hf$sections)))

# Re-sourcing this file must not wipe an index built earlier in the session
if (!exists(".hf", inherits = FALSE) || !is.environment(.hf)) .hf <- new.env(parent = emptyenv())
if (is.null(.hf$funs))     .hf$funs     <- list()
if (is.null(.hf$sections)) .hf$sections <- list()
if (is.null(.hf$files))    .hf$files    <- character(0)

# Banner: "## ---- label ----" or "# ---- 2. Table Helpers ----". Four dashes,
# so ordinary "# --- note ---" comments are not mistaken for sections.
.hf_rx <- "^\\s*#+\\s*-{4,}\\s*(.+?)\\s*-{4,}\\s*$"

# Lookup key: drop case, section numbering and punctuation
.hf_key <- function(x) gsub("[^a-z0-9]", "", sub("^[0-9]+\\s*[.):-]?\\s*", "", tolower(trimws(x))))

.hf_span <- function(from, to) if (from > to) integer(0) else seq(from, to)

.hf_from <- function(index, file) {
  if (is.null(file)) return(index)
  Filter(function(e) identical(basename(e$file), basename(file)), index)
}

.hf_trim <- function(x) {
  keep <- which(nzchar(trimws(x)) & !grepl("^\\s*#+\\s*$", x))
  if (length(keep)) x[seq(min(keep), max(keep))] else character(0)
}

.hf_find <- function(label, index) {
  hit <- names(index)[.hf_key(names(index)) == .hf_key(label)]
  if (length(hit)) hit[1L] else NA_character_
}

.hf_fun_code <- function(label) {
  hit <- .hf_find(label, .hf$funs)
  if (!is.na(hit)) return(c(.hf$funs[[hit]]$code, ""))
  hit <- .hf_find(label, .hf$sections)
  if (is.na(hit)) stop("'", label, "' is not registered. See hf_status().", call. = FALSE)
  fs <- .hf$sections[[hit]]$funs
  if (!length(fs)) stop("Section '", label, "' holds no functions.", call. = FALSE)
  unlist(lapply(fs, function(n) c(.hf$funs[[n]]$code, "")), use.names = FALSE)
}

.hf_sec_code <- function(label) {
  hit <- .hf_find(label, .hf$sections)
  if (!is.na(hit)) return(c(.hf$sections[[hit]]$header, .hf$sections[[hit]]$body, ""))
  hit <- .hf_find(label, .hf$funs)
  if (is.na(hit)) stop("'", label, "' is not registered. See hf_status().", call. = FALSE)
  c(.hf$funs[[hit]]$code, "")
}

.hf_emit <- function(code) {
  while (length(code) && !nzchar(code[length(code)])) code <- code[-length(code)]
  if (!length(code)) return(invisible(NULL))

  # The fence must outrun any backticks in the code itself
  ticks <- nchar(unlist(regmatches(code, gregexpr("`+", code)), use.names = FALSE))
  fence <- strrep("`", max(3L, ticks + 1L))

  .hf$.n <- (if (is.null(.hf$.n)) 0L else .hf$.n) + 1L
  cat(knitr::knit_child(text = sprintf("%s{r hf-src-%d, eval = FALSE}\n%s\n%s",
                                       fence, .hf$.n, paste(code, collapse = "\n"), fence),
                        quiet = TRUE), sep = "\n")
}