## _________________________________________________
##
## Helper Functions: DuckDB
##
## Aim: Include into a single file the functions shared across the "Bin
## Information Content" scripts related to the DuckDB creation, extraction
## and calculations.
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 18/08/2026
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.1 (19/08/2026)
## _________________________________________________
##
## Notes:
##     - VENDORED COPY of HelperFunctions/Phf-ENDome_Pipeline/Phf-DuckDB.R.
##       The pipeline sources its helpers from workflow/scripts/lib/, since
##       reaching across the `R -> ../HelperFunctions` symlink does not survive
##       a `--use-conda` run. Keep the two in sync, or retire the original.
##
## Changelog:
##     - v1.1 (19/08/2026): Vendored for the Snakemake port. Added contentId()
##       and computeRowHash() (the latter moved here from T01, where it was
##       defined but never called).
##     - v1.0 (18/08/2026): Migrated DuckDB functions from Bin_Information_Score.R helper
##
## Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

################################################################################
# ---- 1. Required libraries ----
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(duckdb)))
suppressWarnings(suppressMessages(library(DBI)))

################################################################################
# ---- 2. Table Helpers ----

#' List every table of the transcript database as lazy tibbles
#'
#' Convenience accessor that returns all the tables found in a DuckDB
#' connection as a named list of lazy \code{dplyr} tibbles. It allows referring
#' to any table as \code{db$transcripts}, \code{db$bins}, etc. without calling
#' \code{dplyr::tbl()} every time.
#'
#' @param con DBIConnection to the transcript DuckDB.
#'
#' @returns Named list of lazy tibbles, one per table in the connection.
#' @export
getTbls <- function(con) {
  tbl_names <- DBI::dbListTables(con)
  stats::setNames(lapply(tbl_names, function(tbl_name) dplyr::tbl(con, tbl_name)), tbl_names)
}

#' Primary keys defined in the transcript database
#'
#' @param con DBIConnection to the transcript DuckDB.
#'
#' @returns Dataframe with columns \code{table} and \code{pk_col}, the latter a
#'   comma-separated string of the columns forming the primary key.
#' @export
getDbPrimaryKeys <- function(con) {
  DBI::dbGetQuery(con, "
    SELECT table_name AS \"table\",
           array_to_string(constraint_column_names, ', ') AS pk_col
    FROM duckdb_constraints()
    WHERE constraint_type = 'PRIMARY KEY' AND schema_name = 'main'
    ORDER BY table_name")
}

#' Foreign keys defined in the transcript database
#'
#' @param con DBIConnection to the transcript DuckDB.
#'
#' @returns Dataframe with columns \code{child_table}, \code{child_fk_cols},
#'   \code{parent_table} and \code{parent_key_cols}.
#' @export
getDbForeignKeys <- function(con) {
  DBI::dbGetQuery(con, "
    SELECT table_name AS child_table,
           array_to_string(constraint_column_names, ', ') AS child_fk_cols,
           referenced_table                               AS parent_table,
           array_to_string(referenced_column_names, ', ') AS parent_key_cols
    FROM duckdb_constraints()
    WHERE constraint_type = 'FOREIGN KEY' AND schema_name = 'main'
    ORDER BY child_table, child_fk_cols")
}

################################################################################
# ---- 3. Content-derived identifiers ----

#' Stable content hash per row
#'
#' @param df A data.frame.
#' @param id_cols Character vector of column names whose combined values
#'   define row identity.
#' @param algo \code{digest::digest()} algorithm (default "xxhash64").
#'
#' @returns Character vector of hashes, one per row of \code{df}.
#' @export
computeRowHash <- function(df, id_cols, algo = "xxhash64"){
  stopifnot(all(id_cols %in% names(df)))

  # Add ASCII unit-separator (0x1F) between concatenates columns to ensure
  # different hasing
  content <- do.call(paste, c(as.list(df[id_cols]), sep = "\x1f"))
  vapply(content, digest::digest, character(1), algo = algo, serialize = FALSE, USE.NAMES = FALSE)
}

#' Assign a content-derived surrogate key
#'
#' Builds an identifier of the form \code{<prefix>_<hash>} from the columns that
#' define a row's identity.
#'
#' @param df A data.frame.
#' @param id_cols Character vector of the columns defining row identity.
#' @param prefix String prepended to the hash.
#' @param algo \code{digest::digest()} algorithm (default "xxhash64").
#'
#' @returns Character vector of ids, one per row of \code{df}.
#' @export
assignId <- function(df, id_cols, prefix, algo = "xxhash64"){
  sprintf("%s_%s", prefix, computeRowHash(df, id_cols, algo = algo))
}

################################################################################
# ---- 4. Table generation ----

#' Write a data frame to DuckDB with its primary and foreign keys
#'
#' @param con DBIConnection to the DuckDB.
#' @param name String, name of the table to create.
#' @param spec List with \code{df} (the data), and optionally \code{pk} (primary
#'   key column(s)) and \code{fk} (named vector/list mapping a column of
#'   \code{df} to a "table.column" reference).
#'
#' @returns Invisibly, the number of rows written.
#' @export
writeDuckTable <- function(con, name, spec) {
  pk_clause <- if (!is.null(spec$pk)) paste(sprintf('"%s"', spec$pk), collapse = ", ")

  # If there's no FK, simply create the table
  if (length(spec$fk) == 0) {
    DBI::dbWriteTable(con, name, spec$df, overwrite = TRUE)
    if (!is.null(spec$pk)) {
      DBI::dbExecute(con, sprintf('ALTER TABLE "%s" ADD PRIMARY KEY (%s)', name, pk_clause))
    }
    return(invisible(nrow(spec$df)))
  }

  # Add FK via staging table
  staging <- paste0("staging_", name)
  DBI::dbWriteTable(con, staging, spec$df, overwrite = TRUE)
  on.exit(DBI::dbExecute(con, sprintf('DROP TABLE IF EXISTS "%s"', staging)), add = TRUE)

  col_info <- DBI::dbGetQuery(con, sprintf('DESCRIBE "%s"', staging))
  col_defs <- sprintf('"%s" %s', col_info$column_name, col_info$column_type)

  pk_def <- if (!is.null(spec$pk)) sprintf("PRIMARY KEY (%s)", pk_clause)

  # `spec$fk` may be a named character vector or a named list; unlist() so that
  # str_split_i() sees a character vector either way.
  fk <- unlist(spec$fk)
  ref_table <- stringr::str_split_i(fk, stringr::fixed("."), 1)
  ref_cols  <- stringr::str_split_i(fk, stringr::fixed("."), 2)

  fk_defs <- sprintf('FOREIGN KEY ("%s") REFERENCES "%s" ("%s")', names(fk), ref_table, ref_cols)
  all_defs <- paste(c(col_defs, pk_def, fk_defs), collapse = ",\n  ")

  DBI::dbExecute(con, sprintf('DROP TABLE IF EXISTS "%s"', name))
  DBI::dbExecute(con, sprintf('CREATE TABLE "%s" (\n  %s\n)', name, all_defs))

  # Self-referential FKs (e.g. transcripts.superseded_by -> transcripts.transcript_id)
  # require the referenced rows to exist first. `resolveTerminalNodes()` maps every
  # transcript directly to its terminal node, so a single "targets first, referrers
  # second" pass is enough -- there is never a chain to walk.
  self_cols <- names(fk)[ref_table == name]
  if (length(self_cols) > 0) {
    is_target <- paste(sprintf('"%s" IS NULL', self_cols), collapse = " AND ")
    DBI::dbExecute(con, sprintf('INSERT INTO "%s" SELECT * FROM "%s" WHERE %s', name, staging, is_target))
    DBI::dbExecute(con, sprintf('INSERT INTO "%s" SELECT * FROM "%s" WHERE NOT (%s)', name, staging, is_target))
  } else {
    DBI::dbExecute(con, sprintf('INSERT INTO "%s" SELECT * FROM "%s"', name, staging))
  }

  invisible(nrow(spec$df))
}

################################################################################
# ---- 5. Others ----

#' Parse the ORFannotate FASTA headers into a tibble
#'
#' Headers look like `<transcript_id> gene_id=<...>;coding_prob=<...>`.
#'
#' @param string_set XStringSet read from an ORFannotate FASTA.
#'
#' @returns Tibble with transcript_id, gene_id, coding_prob, width, sequence.
StringSetToTibble <- function(string_set){
  data.frame(
    header = names(string_set),
    sequence = as.character(string_set),
    width = BiocGenerics::width(string_set)
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      transcript_id = str_extract(header, "^[^ ]+"),
      gene_id = str_extract(header, "gene_id=([^;]+)") %>% str_remove("gene_id="),
      coding_prob = str_extract(header, "coding_prob=([0-9.]+)") %>% str_remove("coding_prob=") %>% as.numeric()
    ) %>%
    dplyr::select(transcript_id, gene_id, coding_prob, width, sequence)
}