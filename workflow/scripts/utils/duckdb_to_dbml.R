#!/usr/bin/env Rscript
## _________________________________________________
##
## DuckDB -> DBML (dbdiagram.io)
##
## Aim: Introspect a DuckDB file and emit the DBML that reproduces it as an
## entity-relationship diagram on https://dbdiagram.io/d
##
## Project: ENDome generation - Bin Information Content
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-09-01
##
## Latest Version: v1.0 (2026-09-01)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - Everything is read from the DuckDB catalog: duckdb_tables(),
##      duckdb_views(), duckdb_columns() and duckdb_constraints(). Only the
##      keys the database actually declares are drawn, so a join that is used
##      in the R code but never declared as a FOREIGN KEY will NOT appear.
##    - Tables that only ever exist as TEMPORARY (protein_pairs) or as an R
##      object (bin_pairs_distance), and any table that is not implemented yet
##      (runs), cannot be introspected. Keep those blocks in a sidecar file and
##      pass it with `--extra`; its content is appended verbatim.
##    - The connection is opened read-only, so it is safe to point at the
##      multi-GB DuckDBs under debug_results/ or Results/.
##
## Usage:
##    Rscript workflow/scripts/utils/duckdb_to_dbml.R <database.duckdb> [options]
##
##      -o, --output <file>   Write to <file> instead of stdout
##      --project <name>      Project name (default: the database basename)
##      --extra <file>        Append this file verbatim (hand-written blocks)
##      --exact-counts        COUNT(*) every table (default: catalog estimate)
##      --no-counts           Do not report row counts at all
##      --no-comments         Ignore COMMENT ON table/column metadata
##      --palette             Colour the table headers, cycling a palette
##      -h, --help            Print this message
##
## Examples:
##    Rscript workflow/scripts/utils/duckdb_to_dbml.R \
##      data/test_data/Wood.control/test.pc.w500.3p.duckdb \
##      --extra schema/ENDome_DB.extra.dbml -o schema/ENDome_DB.dbml
##
## Changelog:
##    - v1.0 (2026-09-01): Initial version
##
## Contact: guillermorocamora@gmail.com
## _________________________________________________

############################################################################## #
# ---- 0. Setup ----

suppressWarnings(suppressMessages({
  library(DBI)
  library(duckdb)
}))

#----------------------------------------------------------------------------- #
## 0.1 Command line arguments ----
usage <- function() {
  cat(
    "Usage: duckdb_to_dbml.R <database.duckdb> [options]\n\n",
    "  -o, --output <file>   Write to <file> instead of stdout\n",
    "      --project <name>  Project name (default: the database basename)\n",
    "      --extra <file>    Append this file verbatim (hand-written blocks)\n",
    "      --exact-counts    COUNT(*) every table (default: catalog estimate)\n",
    "      --no-counts       Do not report row counts at all\n",
    "      --no-comments     Ignore COMMENT ON table/column metadata\n",
    "      --palette         Colour the table headers, cycling a palette\n",
    "  -h, --help            Print this message\n",
    sep = ""
  )
}

parseArgs <- function(args) {
  opts <- list(
    db = NA_character_, output = NA_character_, project = NA_character_,
    extra = NA_character_, exact_counts = FALSE, counts = TRUE,
    comments = TRUE, palette = FALSE
  )

  # Options taking a value, mapped to their slot in `opts`
  valued <- c("-o" = "output", "--output" = "output", "--project" = "project",
              "--extra" = "extra")

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0)
    } else if (arg %in% names(valued)) {
      if (i == length(args)) stop("Missing value for ", arg, call. = FALSE)
      opts[[valued[[arg]]]] <- args[[i + 1]]
      i <- i + 2
      next
    } else if (arg == "--exact-counts") {
      opts$exact_counts <- TRUE
    } else if (arg == "--no-counts") {
      opts$counts <- FALSE
    } else if (arg == "--no-comments") {
      opts$comments <- FALSE
    } else if (arg == "--palette") {
      opts$palette <- TRUE
    } else if (startsWith(arg, "-")) {
      stop("Unknown option: ", arg, call. = FALSE)
    } else if (is.na(opts$db)) {
      opts$db <- arg
    } else {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
    i <- i + 1
  }

  if (is.na(opts$db)) {
    usage()
    stop("No DuckDB file given.", call. = FALSE)
  }
  if (!file.exists(opts$db)) stop("No such file: ", opts$db, call. = FALSE)
  if (!is.na(opts$extra) && !file.exists(opts$extra)) {
    stop("No such file: ", opts$extra, call. = FALSE)
  }
  if (is.na(opts$project)) {
    opts$project <- sub("\\.duckdb$", "", basename(opts$db))
  }
  opts
}

############################################################################## #
# ---- 1. DBML emitters ----

### Anything that is not a plain identifier, or that could be read as a DBML
### keyword, has to be quoted. Over-quoting is harmless, under-quoting is not.
dbml_keywords <- c(
  "table", "tablegroup", "ref", "enum", "project", "note", "indexes", "index",
  "pk", "primary", "key", "unique", "increment", "default", "null", "not",
  "as", "type", "name", "color", "headercolor", "start", "end", "in", "out"
)

#' Quote a table/column identifier if DBML needs it
quoteId <- function(x) {
  needs_quote <- !grepl("^[A-Za-z_][A-Za-z0-9_]*$", x) | tolower(x) %in% dbml_keywords
  ifelse(needs_quote, sprintf('"%s"', gsub('"', "", x, fixed = TRUE)), x)
}

#' Format a string as a DBML note body
#'
#' Single-quoted for one-liners; triple-quoted when the text is multi-line or
#' contains a quote of its own.
noteLiteral <- function(x) {
  if (grepl("\n", x, fixed = TRUE) || grepl("'", x, fixed = TRUE)) {
    # A note must not end on a quote, or it would close the delimiter early
    if (endsWith(x, "'")) x <- paste0(x, " ")
    sprintf("'''%s'''", x)
  } else {
    sprintf("'%s'", x)
  }
}

#' Turn the DuckDB types into DBML types, lifting ENUMs into Enum blocks
#'
#' DuckDB reports an anonymous ENUM as `ENUM('chr1', 'chr2', ...)`, which the
#' DBML parser cannot read inline. Every distinct value set becomes one shared
#' `Enum` block named after the first column that uses it.
#'
#' @param columns Data frame from \code{getColumns()}.
#'
#' @returns List with \code{columns} (a \code{dbml_type} column added) and
#'   \code{enums} (named list of value vectors).
resolveTypes <- function(columns) {
  enums <- list()
  enum_keys <- character(0)

  dbml_type <- vapply(seq_len(nrow(columns)), function(i) {
    type <- columns$data_type[[i]]

    ### ENUM('a', 'b', ...) -> a named Enum block
    if (grepl("^ENUM\\(", type)) {
      values <- regmatches(type, gregexpr("'(?:[^']|'')*'", type))[[1]]
      values <- gsub("''", "'", substr(values, 2, nchar(values) - 1), fixed = TRUE)
      key <- paste(values, collapse = "\x1f")

      known <- match(key, enum_keys)
      if (is.na(known)) {
        enum_name <- sprintf("%s_%s", columns$table_name[[i]], columns$column_name[[i]])
        enums[[enum_name]] <<- values
        enum_keys <<- c(enum_keys, key)
        names(enum_keys)[length(enum_keys)] <<- enum_name
        return(enum_name)
      }
      return(names(enum_keys)[known])
    }

    ### Plain types read better lowercased; anything else is quoted verbatim so
    ### that brackets, commas and nested types cannot break the parser
    if (grepl("^[A-Za-z0-9_ ]+(\\([0-9]+(,[0-9]+)?\\))?$", type)) tolower(type)
    else sprintf('"%s"', gsub('"', "", type, fixed = TRUE))
  }, character(1))

  columns$dbml_type <- dbml_type
  list(columns = columns, enums = enums)
}

#' Render the Enum blocks
renderEnums <- function(enums) {
  if (length(enums) == 0) return(character(0))

  unlist(lapply(names(enums), function(nm) {
    values <- enums[[nm]]
    simple <- grepl("^[A-Za-z_][A-Za-z0-9_]*$", values)
    c(sprintf("Enum %s {", quoteId(nm)),
      sprintf("  %s", ifelse(simple, values, sprintf('"%s"', values))),
      "}", "")
  }))
}

#' Collapse the column settings into the `[...]` suffix
settingsSuffix <- function(settings) {
  settings <- settings[!is.na(settings) & nzchar(settings)]
  if (length(settings) == 0) return("")
  sprintf(" [%s]", paste(settings, collapse = ", "))
}

############################################################################## #
# ---- 2. Catalog introspection ----

#' Every non-internal table and view of the main schema
getRelations <- function(con, with_comments = TRUE) {
  tables <- DBI::dbGetQuery(con, "
    SELECT table_name AS name, 'table' AS kind, comment, estimated_size, temporary
    FROM duckdb_tables()
    WHERE NOT internal AND schema_name = current_schema()")

  views <- DBI::dbGetQuery(con, "
    SELECT view_name AS name, 'view' AS kind, comment, NULL AS estimated_size, temporary
    FROM duckdb_views()
    WHERE NOT internal AND schema_name = current_schema()")

  rel <- rbind(tables, views)
  if (!with_comments) rel$comment <- NA_character_
  rel[order(rel$kind, rel$name), , drop = FALSE]
}

#' Every column of the main schema, in declaration order
getColumns <- function(con, with_comments = TRUE) {
  cols <- DBI::dbGetQuery(con, "
    SELECT table_name, column_name, column_index, data_type, is_nullable,
           column_default, comment
    FROM duckdb_columns()
    WHERE NOT internal AND schema_name = current_schema()
    ORDER BY table_name, column_index")
  if (!with_comments) cols$comment <- NA_character_
  cols
}

#' Every constraint, with the column lists flattened to comma-separated strings
getConstraints <- function(con) {
  DBI::dbGetQuery(con, "
    SELECT table_name,
           constraint_type,
           array_to_string(constraint_column_names, ',') AS cols,
           referenced_table,
           array_to_string(referenced_column_names, ',')  AS ref_cols,
           constraint_text
    FROM duckdb_constraints()
    WHERE schema_name = current_schema()
    ORDER BY table_name, constraint_type")
}

#' Exact or estimated row count per relation
getRowCounts <- function(con, relations, mode = c("estimate", "exact", "none")) {
  mode <- match.arg(mode)
  if (mode == "none") return(stats::setNames(rep(NA_real_, nrow(relations)), relations$name))
  if (mode == "estimate") return(stats::setNames(as.numeric(relations$estimated_size), relations$name))

  counts <- vapply(relations$name, function(nm) {
    as.numeric(DBI::dbGetQuery(con, sprintf('SELECT count(*) AS n FROM "%s"', nm))$n)
  }, numeric(1))
  stats::setNames(counts, relations$name)
}

############################################################################## #
# ---- 3. DBML generation ----

#' Render one Table block
renderTable <- function(name, cols, constraints, row_count, kind, comment,
                        header_color = NA_character_, counts_mode = "estimate") {
  tbl_constraints <- constraints[constraints$table_name == name, , drop = FALSE]

  splitCols <- function(x) strsplit(x, ",", fixed = TRUE)[[1]]

  pk <- tbl_constraints[tbl_constraints$constraint_type == "PRIMARY KEY", ]
  pk_cols <- if (nrow(pk) > 0) splitCols(pk$cols[[1]]) else character(0)

  ### Only single-column UNIQUE / NOT NULL can be inlined on the column
  uniques <- tbl_constraints[tbl_constraints$constraint_type == "UNIQUE", , drop = FALSE]
  unique_lists <- lapply(uniques$cols, splitCols)
  unique_single <- unlist(unique_lists[lengths(unique_lists) == 1])

  not_null <- tbl_constraints[tbl_constraints$constraint_type == "NOT NULL", , drop = FALSE]
  not_null_cols <- unlist(lapply(not_null$cols, splitCols))

  ### Foreign keys are drawn as `Ref:` lines, but flagged here so the diagram
  ### reads correctly even before the refs are resolved
  fks <- tbl_constraints[tbl_constraints$constraint_type == "FOREIGN KEY", , drop = FALSE]

  header_setting <- if (!is.na(header_color)) sprintf(" [headercolor: %s]", header_color) else ""
  lines <- sprintf("Table %s%s {", quoteId(name), header_setting)

  for (i in seq_len(nrow(cols))) {
    col <- cols[i, ]
    settings <- character(0)

    if (length(pk_cols) == 1 && col$column_name == pk_cols) settings <- c(settings, "pk")
    if (col$column_name %in% unique_single) settings <- c(settings, "unique")
    ### A PK is implicitly NOT NULL; do not repeat it
    if (col$column_name %in% not_null_cols && !col$column_name %in% pk_cols) {
      settings <- c(settings, "not null")
    }
    if (!is.na(col$column_default) && nzchar(col$column_default)) {
      settings <- c(settings, sprintf("default: `%s`", col$column_default))
    }
    if (!is.na(col$comment) && nzchar(col$comment)) {
      settings <- c(settings, sprintf("note: %s", noteLiteral(col$comment)))
    }

    lines <- c(lines, sprintf("  %s %s%s",
      quoteId(col$column_name), col$dbml_type, settingsSuffix(settings)))
  }

  ### Composite primary keys and every unique constraint go to the index block
  index_lines <- character(0)
  if (length(pk_cols) > 1) {
    index_lines <- c(index_lines, sprintf("    (%s) [pk]",
      paste(quoteId(pk_cols), collapse = ", ")))
  }
  for (u in unique_lists[lengths(unique_lists) > 1]) {
    index_lines <- c(index_lines, sprintf("    (%s) [unique]",
      paste(quoteId(u), collapse = ", ")))
  }
  if (length(index_lines) > 0) {
    lines <- c(lines, "", "  indexes {", index_lines, "  }")
  }

  ### Table note: provenance the catalog knows about, plus any COMMENT ON
  note_bits <- character(0)
  if (kind == "view") note_bits <- c(note_bits, "DuckDB VIEW.")
  if (!is.na(row_count)) {
    note_bits <- c(note_bits, sprintf("%s rows%s.",
      format(row_count, big.mark = ",", scientific = FALSE, trim = TRUE),
      if (counts_mode == "estimate") " (catalog estimate)" else ""))
  }
  if (nrow(fks) > 0) {
    note_bits <- c(note_bits, sprintf("Foreign keys: %s.",
      paste(sprintf("%s -> %s.%s", fks$cols, fks$referenced_table, fks$ref_cols),
            collapse = "; ")))
  }
  if (!is.na(comment) && nzchar(comment)) note_bits <- c(note_bits, comment)

  if (length(note_bits) > 0) {
    lines <- c(lines, "", sprintf("  Note: %s", noteLiteral(paste(note_bits, collapse = " "))))
  }

  c(lines, "}", "")
}

#' Render every `Ref:` line
#'
#' A foreign key whose child columns are exactly the child's primary key is a
#' one-to-one relationship (`-`); anything else is many-to-one (`>`).
renderRefs <- function(constraints) {
  fks <- constraints[constraints$constraint_type == "FOREIGN KEY", , drop = FALSE]
  if (nrow(fks) == 0) return(character(0))

  pks <- constraints[constraints$constraint_type == "PRIMARY KEY", , drop = FALSE]
  pk_of <- stats::setNames(as.list(pks$cols), pks$table_name)

  colRef <- function(table, cols) {
    cols <- strsplit(cols, ",", fixed = TRUE)[[1]]
    if (length(cols) == 1) {
      sprintf("%s.%s", quoteId(table), quoteId(cols))
    } else {
      sprintf("%s.(%s)", quoteId(table), paste(quoteId(cols), collapse = ", "))
    }
  }

  vapply(seq_len(nrow(fks)), function(i) {
    fk <- fks[i, ]
    child_pk <- pk_of[[fk$table_name]]
    cardinality <- if (!is.null(child_pk) && identical(child_pk, fk$cols)) "-" else ">"
    self_ref <- if (identical(fk$table_name, fk$referenced_table)) "  // self-reference" else ""

    sprintf("Ref: %s %s %s%s",
      colRef(fk$table_name, fk$cols), cardinality,
      colRef(fk$referenced_table, fk$ref_cols), self_ref)
  }, character(1))
}

#' Order the tables so that the referenced ones come first
orderRelations <- function(relations, constraints) {
  fks <- constraints[constraints$constraint_type == "FOREIGN KEY", , drop = FALSE]
  ### Self-references do not make a table depend on anything else
  fks <- fks[fks$table_name != fks$referenced_table, , drop = FALSE]

  n_parents <- vapply(relations$name, function(nm) sum(fks$table_name == nm), integer(1))
  relations[order(n_parents, relations$name), , drop = FALSE]
}

############################################################################## #
# ---- 4. Main ----
opts <- parseArgs(commandArgs(trailingOnly = TRUE))

message("Reading the catalog of ", opts$db, " ...")
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = opts$db, read_only = TRUE)
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

relations <- getRelations(con, with_comments = opts$comments)
if (nrow(relations) == 0) stop("No user tables found in ", opts$db, call. = FALSE)

columns <- getColumns(con, with_comments = opts$comments)
constraints <- getConstraints(con)

resolved <- resolveTypes(columns)
columns <- resolved$columns

counts_mode <- if (!opts$counts) "none" else if (opts$exact_counts) "exact" else "estimate"
row_counts <- getRowCounts(con, relations, counts_mode)

relations <- orderRelations(relations, constraints)

### Colours are cosmetic: cycled over the tables in the order they are drawn
palette <- c("#1F77B4", "#2CA02C", "#FF7F0E", "#9467BD", "#C0392B", "#17BECF")
header_colors <- if (opts$palette) {
  stats::setNames(rep_len(palette, nrow(relations)), relations$name)
} else {
  stats::setNames(rep(NA_character_, nrow(relations)), relations$name)
}

#----------------------------------------------------------------------------- #
## 4.1 Assemble the document ----
duckdb_version <- DBI::dbGetQuery(con, "SELECT version() AS v")$v

header <- c(
  "// ============================================================================",
  sprintf("// %s", opts$project),
  "//",
  "// GENERATED by workflow/scripts/utils/duckdb_to_dbml.R - do not edit by hand.",
  sprintf("//   source     : %s", normalizePath(opts$db)),
  sprintf("//   generated  : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("//   duckdb     : %s", duckdb_version),
  if (!is.na(opts$extra)) sprintf("//   extra block: %s", opts$extra),
  "//",
  "// Paste into https://dbdiagram.io/d",
  "// ============================================================================",
  "",
  sprintf("Project %s {", quoteId(opts$project)),
  "  database_type: 'DuckDB'",
  sprintf("  Note: %s", noteLiteral(sprintf(
    "Introspected from %s. Only the keys DuckDB actually declares are drawn.",
    basename(opts$db)))),
  "}",
  ""
)

enum_block <- if (length(resolved$enums) > 0) {
  c("// ============================================================================",
    "// Enumerated types",
    "// ============================================================================",
    "",
    renderEnums(resolved$enums))
} else character(0)

table_blocks <- unlist(lapply(seq_len(nrow(relations)), function(i) {
  rel <- relations[i, ]
  renderTable(
    name = rel$name,
    cols = columns[columns$table_name == rel$name, , drop = FALSE],
    constraints = constraints,
    row_count = row_counts[[rel$name]],
    kind = rel$kind,
    comment = rel$comment,
    header_color = header_colors[[rel$name]],
    counts_mode = counts_mode
  )
}))

ref_lines <- renderRefs(constraints)
ref_block <- c(
  "// ============================================================================",
  "// Relationships (declared FOREIGN KEYs)",
  "// ============================================================================",
  "",
  if (length(ref_lines) > 0) ref_lines else "// (none declared)",
  ""
)

extra_block <- if (!is.na(opts$extra)) {
  c("// ============================================================================",
    sprintf("// Appended verbatim from %s", opts$extra),
    "// ============================================================================",
    "",
    readLines(opts$extra, warn = FALSE),
    "")
} else character(0)

dbml <- c(header, enum_block, table_blocks, ref_block, extra_block)

#----------------------------------------------------------------------------- #
## 4.2 Write it out ----
if (is.na(opts$output)) {
  cat(dbml, sep = "\n")
} else {
  dir.create(dirname(opts$output), recursive = TRUE, showWarnings = FALSE)
  writeLines(dbml, opts$output)
  message("Wrote ", opts$output)
}

message(sprintf("%i table(s), %i view(s), %i column(s), %i foreign key(s)",
  sum(relations$kind == "table"), sum(relations$kind == "view"), nrow(columns),
  sum(constraints$constraint_type == "FOREIGN KEY")))
