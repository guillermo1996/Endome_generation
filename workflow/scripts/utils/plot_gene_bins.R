## _________________________________________________
##
## Bin structure plots
##
## Aim: Draw any set of ENDome bins - the coding structure of every member over
## the region kept by the truncation - annotated with the bin size and its
## information score
##
## Project: ENDome generation - Bin Information Content
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Contributors:
##
## Date Created: 2026-09-01
##
## Latest Version: v3.0 (2026-09-01)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - Section 1 opens the DuckDB and leaves `db_bins`, `db_transcripts` and
##      `db_bin_information` in the session as lazy `dplyr` tables. Query them
##      like any other data frame, `collect()` what you want, and hand the
##      `bin_id`s to `plotBins()`. Nothing else is needed - no gene, no config.
##    - The bins are numbered `Bin #1`, `Bin #2`, ... in the order they are
##      passed in, so the ordering is whatever the query arranged them by. They
##      do not have to belong to the same gene; the strip names the gene.
##    - Every transcript is drawn as two half-height tracks. The TOP one is its
##      coding structure - 5' UTR, CDS and 3' UTR, one colour each - and the
##      BOTTOM one is the region the truncation keeps. Reading a row vertically
##      therefore says which part of the ORF survives the truncation, which is
##      what the information score of the bin is computed over.
##    - Each feature is drawn at its own height (`feature_heights`), with the CDS
##      standing taller than the UTRs and the truncated region.
##    - Give two features the same entry in `feature_labels` AND the same entry
##      in `feature_colours` and they share one legend key - labelling both UTRs
##      `"UTR"` in one grey draws a single `UTR` key. Same label, different
##      colours keeps a key each, since one key could only describe one of them.
##    - The exon and CDS structure comes from the GTF, which has to be the
##      PRE-truncation ORF-filtered GTF the DuckDB was built from
##      (`04-ORF_Identification/ORF_Filtration/*.gtf`). The UTRs are not stored
##      exon-by-exon anywhere, so they are derived: the CDS is contiguous in
##      transcript space, so every exonic base below the first CDS base and above
##      the last one is UTR, and the strand says which of the two it is.
##    - The truncated (i.e. retained) region is NOT read from the txendcutr GTF.
##      It is recomputed: the terminal `width` transcript bases of the bin
##      REPRESENTATIVE - the transcript whose id names the bin - are walked back
##      to a genomic boundary, and every member of the bin is tagged against that
##      same boundary. This mirrors what txendcutr collapses the bin on, and it
##      is what makes the members of a bin look alike in the truncated end.
##      `width` and the truncated end are read off `run_id`.
##    - The DuckDB is opened read-only, so it is safe to point at the multi-GB
##      files under debug_results/ or Results/.
##    - This script was documented with the assistance of generative AI
##      (Claude Opus 5).
##
## Changelog:
##    - v3.1 (2026-09-01): per-feature heights.
##    - v3.0 (2026-09-01): drive everything off a plain query against the DuckDB
##      and one `plotBins()` call; dropped the per-gene configuration.
##    - v2.1 (2026-09-01): split every transcript into a coding-structure track
##      and a truncation track.
##    - v2.0 (2026-09-01): reworked as an interactive script.
##    - v1.0 (2026-09-01): Initial version, command line only.
##
## Contact: guillermorocamora@gmail.com
## _________________________________________________

############################################################################## #
# ---- 0. Setup ----

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages # Shortcut to hide package start-up messages

if(!requireNamespace("ggtranscript", quietly = TRUE)){
  devtools::install_github("dzhang32/ggtranscript", upgrade = "never")
}

shhh({
  library(DBI)
  library(duckdb)
  library(rtracklayer)
  library(ggtranscript)
  library(tidyverse)
})

### Package options - set package-level behavior for the session
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----

### 0.2.1 User-defined Parameters ----

### Name of the figure this run writes
figure_label <- "top_information_score"

### Figure aesthetics
figure_width <- 10 # inches; the height is scaled to the content

### Height of every feature drawn. This is the `height` aesthetic of
### `geom_half_range()`, which draws half of it on its side of the transcript,
### so the band is half as thick as the number here. The CDS stands taller than
### the UTRs and the truncated region so the coding core reads at a glance.
feature_heights <- c("5' UTR" = 0.15,
                     "CDS" = 0.25,
                     "3' UTR" = 0.15,
                     "Truncated" = 0.15,
                     "Other" = 0.15)

### Okabe-Ito throughout: blues and green for the coding structure, orange for
### the region the truncation keeps, neutral grey for anything unassigned
feature_colours <- c("5' UTR" = "#bebebe",
                     "CDS" = "#ac39ff",
                     "3' UTR" = "#bebebe",
                     "Truncated" = "#E69F00",
                     "Other" = "#BEBEBE")
feature_labels <- c("5' UTR" = "UTR",
                    "CDS" = "CDS",
                    "3' UTR" = "UTR",
                    "Truncated" = "Truncated Regions",
                    "Other" = "Other")

### 0.2.3 Parameter Checks and Derived Parameters ----

#### Checks
stopifnot(is.character(figure_label), length(figure_label) == 1)
stopifnot(all(c("5' UTR", "CDS", "3' UTR", "Truncated", "Other") %in% names(feature_colours)))
stopifnot(setequal(names(feature_colours), names(feature_labels)))
stopifnot(setequal(names(feature_colours), names(feature_heights)))
stopifnot(all(feature_heights > 0))

#----------------------------------------------------------------------------- #
## 0.3 Script Paths ----
main_path <- "/home/MRGuillermoPerez/RytenLab-Research/40-ENDome_generation"

### The pipeline run to inspect: one iteration of one prefix
analysis_path <- file.path(main_path, "debug_results/Wood.control")

### Input Paths
duckdb_path <- file.path(analysis_path, "05-Truncation-03688/DuckDB/Wood.control.stringtie.pc.w500.3p.duckdb")
gtf_path <- file.path(analysis_path, "04-ORF_Identification-d8708/ORF_Filtration/Wood.control.stringtie.pc.orf_filter.gtf")

### Output Paths
results_path <- file.path(analysis_path, "Gene_bins/Results")
figure_path <- file.path(analysis_path, "Gene_bins/Figures")

bins_fig_path <- \(label, ext = "png") file.path(figure_path, sprintf("%s_bins.%s", label, ext))
bin_summary_output_path <- \(label) file.path(results_path, sprintf("%s_bins.tsv", label))

#### Create output directories
dir.create(results_path, showWarnings = F, recursive = T)
dir.create(figure_path, showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.4 Helper Functions ----

#' Save a plot as both PNG and PDF at the same size/resolution
saveBothExt <- function(fig_path, plot, label, width, height){
  ggsave(filename = fig_path(label, "png"), plot = plot, width = width, height = height,
         units = "in", dpi = 300, limitsize = FALSE)
  ggsave(filename = fig_path(label, "pdf"), plot = plot, width = width, height = height,
         units = "in", dpi = 300, device = cairo_pdf, limitsize = FALSE)
}

#' The run this DuckDB describes, and the truncation it was built with
#'
#' Every table carries the same scalar `run_id`, and `run_id` is
#' `{prefix}.{orf_filter}.w{width}.{txEnd}`, so the last two dot fields carry
#' everything needed to recompute the truncation boundary.
#'
#' @param con Open DuckDB connection.
#'
#' @return A list with `run_id`, `trunc_width` and `tx_end`.
runSettings <- function(con){
  run_id <- DBI::dbGetQuery(con, "SELECT DISTINCT run_id FROM bins")$run_id
  stopifnot(length(run_id) == 1)

  m <- stringr::str_match(run_id, "\\.w([0-9]+)\\.(3p|5p)$")
  if(is.na(m[1, 1])) stop("Cannot read the truncation off run_id '", run_id, "'.", call. = FALSE)

  list(run_id = run_id, trunc_width = as.integer(m[1, 2]), tx_end = m[1, 3])
}

#' Every member of a set of bins, with the bin size and the information score
#'
#' @param con Open DuckDB connection.
#' @param bin_ids Bins to pull.
#'
#' @return A `tibble::tibble` with one row per transcript.
queryBinMembers <- function(con, bin_ids){
  dplyr::tbl(con, "transcripts") %>%
    dplyr::filter(bin_id %in% !!bin_ids) %>%
    dplyr::select(transcript_id, bin_id, gene_id, gene_name, strand, superseded_by) %>%
    dplyr::left_join(dplyr::tbl(con, "bins") %>%
                       dplyr::select(bin_id, n_members, n_superseded), by = "bin_id") %>%
    dplyr::left_join(dplyr::tbl(con, "bin_information") %>%
                       dplyr::select(bin_id, info_score_norm, info_score), by = "bin_id") %>%
    dplyr::collect() %>%
    dplyr::mutate(is_rep = transcript_id == bin_id,
                  is_superseded = !is.na(superseded_by))
}

#' Read the exons of a set of transcripts out of the annotation
#'
#' @param gtf_gr `GRanges` of the full annotation.
#' @param transcript_ids Transcripts to keep.
#'
#' @return A `tibble::tibble` with `tx_id`, `seqnames`, `start`, `end`, `strand`
#'   and `type`.
extractExons <- function(gtf_gr, transcript_ids){
  tibble::as_tibble(gtf_gr) %>%
    dplyr::filter(type == "exon", transcript_id %in% transcript_ids) %>%
    dplyr::transmute(tx_id = transcript_id, seqnames = as.character(seqnames),
                     start = as.integer(start), end = as.integer(end),
                     strand = as.character(strand), type = "exon")
}

#' First and last coding base of a set of transcripts
#'
#' The CDS is contiguous in transcript space, so its genomic span is all that is
#' needed to tell the coding part of an exon from the two UTRs around it.
#'
#' @param gtf_gr `GRanges` of the full annotation.
#' @param transcript_ids Transcripts to keep.
#'
#' @return A `tibble::tibble` with `tx_id`, `cds_min` and `cds_max`. Transcripts
#'   without a CDS are simply absent.
extractCdsBounds <- function(gtf_gr, transcript_ids){
  tibble::as_tibble(gtf_gr) %>%
    dplyr::filter(type == "CDS", transcript_id %in% transcript_ids) %>%
    dplyr::group_by(tx_id = transcript_id) %>%
    dplyr::summarise(cds_min = as.integer(min(start)), cds_max = as.integer(max(end)))
}

#' Genomic boundary of the region a terminal truncation keeps
#'
#' Walks `trunc_width` bases inwards from the truncated end of one transcript,
#' skipping over the introns, and reports the genomic coordinate where the
#' retained region stops. Everything on the `side` of that coordinate is what
#' txendcutr keeps.
#'
#' @param starts,ends Exon coordinates of a single transcript.
#' @param strand `"+"` or `"-"`.
#' @param trunc_width Truncation width, in transcript bases.
#' @param tx_end `"3p"` or `"5p"`.
#'
#' @return A list with `side` (`"high"` = the retained region is at the higher
#'   coordinates, `"low"` = at the lower ones) and the inclusive `boundary`.
retainedBoundary <- function(starts, ends, strand, trunc_width, tx_end = "3p"){
  ord <- order(starts)
  starts <- starts[ord]
  ends <- ends[ord]
  widths <- ends - starts + 1

  ## The truncated end sits at the higher coordinates for a 3' truncation of a
  ## + strand transcript, and for a 5' truncation of a - strand one.
  at_high <- identical(tx_end == "3p", strand == "+")

  ## A transcript shorter than the truncation is kept whole
  if(sum(widths) <= trunc_width){
    return(list(side = if(at_high) "high" else "low",
                boundary = if(at_high) min(starts) else max(ends)))
  }

  if(at_high){
    ## Bases from the start of each exon to the end of the transcript
    to_end <- rev(cumsum(rev(widths)))
    i <- max(which(to_end >= trunc_width))
    covered <- if(i < length(widths)) sum(widths[(i + 1):length(widths)]) else 0
    list(side = "high", boundary = ends[i] - (trunc_width - covered) + 1)
  }else{
    from_start <- cumsum(widths)
    i <- min(which(from_start >= trunc_width))
    covered <- if(i > 1) sum(widths[1:(i - 1)]) else 0
    list(side = "low", boundary = starts[i] + (trunc_width - covered) - 1)
  }
}

#' Truncation boundary of every bin, taken from its representative
#'
#' txendcutr collapses a bin on the truncated form of one transcript - the one
#' whose id names the bin - so that transcript defines the boundary the whole
#' bin is tagged against.
#'
#' @param exons Exon table of every drawn transcript, from `extractExons()`.
#' @param bin_members Bin membership table, from `queryBinMembers()`.
#' @param trunc_width Truncation width, in transcript bases.
#' @param tx_end `"3p"` or `"5p"`.
#'
#' @return A `tibble::tibble` with `bin_id`, `side` and `boundary`.
binBoundaries <- function(exons, bin_members, trunc_width, tx_end){
  reps <- bin_members %>% dplyr::filter(is_rep)

  purrr::map_dfr(seq_len(nrow(reps)), function(i){
    rep_exons <- exons %>% dplyr::filter(tx_id == reps$transcript_id[[i]])
    if(nrow(rep_exons) == 0) return(NULL)

    boundary <- retainedBoundary(rep_exons$start, rep_exons$end, rep_exons$strand[[1]],
                                 trunc_width, tx_end)

    tibble::tibble(bin_id = reps$bin_id[[i]], side = boundary$side, boundary = boundary$boundary)
  })
}

#' Every coordinate at which the exons of a transcript have to be cut
#'
#' Both tracks are drawn from one fragment table, so the exons are cut once, at
#' the union of the boundaries either track needs: where the retained region
#' starts and where the CDS starts and ends.
#'
#' @param bin_members Bin membership table, from `queryBinMembers()`.
#' @param boundaries Output of `binBoundaries()`.
#' @param cds_bounds Output of `extractCdsBounds()`.
#'
#' @return A `tibble::tibble` with `tx_id` and `cut`, the first base of a block.
collectCuts <- function(bin_members, boundaries, cds_bounds){
  truncation_cuts <- bin_members %>%
    dplyr::select(tx_id = transcript_id, bin_id) %>%
    dplyr::left_join(boundaries, by = "bin_id") %>%
    dplyr::filter(!is.na(boundary)) %>%
    ## For "high" the retained block starts on the boundary; for "low" it ends
    ## there, so it is the discarded block that starts one base further on
    dplyr::transmute(tx_id, cut = as.integer(ifelse(side == "high", boundary, boundary + 1)))

  dplyr::bind_rows(
    truncation_cuts,
    cds_bounds %>% dplyr::transmute(tx_id, cut = cds_min),
    cds_bounds %>% dplyr::transmute(tx_id, cut = cds_max + 1L)
  ) %>%
    dplyr::distinct()
}

#' Split exons at a set of per-transcript coordinates
#'
#' @param exons Exon table with `tx_id`, `start` and `end`.
#' @param cuts Output of `collectCuts()`; a cut is the first base of a block, and
#'   the ones falling outside their exon are ignored.
#'
#' @return The same table with the straddling exons replaced by their fragments.
splitExonsAtCuts <- function(exons, cuts){
  ex <- exons %>% dplyr::mutate(exon_id = dplyr::row_number())

  inner_cuts <- ex %>%
    dplyr::select(exon_id, tx_id, start, end) %>%
    dplyr::inner_join(cuts, by = "tx_id", relationship = "many-to-many") %>%
    dplyr::filter(cut > start, cut <= end) %>%
    dplyr::distinct(exon_id, frag_start = cut)

  ex %>%
    dplyr::select(exon_id, frag_start = start) %>%
    dplyr::bind_rows(inner_cuts) %>%
    dplyr::left_join(ex, by = "exon_id") %>%
    dplyr::arrange(exon_id, frag_start) %>%
    dplyr::group_by(exon_id) %>%
    ## Every fragment runs up to the base before the next cut, the last one up to
    ## the end of the exon it came from
    dplyr::mutate(frag_end = dplyr::coalesce(dplyr::lead(frag_start) - 1L, end)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(start = frag_start, end = frag_end) %>%
    dplyr::select(-exon_id, -frag_start, -frag_end)
}

#' Label every exon fragment for both tracks
#'
#' @param fragments Output of `splitExonsAtCuts()`, carrying `bin_id`, `tx_id`,
#'   `start`, `end` and `strand`.
#' @param boundaries Output of `binBoundaries()`.
#' @param cds_bounds Output of `extractCdsBounds()`.
#'
#' @return The same table with a `tag` column (`"Truncated"` / `"Other"`, the
#'   bottom track) and a `feature` column (`"5' UTR"` / `"CDS"` / `"3' UTR"` /
#'   `"Other"`, the top track).
annotateFragments <- function(fragments, boundaries, cds_bounds){
  fragments %>%
    dplyr::left_join(boundaries, by = "bin_id") %>%
    dplyr::left_join(cds_bounds, by = "tx_id") %>%
    dplyr::mutate(
      ## Bins without a usable representative keep their exons untagged
      tag = dplyr::case_when(
        is.na(boundary) ~ "Other",
        side == "high" & start >= boundary ~ "Truncated",
        side == "low" & end <= boundary ~ "Truncated",
        .default = "Other"
      ),
      feature = dplyr::case_when(
        is.na(cds_min) ~ "Other",
        start >= cds_min & end <= cds_max ~ "CDS",
        end < cds_min ~ ifelse(strand == "+", "5' UTR", "3' UTR"),
        start > cds_max ~ ifelse(strand == "+", "3' UTR", "5' UTR"),
        .default = "Other"
      )
    ) %>%
    dplyr::select(-side, -boundary, -cds_min, -cds_max)
}

#' Glue back the neighbouring fragments a track draws in the same colour
#'
#' Both tracks are cut at the union of their boundaries, so without this every
#' cut made for the OTHER track would show up as an internal border.
#'
#' @param ranges Rescaled exon fragments.
#' @param label_col Name of the column that colours the track.
#'
#' @return The same table with the adjacent same-coloured fragments merged.
mergeAdjacentRanges <- function(ranges, label_col){
  ranges %>%
    dplyr::arrange(tx_id, start) %>%
    dplyr::group_by(tx_id) %>%
    dplyr::mutate(block = cumsum(dplyr::coalesce(
      .data[[label_col]] != dplyr::lag(.data[[label_col]]) | start != dplyr::lag(end) + 1, TRUE
    ))) %>%
    dplyr::group_by(tx_id, block) %>%
    dplyr::summarise(dplyr::across(-c(start, end), dplyr::first),
                     start = min(start), end = max(end), .groups = "drop") %>%
    dplyr::select(-block)
}

#' Rescale the structures so the introns stop dominating the picture
#'
#' @param fragments Labelled exon fragments, from `annotateFragments()`.
#' @param exons_raw The same exons BEFORE they were cut.
#' @param gap_width Target intron width after rescaling.
#'
#' @return A `tibble::tibble` of exon fragments and introns in rescaled
#'   coordinates.
rescaleStructure <- function(fragments, exons_raw, gap_width = 100L){
  coord_cols <- c("tx_id", "seqnames", "start", "end", "strand", "type")

  ## The introns have to come from the UNCUT exons. Deriving them from the
  ## fragments would put a zero-width "intron" between the two halves of every
  ## exon that was cut, and would drop the intron a boundary falls into, leaving
  ## the block behind it floating away from the rest of the transcript.
  introns <- ggtranscript::to_intron(exons_raw[coord_cols], "tx_id")

  ## Grouping the rescaling by transcript alone keeps its exons and introns in
  ## one cumulative sum, so consecutive fragments stay flush against each other.
  ## `tag` / `feature` ride along as plain columns (NA on the introns, which are
  ## not filled). Both tracks are rescaled here at once, so they cannot drift.
  ggtranscript::shorten_gaps(fragments[c(coord_cols, "tag", "feature")], introns, "tx_id",
                             target_gap_width = as.integer(gap_width))
}

#' Draw a set of ENDome bins
#'
#' The one function this script is about. Everything else above is the machinery
#' it calls.
#'
#' @param bin_ids Bins to draw. They are numbered `Bin #1`, `Bin #2`, ... in this
#'   order, and do not have to belong to the same gene.
#' @param con Open DuckDB connection.
#' @param gtf_gr `GRanges` of the annotation, from section 1.2.
#' @param drop_superseded Hide the transcripts absorbed into a bin through an
#'   overlap.
#' @param show_truncation Draw the truncation track under every transcript.
#' @param gap_width Target intron width after rescaling.
#' @param trunc_width,tx_end Truncation settings; read off `run_id` when `NULL`.
#' @param title,subtitle Plot annotations. The default subtitle counts what was
#'   drawn.
#'
#' @return A `ggplot` object.
plotBins <- function(bin_ids, con, gtf_gr,
                     drop_superseded = FALSE, show_truncation = TRUE, gap_width = 100L,
                     trunc_width = NULL, tx_end = NULL,
                     title = NULL, subtitle = NULL){
  bin_ids <- unique(as.character(bin_ids))
  stopifnot(length(bin_ids) >= 1)

  run <- runSettings(con)
  if(is.null(trunc_width)) trunc_width <- run$trunc_width
  if(is.null(tx_end)) tx_end <- run$tx_end

  ## ---- Members of the bins -------------------------------------------------
  bin_members <- queryBinMembers(con, bin_ids)

  missing_bins <- setdiff(bin_ids, bin_members$bin_id)
  if(length(missing_bins) > 0){
    warning(sprintf("%i bin(s) not in the database: %s", length(missing_bins),
                    paste(head(missing_bins, 5), collapse = ", ")), call. = FALSE)
    bin_ids <- setdiff(bin_ids, missing_bins)
  }
  if(nrow(bin_members) == 0) stop("None of the bins is in the database.", call. = FALSE)

  if(drop_superseded){
    n_before <- nrow(bin_members)
    bin_members <- bin_members %>% dplyr::filter(!is_superseded)
    message(sprintf("Dropped %i superseded transcript(s).", n_before - nrow(bin_members)))
  }

  ## ---- Exon and CDS structure ----------------------------------------------
  exons <- extractExons(gtf_gr, bin_members$transcript_id)
  cds_bounds <- extractCdsBounds(gtf_gr, bin_members$transcript_id)

  missing_tx <- setdiff(bin_members$transcript_id, exons$tx_id)
  if(length(missing_tx) > 0){
    warning(sprintf("%i transcript(s) absent from the GTF and not drawn: %s", length(missing_tx),
                    paste(head(missing_tx, 5), collapse = ", ")), call. = FALSE)
    bin_members <- bin_members %>% dplyr::filter(!transcript_id %in% missing_tx)
  }
  if(nrow(exons) == 0) stop("None of the bin members is in the GTF.", call. = FALSE)

  exons <- exons %>%
    dplyr::left_join(bin_members %>% dplyr::select(tx_id = transcript_id, bin_id), by = "tx_id")

  ## ---- Coding structure and truncated regions ------------------------------
  boundaries <- if(show_truncation){
    binBoundaries(exons, bin_members, trunc_width, tx_end)
  }else{
    tibble::tibble(bin_id = character(), side = character(), boundary = integer())
  }

  fragments <- exons %>%
    splitExonsAtCuts(collectCuts(bin_members, boundaries, cds_bounds)) %>%
    annotateFragments(boundaries, cds_bounds)

  ## ---- Bin and transcript labels -------------------------------------------
  bin_labels <- bin_members %>%
    dplyr::distinct(bin_id, gene_name, n_members, info_score_norm) %>%
    dplyr::mutate(bin_index = match(bin_id, bin_ids)) %>%
    dplyr::arrange(bin_index) %>%
    dplyr::mutate(bin_label = sprintf("Bin #%i - %s\n(n = %i, I = %s)", bin_index,
                                      ifelse(is.na(gene_name), "?", gene_name), n_members,
                                      ifelse(is.na(info_score_norm), "NA",
                                             formatC(info_score_norm, format = "f", digits = 2))))

  tx_meta <- bin_members %>%
    dplyr::left_join(bin_labels %>% dplyr::select(bin_id, bin_index, bin_label), by = "bin_id") %>%
    dplyr::mutate(tx_label = ifelse(is_rep, paste0(transcript_id, " *"), transcript_id)) %>%
    dplyr::arrange(bin_index, dplyr::desc(is_rep), transcript_id)

  ## ---- Rescaled structure --------------------------------------------------
  ## `shorten_gaps()` rescales against the gaps between ALL the exons it is
  ## given, so it only handles one locus at a time - and a shared x axis only
  ## says anything when the bins sit on one anyway. Bins from different loci are
  ## therefore rescaled one by one and get an axis each.
  shared_x <- dplyr::n_distinct(paste(exons$seqnames, exons$strand)) == 1

  rescaled <- if(shared_x){
    rescaleStructure(fragments, exons, gap_width)
  }else{
    purrr::map_dfr(unique(fragments$bin_id), function(bin){
      rescaleStructure(dplyr::filter(fragments, bin_id == bin),
                       dplyr::filter(exons, bin_id == bin), gap_width)
    })
  }

  ## `to_intron()` / `shorten_gaps()` only carry the grouping columns through, so
  ## the bin annotation is joined back on afterwards.
  structure_rs <- rescaled %>%
    dplyr::inner_join(tx_meta %>% dplyr::select(tx_id = transcript_id, bin_label, tx_label),
                      by = "tx_id") %>%
    dplyr::mutate(
      bin_label = factor(bin_label, levels = bin_labels$bin_label),
      ## The first level of a discrete y axis is drawn at the bottom
      tx_label = factor(tx_label, levels = rev(tx_meta$tx_label))
    )

  ## ---- Plot ----------------------------------------------------------------
  if(is.null(subtitle)){
    subtitle <- sprintf("%i bins, %i transcripts", nrow(bin_labels), nrow(tx_meta))
  }
  caption <- paste0(
    run$run_id,
    if(show_truncation) sprintf(" | truncation: %i bp from the %s end", trunc_width, tx_end) else "",
    if(show_truncation) "\nTop half: coding structure. Bottom half: the region the truncation keeps." else "",
    "\n* bin representative (the transcript the bin is named after)",
    ## `facet_grid()` shares one x axis down a column, so say so rather than
    ## letting the reader read across bins that were rescaled separately
    if(!shared_x) "\nBins sit on different loci: each is rescaled on its own, so the x axis only reads within a bin." else ""
  )
  caption = ""

  drawBinStructurePlot(structure_rs, show_truncation, title, subtitle, caption)
}

#' Draw the two tracks of every transcript
#'
#' Each transcript gets two half-height tracks: its coding structure on top and
#' the region the truncation keeps underneath, with the introns running between
#' them.
#'
#' @param structure_rs Rescaled structure, from `plotBins()`.
#' @param show_truncation Draw the truncation track.
#' @param title,subtitle,caption Plot annotations.
#' @param colours,labels,heights Named vectors giving the fill, the legend label
#'   and the height of every category drawn by either track.
#'
#' @return A `ggplot` object.
drawBinStructurePlot <- function(structure_rs, show_truncation = TRUE,
                                 title = NULL, subtitle = NULL, caption = NULL,
                                 colours = feature_colours, labels = feature_labels,
                                 heights = feature_heights){
  exon_data <- dplyr::filter(structure_rs, type == "exon")
  intron_data <- dplyr::filter(structure_rs, type == "intron")

  ### Without the truncation track the coding structure takes the whole row
  height_scale <- if(show_truncation) 1 else 2

  ## `height` is an aesthetic of the geom, not a parameter, so it can be mapped
  ## from the data and vary from one feature to the next. It is resolved after
  ## the merge below, which is what guarantees a block has only one of them.
  scaleHeights <- function(ranges, label_col){
    ranges %>% dplyr::mutate(height = height_scale * dplyr::coalesce(
      unname(heights[as.character(.data[[label_col]])]), unname(heights[["Other"]])
    ))
  }

  feature_data <- mergeAdjacentRanges(exon_data, "feature") %>% scaleHeights("feature")
  ## Only the retained region is drawn underneath: the discarded part of the
  ## transcript is already there, in full, on the coding track above it
  trunc_data <- mergeAdjacentRanges(exon_data, "tag") %>%
    dplyr::filter(tag == "Truncated") %>%
    scaleHeights("tag")

  ### Only the categories that ended up on the plot reach the legend
  drawn <- c(as.character(feature_data$feature),
             if(show_truncation) as.character(trunc_data$tag))
  legend_breaks <- names(labels)[names(labels) %in% drawn]

  ## The scale draws one key per break, so two categories given the same label -
  ## `"5' UTR"` and `"3' UTR"` both labelled `"UTR"`, say - would otherwise get a
  ## key each. They collapse into one whenever they also share a fill; with
  ## different fills they have to stay apart, or the single key would lie about
  ## one of them. `col2rgb()` is what makes "#bebebe", "#BEBEBE" and "grey"
  ## count as the same colour.
  if(length(legend_breaks) > 0){
    key <- paste(labels[legend_breaks],
                 apply(grDevices::col2rgb(colours[legend_breaks]), 2, paste, collapse = "-"))
    unmergeable <- duplicated(labels[legend_breaks]) & !duplicated(key)
    if(any(unmergeable)){
      message("Kept separate legend keys for the same label with different fills: ",
              paste(unique(labels[legend_breaks][unmergeable]), collapse = ", "))
    }
    legend_breaks <- legend_breaks[!duplicated(key)]
  }

  bin_plot <- structure_rs %>%
    ggplot(aes(xstart = start, xend = end, y = tx_label)) +
    ## Every layer below draws a subset of the transcripts, and a free discrete
    ## y scale re-sorts its range ALPHABETICALLY as soon as a layer introduces a
    ## level the previous ones did not carry. Training the scale on the complete
    ## data first pins the panel order to the factor levels.
    geom_blank() +
    ## Top half: 5' UTR / CDS / 3' UTR
    ggtranscript::geom_half_range(data = feature_data, aes(fill = feature, height = height),
                                  range.orientation = "top")

  if(show_truncation){
    ## Bottom half: the region the truncation keeps
    bin_plot <- bin_plot +
      ggtranscript::geom_half_range(data = trunc_data, aes(fill = tag, height = height),
                                    range.orientation = "bottom")
  }

  bin_plot +
    ggtranscript::geom_intron(data = intron_data, aes(strand = strand),
                              arrow = grid::arrow(ends = "last", length = grid::unit(0.07, "inches"))) +
    scale_fill_manual(values = colours, breaks = legend_breaks, labels = labels[legend_breaks]) +
    facet_grid(rows = vars(bin_label), scales = "free_y", space = "free_y") +
    labs(x = "Rescaled Genomic Coordinates", y = "", fill = "",
         title = title, subtitle = subtitle, caption = caption) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0),
      panel.spacing.y = unit(0.8, "lines"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.caption = element_text(hjust = 0, colour = "grey30"),
      legend.position = "bottom"
    )
}

#' Figure height that fits the bins and transcripts a plot ended up with
#'
#' @param bin_plot Output of `plotBins()`.
#'
#' @return A height in inches.
binFigureHeight <- function(bin_plot){
  n_bins <- dplyr::n_distinct(bin_plot$data$bin_label)
  n_tx <- dplyr::n_distinct(bin_plot$data$tx_label)

  max(3, 1.1 * n_bins + 0.3 * n_tx + 1.2)
}

############################################################################## #
# ---- 1. Load data ----

#----------------------------------------------------------------------------- #
## 1.1 Open the ENDome DuckDB ----
if(!exists("con")) con <- DBI::dbConnect(duckdb::duckdb(), dbdir = duckdb_path, read_only = TRUE)

### Lazy tables to query in section 2. They behave like data frames until
### `collect()`, and the filtering happens inside DuckDB.
db_bins <- dplyr::tbl(con, "bins")
db_transcripts <- dplyr::tbl(con, "transcripts")
db_bin_information <- dplyr::tbl(con, "bin_information")

message(runSettings(con)$run_id)

#----------------------------------------------------------------------------- #
## 1.2 Load the annotation ----
### The exons draw the structure, the CDS separates it into UTRs and coding
if(!exists("gtf_gr")) gtf_gr <- rtracklayer::import(gtf_path, feature.type = c("exon", "CDS"))

############################################################################## #
# ---- 2. Select the bins to draw ----

### Bins of at least four members, most consistent first
bin_selection <- db_bin_information %>%
  dplyr::filter(n_members >= 3, n_members <= 3) %>%
  dplyr::arrange(-info_score_norm, -n_members) %>% 
  head(6) %>%
  dplyr::collect()

bin_selection

# Useful bins!
bin_selection <- db_bin_information %>%
dplyr::filter(n_members >= 3, n_members <= 3) %>%
  dplyr::arrange(-info_score_norm, -n_members) %>%
  head(6) %>%
  dplyr::collect() # Keep 3

bin_selection <- db_bin_information %>%
  dplyr::filter(n_members >= 4, n_members <= 5) %>%
  dplyr::arrange(info_score_norm, -n_members) %>%
  head(6) %>%
  dplyr::collect() # Keep 2

### Other ways of getting to a `bin_id`, all of them ending in a character vector
### handed to `plotBins()`:
###
### ## The bins whose members disagree the most
### db_bin_information %>% filter(n_members >= 4) %>% arrange(info_score_norm) %>% head(6) %>% collect()
###
### ## Every bin of one gene, largest first
### db_transcripts %>% filter(gene_name == "CYRIB") %>% distinct(bin_id) %>%
###   inner_join(db_bin_information, by = "bin_id") %>% arrange(desc(n_members)) %>% collect()
###
### ## One bin, by hand
### tibble::tibble(bin_id = "ENST00000519110.5")

############################################################################## #
# ---- 3. Visualizations ----

bin_plot <- plotBins(bin_selection$bin_id[[2]], con, gtf_gr, title = "", subtitle = "")
bin_plot

saveBothExt(bins_fig_path, bin_plot, figure_label,
            width = figure_width, height = binFigureHeight(bin_plot))

############################################################################## #
# ---- 4. Export results ----

### The bins drawn in the figure, in the order they were numbered
bin_selection %>%
  dplyr::select(run_id, bin_id, n_members, info_score_norm, info_score) %>%
  dplyr::mutate(bin_index = dplyr::row_number(), .before = bin_id) %>%
  readr::write_tsv(bin_summary_output_path(figure_label))

############################################################################## #
# ---- 5. Clean up ----

### Only once the session is done with the database - every section above reuses
### the open connection.
DBI::dbDisconnect(con, shutdown = TRUE); rm(con)
