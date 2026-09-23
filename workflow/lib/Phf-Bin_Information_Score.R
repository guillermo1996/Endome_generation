## _________________________________________________
##
## Helper Functions: Bin Information Content
##
## Aim: Include into a single file the functions shared across the "Bin
## Information Content" scripts related to the Bin Information Content construction
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 05/08/2026
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.2 (19/08/2026)
## _________________________________________________
##
## Notes:
##     - VENDORED COPY of HelperFunctions/Phf-ENDome_Pipeline/Phf-Bin_Information_Score.R.
##       See the note in lib/Phf-DuckDB.R.
##     - The two plotting functions of the original are NOT vendored: they depend
##       on prettier_limits()/prettier_breaks() from hf_graph_and_themes.R, which
##       would drag a third helper file into the pipeline env. S04-ENDome_Report.Rmd
##       plots inline instead.
##
## Changelog:
##     - v1.2 (19/08/2026): Vendored for the Snakemake port. `bin_members` is now a
##       materialised table, so createBinMembership() reads it instead of
##       reconstructing the mapping. Added resolveTerminalNodes() (moved from T01,
##       with a diagnostic instead of a bare stopifnot).
##     - v1.1 (18/08/2026): Removed DuckDB functions to its own helper file
##     - v1.0 (05/08/2026): Initial release
##
## Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

################################################################################
# ---- 1. Required libraries ----
suppressWarnings(suppressMessages(library(dplyr)))

################################################################################
# ---- 2. Transcript graph resolution ----

#' Resolve every node of a directed graph to the terminal node of its component
#'
#' Used on the txendcutr overlaps table, where an edge means "this transcript was
#' superseded by that one". Each weakly connected component is expected to funnel
#' into exactly one sink.
#'
#' @param g igraph object.
#' @param report_self boolean, whether to keep the rows where a node maps to
#'   itself (i.e. the terminal nodes). Defaults to TRUE.
#'
#' @returns Tibble with columns \code{node} and \code{terminal_node}.
#' @export
resolveTerminalNodes <- function(g, report_self = TRUE){
  # Calculate the node component and degree
  comp <- igraph::components(g, mode = "weak")
  out_deg <- igraph::degree(g, mode = "out")
  node_names <- igraph::V(g)$name

  # Extract the terminal nodes and ensure there is a single terminal node per component
  terminal_idx <- which(out_deg == 0)

  # The prototype asserted this with a bare stopifnot(). A component with two
  # sinks (A->B, A->C) is not obviously impossible on real data, and when it
  # happens the caller needs to know WHICH transcripts are involved.
  duplicated_comps <- unique(comp$membership[terminal_idx][duplicated(comp$membership[terminal_idx])])
  if (length(duplicated_comps) > 0) {
    offenders <- node_names[terminal_idx][comp$membership[terminal_idx] %in% duplicated_comps]
    stop(
      sprintf(
        paste("%i overlap component(s) resolve to more than one terminal transcript, so the",
              "superseded-by mapping is ambiguous.\n  Offending terminal nodes: %s"),
        length(duplicated_comps),
        paste(utils::head(offenders, 20), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Match each node to a terminal node and create the tibble
  match_pos <- match(comp$membership, comp$membership[terminal_idx])

  in_out_nodes <- tibble::tibble(
    node = node_names,
    terminal_node = node_names[terminal_idx][match_pos]
  )

  # Remove self-maps if expressed
  if(!report_self) in_out_nodes <- in_out_nodes %>% dplyr::filter(node != terminal_node)
  in_out_nodes
}

################################################################################
# ---- 3. Bin membership ----

#' Map every bin to the transcripts it contains
#'#'
#' @param db named list of lazy tibbles, as returned by \code{getTbls()}.
#' @param include_superseded boolean, whether to keep the superseded transcripts.
#'   Defaults to TRUE.
#'
#' @returns Lazy tibble with \code{bin_id}, \code{transcript_id},
#'   \code{is_superseded} and \code{superseded_by}.
#' @export
createBinMembership <- function(db, include_superseded = TRUE){
  bin_members <- db$transcripts %>% 
    dplyr::select(bin_id, transcript_id, superseded_by) %>% 
    dplyr::mutate(is_superseded = if_else(!is.na(superseded_by), TRUE, FALSE), .before = superseded_by)

  if(!include_superseded){
    bin_members <- bin_members %>% dplyr::filter(!is_superseded)
  }

  bin_members
}

#' Map every bin to the protein sequences of its transcripts
#'
#' Extends \code{createBinMembership()} by walking the
#' transcripts -> cds -> proteins relationship. Transcripts without an ORF (i.e.
#' with no \code{aa_id}) are dropped.
#'
#' @param db named list of lazy tibbles, as returned by \code{getTbls()}.
#' @param include_superseded boolean, see \code{createBinMembership()}.
#'
#' @returns Lazy tibble with columns \code{bin_id}, \code{transcript_id},
#'   \code{aa_id} and \code{aa_seq}.
#' @export
createBinProteins <- function(db, include_superseded = FALSE){
  createBinMembership(db, include_superseded = include_superseded) %>%
    dplyr::inner_join(dplyr::select(db$transcripts, transcript_id, cds_id), by = "transcript_id") %>%
    dplyr::inner_join(dplyr::select(db$cds, cds_id, aa_id), by = "cds_id") %>%
    dplyr::filter(!is.na(aa_id)) %>%
    dplyr::inner_join(dplyr::select(db$proteins, aa_id, aa_seq), by = "aa_id") %>%
    dplyr::select(bin_id, transcript_id, aa_id, aa_seq)
}

#' Export protein sequences to a FASTA file keyed by their content id
#'
#' @param proteins dataframe or lazy tibble with \code{aa_id} and \code{aa_seq}.
#' @param fasta_path string, path of the FASTA file to write.
#'
#' @returns Invisibly, the \code{AAStringSet} that was written to disk.
#' @export
writeBinProteinFasta <- function(proteins, fasta_path){
  protein_fasta <- proteins %>%
    dplyr::select(aa_id, aa_seq) %>%
    dplyr::filter(!is.na(aa_seq), nzchar(aa_seq)) %>%
    dplyr::distinct() %>%
    dplyr::collect() %>%
    dplyr::mutate(aa_seq = sub("\\*+$", "", aa_seq))

  seqs <- Biostrings::AAStringSet(stats::setNames(protein_fasta$aa_seq, protein_fasta$aa_id))

  dir.create(dirname(fasta_path), showWarnings = FALSE, recursive = TRUE)
  Biostrings::writeXStringSet(seqs, fasta_path)

  return(invisible(seqs))
}

################################################################################
# ---- 4. Bin PID similarity ----

#' Intra-bin protein similarity from pairwise alignments
#'
#' Computes the sequence similarity of a single bin by aligning every pair of
#' its protein sequences with \code{pwalign::pairwiseAlignment()} and averaging
#' the resulting percentage identity (PID) and alignment scores.
#'
#' Sequences outside the requested length range are discarded before the
#' alignment, since very short sequences inflate the PID and very long ones make
#' the alignment prohibitively expensive.
#'
#' @param df dataframe with the proteins of a single bin. Requires the
#'   \code{aa_id} and \code{aa_seq} columns, and optionally \code{bin_id}.
#' @param min_width numeric, minimum sequence length (aa) to be comparable.
#' @param max_width numeric, maximum sequence length (aa) to be comparable.
#' @param type string, alignment type passed to \code{pwalign::pairwiseAlignment()}.
#'
#' @returns Single-row tibble with the similarity metrics of the bin.
#' @export
measureSeqSimilarity <- function(df, min_width = 30, max_width = 10000, type = "global"){
  # Create the AAStringSet
  seqs <- Biostrings::AAStringSet(stats::setNames(df$aa_seq, df$aa_id))

  # Filter short and long sequences
  short_seqs <- BiocGenerics::width(seqs) < min_width
  long_seqs <- BiocGenerics::width(seqs) > max_width
  good_seqs <- seqs[!short_seqs & !long_seqs]

  # Return if bin similarity cannot be computed
  if(length(good_seqs) == 1){
    return(generateSimOutput(df, 1, 0, NA, short_seqs, long_seqs, type))
  }else if(length(good_seqs) == 0){
    return(generateSimOutput(df, NA, NA, NA, short_seqs, long_seqs, type))
  }

  # Create the combination of sequences
  seq_comb <- utils::combn(seq_along(good_seqs), 2)

  # Run the pairwise alignment and extract metrics
  pw_alignment <- pwalign::pairwiseAlignment(
    good_seqs[seq_comb[1, ]], good_seqs[seq_comb[2, ]],
    gapOpening = 11, gapExtension = 1, type = type, substitutionMatrix = "BLOSUM62")
  pw_pid <- pwalign::pid(pw_alignment)

  sim_pid <- mean(pw_pid)/100
  sd_pid <- stats::sd(pw_pid)/100
  sim_score <- mean(pwalign::score(pw_alignment))

  # Generate the output data.frame
  sim_output <- generateSimOutput(df, sim_pid, sd_pid, sim_score, short_seqs, long_seqs, type)
  return(sim_output)
}

#' Assemble the similarity report of a single bin
#'
#' @param df dataframe with the proteins of a single bin.
#' @param sim_pid numeric, mean percentage identity, as a proportion.
#' @param sd_pid numeric, standard deviation of the percentage identity.
#' @param sim_score numeric, mean alignment score.
#' @param short_seqs logical vector (or count) of sequences dropped as too short.
#' @param long_seqs logical vector (or count) of sequences dropped as too long.
#' @param type string, alignment type employed.
#'
#' @returns Single-row tibble.
#' @export
generateSimOutput <- function(df, sim_pid, sd_pid, sim_score, short_seqs = 0, long_seqs = 0, type = "global"){
  # Get the bin ID
  bin_id <- if (tibble::has_name(df, "bin_id")) unique(df$bin_id) else ""

  # Add comments about the results
  if(is.na(sim_pid)){
    comment <- "No sequences left after filtering by length."
  }else if(is.na(sim_score)){
    if(sum(short_seqs) != 0 || sum(long_seqs) != 0){
      comment <- "One sequence left after filtering by length."
    }else{
      comment <- "Single sequence in the bin."
    }
  }else{
    comment <- ""
  }

  # Generate the output data.frame
  tibble::tibble(
    bin_id = bin_id,
    n_orf = nrow(df),
    n_comp = n_orf - sum(short_seqs) - sum(long_seqs),
    mean_orf_width = mean(nchar(df$aa_seq)),
    sim_pid = sim_pid,
    sd_pid = sd_pid,
    sim_score = sim_score,
    short_seqs = sum(short_seqs),
    long_seqs = sum(long_seqs),
    type = type,
    comments = comment
  )
}

#' Intra-bin protein similarity across every bin
#'
#' Splits the bin proteins per bin and applies \code{measureSeqSimilarity()} to
#' each of them. Bins are grouped into clusters before being dispatched to
#' \code{BiocParallel}.
#'
#' @param bin_proteins dataframe or lazy tibble with \code{bin_id}, \code{aa_id}
#'   and \code{aa_seq}, as returned by \code{createBinProteins()}.
#' @param BPPARAM BiocParallelParam to dispatch the bins with.
#' @param clusters_per_worker numeric, clusters of bins to generate per worker.
#' @param ... further arguments passed to \code{measureSeqSimilarity()}.
#'
#' @returns Tibble with one row per bin. See \code{generateSimOutput()}.
#' @export
measureBinSimilarity <- function(bin_proteins, BPPARAM = BiocParallel::bpparam(), clusters_per_worker = 2, ...){
  # Split the proteins per bin
  bin_list <- bin_proteins %>%
    dplyr::select(bin_id, aa_id, aa_seq) %>%
    dplyr::collect() %>%
    S4Vectors::split(., .$bin_id)

  # Group the bins into clusters to reduce the dispatching overhead
  n_clusters <- max(1, BiocParallel::bpnworkers(BPPARAM) * clusters_per_worker)
  bin_clusters <- split(bin_list, ceiling(seq_along(bin_list)/(length(bin_list)/n_clusters)))

  # Measure the similarity of every bin of every cluster
  bin_sim <- BiocParallel::bplapply(bin_clusters, function(bin_cluster, ...){
    purrr::map(bin_cluster, measureSeqSimilarity, ...) %>% dplyr::bind_rows()
  }, ..., BPPARAM = BPPARAM) %>% dplyr::bind_rows()

  return(bin_sim)
}

################################################################################
# ---- 5. Bin Rao's entropy calculation ----

#' Distance between two ratio-scale variables, bounded in [0, 1]
#'
#' @param x,y numeric vectors of the feature
#' @param scaling string, "ratio" for 1 - min/max or "canberra" for
#'   |x-y|/(x+y). Defaults to "ratio".
#'
#' @returns Numeric vector in [0, 1], NA whenever either value is NA.
#' @export
relativeDelta <- function(x, y, scaling = "ratio"){
  switch(scaling,
    ratio = abs(x - y) / pmax(x, y, 1),
    canberra = abs(x - y) / pmax(x + y, 1),
    stop("Unknown `scaling`: ", scaling)
  )
}

#' Graded distance between the initiation/termination codons of two ORFs
#'
#' Three levels: the ORFs disagree completely (1) when only one of them
#' initiates/terminates canonically, partially (0.5) when both are in the same
#' class but use a different codon (e.g. TGA vs TAA), and not at all (0) when
#' they use the same codon.
#'
#' @param canonical_i,canonical_j logical vectors, TRUE when the codon is
#'   canonical (ATG to start; TAA/TAG/TGA to stop).
#' @param sequence_i,sequence_j character vectors with the codon itself.
#'
#' @returns Numeric vector in {0, 0.5, 1}, NA whenever either member has no ORF.
#' @export
codonDissimilarity <- function(canonical_i, canonical_j, sequence_i, sequence_j){
  dplyr::case_when(
    is.na(canonical_i) | is.na(canonical_j) ~ NA_real_,
    canonical_i != canonical_j ~ 1, # only one ORF initiates/terminates properly
    sequence_i != sequence_j ~ 0.5, # same class, different codon
    .default = 0 # same codon
  )
}

#' Discordance of a symmetric qualitative variable
#'
#' NA propagates, which is what \code{combineGower()} wants: a variable that
#' cannot be compared for a pair drops out of that pair's weighting.
#' @export
discordance <- function(x, y) as.numeric(x != y)

#' Gower's general coefficient of similarity, expressed as a distance
#'
#' d_ij = sum_k (w_k * delta_ijk * d_ijk) / sum_k (w_k * delta_ijk), with
#' delta_ijk = 0 when variable k cannot be compared for the pair (missing value,
#' or an asymmetric variable in its reference state in both members). The
#' weights of the non-comparable variables are therefore set to 0 and the rest
#' renormalised, so a pair scored on fewer variables is neither penalised nor
#' rewarded. Pairs with no comparable variable return NA.
#'
#' @param df dataframe holding one \code{d_<variable>} column per variable.
#' @param weights named numeric vector of variable weights. Names must match the
#'   \code{d_<variable>} columns of \code{df}; zero-weighted variables are ignored.
#'
#' @returns Numeric vector of length \code{nrow(df)} with the combined d_ij.
#' @export
combineGower <- function(df, weights){
  weights <- weights[weights > 0]
  missing_cols <- setdiff(paste0("d_", names(weights)), colnames(df))
  if (length(missing_cols) > 0) {
    stop("Weighted variable(s) with no distance column: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  components <- as.matrix(df[paste0("d_", names(weights))])
  stopifnot(all(components >= 0 & components <= 1, na.rm = TRUE))

  # Remove weights of non-comparable variables
  weight_mat <- matrix(weights, nrow = nrow(components), ncol = length(weights), byrow = TRUE)
  weight_mat[is.na(components)] <- 0
  total_weight <- rowSums(weight_mat)

  dplyr::if_else(
    total_weight > 0,
    rowSums(components*weight_mat, na.rm = TRUE) / total_weight,
    NA_real_
  )
}
