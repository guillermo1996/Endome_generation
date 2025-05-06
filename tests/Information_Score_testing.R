## _________________________________________________
##
## Script title
##
## Aim:
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-04-10
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-04-10)
##
## _________________________________________________
##
## - Notes:
##
## - Changelog:
##
## - Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

############################################################################## #
# ---- 0. Setup ----

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(rtracklayer))
shhh(library(GenomicFeatures))
shhh(library(GenomicRanges))
shhh(library(plyranges))
shhh(library(tidyverse))

#----------------------------------------------------------------------------- #
## 0.2 Script Paths ----
main_path <- here::here()

### Input Paths
gtf_path <- file.path(main_path, "results/ORF_Filter/sq3.annotated_orf.filter.gtf")
gtf_3p_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.gtf")
overlap_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.overlaps.tsv")
merge_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.merge.tsv")
sq3_class_path<- file.path(main_path, "results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt")

### Output Paths

#### Create output directory
# dir.create(results_path, showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.3 Script Parameters ----

#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----

### G: All parameters should be changeable in Section '0.3 Script Parameters',
### but the user may want to pass arguments when calling the script using
### 'Rscript script.R --argument value'. To do so, use the following 'optparse'
### template.
#option_list = list(
#  make_option(c("-t", "--test"), default = NULL, type = "double", help = "test")
#)
#
#opt_parser = OptionParser(option_list = option_list)
#opt = parse_args(opt_parser)
#
#if(!is.null(opt)) test = opt

#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----

############################################################################## #
# ---- 1. Load all the needed files ----
gtf <- rtracklayer::import(gtf_path)
gtf_3p <- rtracklayer::import(gtf_3p_path)
overlap_table <- vroom::vroom(overlap_table_path)
merge_table <- vroom::vroom(merge_table_path) %>% dplyr::select(tx_in, tx_out) %>% dplyr::filter(tx_in != tx_out)
sq3_class <- vroom::vroom(sq3_class_path)


############################################################################## #
# ---- 2. Map Transcripts - Truncated - Binned ----

#----------------------------------------------------------------------------- #
## Generate the initial tables ----

tx_tbl <- gtf %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(seqnames, start, end, strand, type, transcript_id, gene_id, ORF_cat, hasCDS) %>%
  dplyr::mutate(type = "transcript")

tr_tbl <- gtf_3p %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(truncation_id = paste0("txcutr_", row_number())) %>%
  dplyr::select(truncation_id, transcript_id, gene_id, seqnames, start, end, strand, type) %>%
  dplyr::mutate(type = "truncation")

#----------------------------------------------------------------------------- #
## Map transcripts to truncated ----
mapping_TxTr <- function(tr_tbl, overlap_table = "", n_threads = 32, output_path = "", overwrite = F){
  if(file_test("-f", output_path) & !overwrite){
    map_tx_tr <- readRDS(output_path)
    return(map_tx_tr)
  }

  # We aim to genereate a mapping table between transcripts IDs and truncations
  # IDs. Two steps:
  #
  # Step 1. From the truncated table we extract the mapping between trancripts
  # and truncations.
  #
  # Step 2. There are transcripts that, when truncated, are merged into other
  # truncations. To extract the mapping between these transcripts and the final
  # truncation, we use the "overlap_table" reference generated from gtxcutr.

  ## Step 1: From the truncation table:
  map_tx_tr_s1 <- tr_tbl %>% dplyr::select(transcript_id, truncation_id)

  ## Step 2: From the "overlap_table":
  ##
  ## The approach proposed is the following:
  ##
  ## 1. Convert the "overlap_table" into a directed graph. Each transcript is a
  ## node and their overlaps are edges.
  ##
  ## 2. Identify "terminal_nodes" (those that are found in the truncation table)
  ## and "input_nodes" (those that get merged into a terminal node).
  ##
  ## 3. Calculate the distance between all input nodes and temrinal nodes. Each
  ## input node will match to a single terminal node with a distance of 1.
  ##
  ## 4. Split the input node in batches to parallelize the terminal node
  ## extraction. The idea is the identify which terminal node correspond to each
  ## input node.
  ##
  ## 5. Merge the final results into a single dataframe that complements the
  ## results from Step 1.

  ### 1. Create a directed graph
  library(igraph)
  g <- graph_from_data_frame(overlap_table, directed = T)

  ### 2. Identify terminal ndoes and input nodes
  all_nodes <- V(g)$name
  is_terminal <- degree(g, mode = "out") == 0
  terminal_nodes <- all_nodes[is_terminal]
  input_nodes <- all_nodes[!is_terminal]

  ### 3. Calculate the distance
  distances <- distances(g, v = input_nodes, to = terminal_nodes, mode = "out")

  ### 4. Split the input nodes into batches and run the input to output mapping
  input_nodes_batches <- split(input_nodes, ceiling(seq_along(input_nodes)/(length(input_nodes)/n_threads)))
  bpparam = BiocParallel::MulticoreParam(workers = n_threads, progressbar = T)

  map_tx_tr_s2 <- BiocParallel::bplapply(input_nodes_batches, BPPARAM = bpparam, function(in_batch){
    batch_results <- lapply(in_batch, function(input_node){
      i <- which(input_nodes == input_node)
      output_node <- terminal_nodes[which(distances[i, ] == 1)]

      return(output_node)
    }) %>% setNames(., in_batch)
    batch_df <- tibble::enframe(purrr::map_chr(batch_results, identity), name = "input_nodes", value = "output_nodes")

    return(batch_df)
  }) %>% dplyr::bind_rows()

  ### 5. Merge the final results
  map_tx_tr <- map_tx_tr_s2 %>%
    dplyr::left_join(tr_tbl %>% dplyr::select(transcript_id, truncation_id), by = c("output_nodes" = "transcript_id")) %>%
    dplyr::select(transcript_id = input_nodes, truncation_id) %>%
    dplyr::bind_rows(map_tx_tr_s1) %>%
    dplyr::mutate(truncation_id = factor(truncation_id, levels = unique(tr_tbl$truncation_id))) %>%
    dplyr::arrange(truncation_id) %>%
    dplyr::mutate(truncation_id = as.character(truncation_id))

  if(output_path != "") map_tx_tr %>% saveRDS(output_path)
  return(map_tx_tr)
}

map_tx_tr_path <- file.path(main_path, "tests/map_tx_tr.rds")
map_tx_tr <- mapping_TxTr(tr_tbl = tr_tbl, overlap_table = overlap_table, n_threads = 32, output_path = map_tx_tr_path)

# Update transcript table with the reference to truncation table
tx_tbl <- tx_tbl %>%
  dplyr::left_join(map_tx_tr, by = "transcript_id") %>%
  dplyr::relocate(transcript_id, truncation_id, gene_id)

#----------------------------------------------------------------------------- #
## Map binned to transcripts and truncated ----
mapping_TrBin <- function(bin_tbl, overlap_table = "", n_threads = 32, output_path = "", overwrite = F){
  if(file_test("-f", output_path) & !overwrite){
    map_tr_bin <- readRDS(output_path)
    return(map_tr_bin)
  }

  ## Step 1: From the truncation table:
  map_tr_bin_s1 <- bin_tbl %>% dplyr::select(truncation_id, bin_id)

  merge_table2 <- merge_table %>%
    dplyr::left_join(map_tx_tr, by = c("tx_in" = "transcript_id")) %>%
    dplyr::left_join(map_tx_tr, by = c("tx_out" = "transcript_id")) %>%
    dplyr::select(tr_in = truncation_id.x, tr_out = truncation_id.y)


  library(igraph)
  g <- graph_from_data_frame(merge_table2, directed = T)

  ### 2. Identify terminal ndoes and input nodes
  all_nodes <- V(g)$name
  is_terminal <- degree(g, mode = "out") == 0
  terminal_nodes <- all_nodes[is_terminal]
  input_nodes <- all_nodes[!is_terminal]

  ### 3. Calculate the distance
  distances <- distances(g, v = input_nodes, to = terminal_nodes, mode = "out")

  ### 4. Split the input nodes into batches and run the input to output mapping
  input_nodes_batches <- split(input_nodes, ceiling(seq_along(input_nodes)/(length(input_nodes)/n_threads)))
  bpparam = BiocParallel::MulticoreParam(workers = n_threads, progressbar = T)

  map_tr_bin_s2 <- BiocParallel::bplapply(input_nodes_batches, BPPARAM = bpparam, function(in_batch){
    batch_results <- lapply(in_batch, function(input_node){
      i <- which(input_nodes == input_node)
      output_node <- terminal_nodes[which(distances[i, ] == 1)]

      return(output_node)
    }) %>% setNames(., in_batch)
    batch_df <- tibble::enframe(purrr::map_chr(batch_results, identity), name = "input_nodes", value = "output_nodes")

    return(batch_df)
  }) %>% dplyr::bind_rows()

  ### 5. Merge the final results
  map_tx_tr <- map_tx_tr_s2 %>%
    dplyr::left_join(tr_tbl %>% dplyr::select(transcript_id, truncation_id), by = c("output_nodes" = "transcript_id")) %>%
    dplyr::select(transcript_id = input_nodes, truncation_id) %>%
    dplyr::bind_rows(map_tx_tr_s1) %>%
    dplyr::mutate(truncation_id = factor(truncation_id, levels = unique(tr_tbl$truncation_id))) %>%
    dplyr::arrange(truncation_id) %>%
    dplyr::mutate(truncation_id = as.character(truncation_id))

  if(output_path != "") map_tx_tr %>% saveRDS(output_path)
  return(map_tx_tr)
}

bin_tbl <- tr_tbl %>%
  dplyr::filter(transcript_id %in% unique(merge_table$tx_out)) %>%
  dplyr::mutate(bin_id = paste0("bin_", row_number())) %>%
  dplyr::relocate(bin_id) %>%
  dplyr::mutate(type = "bin")


map_tr_bin_s1 <- bin_tbl %>% dplyr::select(truncation_id, bin_id)

merge_table2 <- merge_table %>%
  dplyr::mutate(tx_in = recode(tx_in, !!!lv))

g <- graph_from_data_frame(merge_table)
all_nodes <- V(g)$name
is_terminal <- degree(g, mode = "out") == 0
terminal_nodes <- all_nodes[is_terminal]
input_nodes <- all_nodes[!is_terminal]

distances <- distances(g, v = input_nodes, to = terminal_nodes, mode = "out")
n_threads <- 32
input_nodes_batches <- split(input_nodes, ceiling(seq_along(input_nodes)*n_threads/length(input_nodes)))

map_tx_tr_s2 <- BiocParallel::bplapply(input_nodes_batches, function(in_batch){
  batch_results <- lapply(in_batch, function(input_node){
    i <- which(input_nodes == input_node)
    output_node <- terminal_nodes[which(is.finite(distances[i, ]))]

    return(output_node)
  }) %>% setNames(., in_batch)
  batch_df <- tibble::enframe(purrr::map_chr(batch_results, identity), name = "input_nodes", value = "output_nodes")

  return(batch_df)
}, BPPARAM = BiocParallel::MulticoreParam(workers = n_threads, progressbar = T)) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(map_tr_bin_s1, by = c("output_nodes" = ""))

map_tx_tr
map_tx_tr_s2 %>% dplyr::filter(input_nodes != output_nodes)









map_tr_bin_1 <- bin_tbl %>% dplyr::select(truncation_id, bin_id)
map_tx_bin_1 <- bin_tbl %>% dplyr::select(transcript_id, bin_id)

map_tr_bin <- merge_table %>% dplyr::left_join(map_tx_bin_1, by = c("tx_out" = "transcript_id")) %>%
  dplyr::left_join(map_tx_tr, by = c("tx_in" = "transcript_id")) %>%
  dplyr::select(truncation_id, bin_id)

merge_table %>%
  dplyr::left_join(bin_tbl %>% dplyr::select())







# input_transcripts <- gtf %>% plyranges::filter(type == "transcript")
# output_transcripts <- truncted_transcripts %>% plyranges::filter(!transcript_id %in% unique(merge_table$tx_out))


# Table between Raw Transcript and ORF
orf_sequence_tbl <- sq3_class %>%
  dplyr::filter(!is.na(ORF_seq)) %>%
  dplyr::distinct(ORF_seq) %>%
  dplyr::mutate(ORF_id = paste0("orf_", row_number()), .before = ORF_seq)

map_transcript_orf <- sq3_class %>%
  dplyr::filter(!is.na(ORF_seq)) %>%
  dplyr::distinct(isoform, ORF_seq) %>%
  dplyr::left_join(orf_sequence_tbl, by = "ORF_seq") %>%
  dplyr::select(-ORF_seq)

raw_transcript_tbl <- gtf %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(seqnames, start, end, strand, type, transcript_id, gene_id, ORF_cat, hasCDS, ORF_length) %>%
  dplyr::left_join(map_transcript_orf, by = c("transcript_id" = "isoform")) %>%
  dplyr::relocate(transcript_id, ORF_id, gene_id, seqnames, start, end, strand, type, ORF_cat, hasCDS, ORF_length)


# Table between Raw Transcripts and Truncated Transcripts
truncated_transcripts <- gtf_3p %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(tx_transcript_id = paste0("txcutr_", row_number())) %>%
  dplyr::select(tx_transcript_id, parent_transcript_id = transcript_id, gene_id, seqnames, start, end, strand, type)

# Reduce the Overlap Table
library(igraph)
g <- graph_from_data_frame(overlap_table, directed = T)
all_nodes <- V(g)$name
is_terminal <- degree(g, mode = "out") == 0
terminal_nodes <- all_nodes[is_terminal] # tn2 <- intersect(unique(overlap_table$subjectHits), truncated_transcripts$parent_transcript_id )

bpparam = BiocParallel::SerialParam(progressbar = T)
map_overlaps <- BiocParallel::bplapply(terminal_nodes, BPPARAM = bpparam, function(end_node){
  input_nodes <- igraph::subcomponent(g, end_node, mode = "in")
  input_nodes <- setdiff(input_nodes$name, end_node)

  return(
    tibble::tibble(input_nodes = input_nodes,
                   output_nodes = rep(end_node, length(input_nodes)))
  )
}) %>% dplyr::bind_rows()

map_transcript_truncated <- truncated_transcripts %>%
  dplyr::select(transcript_id = parent_transcript_id, tx_transcript_id)

map_transcript_overlaps <- map_overlaps %>%
  dplyr::left_join(map_transcript_truncated, by = c("output_nodes" = "transcript_id")) %>%
  dplyr::select(transcript_id = input_nodes, tx_transcript_id)

map_transcript_truncated <- dplyr::bind_rows(
  map_transcript_truncated,
  map_transcript_overlaps
)

raw_transcript_tbl <- raw_transcript_tbl %>%
  dplyr::left_join(map_transcript_truncated, by = "transcript_id") %>%
  dplyr::relocate(transcript_id, ORF_id, tx_transcript_id, gene_id)




