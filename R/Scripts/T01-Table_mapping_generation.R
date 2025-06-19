## _________________________________________________
##
## Table Mapping Script
##
## Aim: To generate the ORF, Transcript, Truncation and Bin tables that are
## employed in the study. This will not be the final implementation (hopefully).
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
shhh(library(Biostrings))
shhh(library(pwalign))
shhh(library(rtracklayer))
shhh(library(GenomicFeatures))
shhh(library(GenomicRanges))
shhh(library(plyranges))
shhh(library(tidyverse))

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
utr_site <- "3p"

#----------------------------------------------------------------------------- #
## 0.3 Script Paths ----
main_path <- "~/RytenLab-Research/38-Endome_generation"
script_path <- file.path(main_path, "R")

### Input Paths
gtf_path <- file.path(main_path, "results/ORF_Filter/sq3.annotated_orf.filter.gtf")
sq3_class_path <- file.path(main_path, "results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt")

if(utr_site == "3p"){
  gtf_trunc_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.gtf")
  overlap_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.overlaps.tsv")
  merge_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.merge.tsv")
}else if (utr_site == "5p"){
  gtf_trunc_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w300.5p.gtf")
  overlap_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w300.5p.overlaps.tsv")
  merge_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w300.5p.merge.tsv")
}else{
  stop("No valid truncation direction")
}

### Output Paths
# results_path <- file.path(script_path, "Results/Table_mapping_3p")
results_path <- file.path(script_path, paste0("Results/Table_mapping_", utr_site))

tx_tbl_path <- file.path(results_path, "tx_tbl.rds")
tr_tbl_path <- file.path(results_path, "tr_tbl.rds")
bin_tbl_path <- file.path(results_path, "bin_tbl.rds")
orf_tbl_path <- file.path(results_path, "orf_tbl.rds")

map_tx_tr_path <- file.path(results_path, "map_tx_tr.rds")
map_tr_bin_path <- file.path(results_path, "map_tr_bin.rds")

#### Create output directory
dir.create(results_path, showWarnings = F, recursive = T)

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
sq3_class <- vroom::vroom(sq3_class_path)

gtf_trunc <- rtracklayer::import(gtf_trunc_path)
overlap_table <- vroom::vroom(overlap_table_path)
merge_table <- vroom::vroom(merge_table_path) %>% dplyr::select(tx_in, tx_out) %>% dplyr::filter(tx_in != tx_out)


############################################################################## #
# ---- 2. Map Transcripts - Truncated - Binned ----

#----------------------------------------------------------------------------- #
## Generate the initial tables ----

tx_tbl <- gtf %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(seqnames, start, end, strand, type, transcript_id, gene_id, ORF_cat, hasCDS) %>%
  dplyr::mutate(type = "transcript")

tr_tbl <- gtf_trunc %>%
  tibble::as_tibble() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(truncation_id = paste0("txcutr_", row_number())) %>%
  dplyr::select(truncation_id, transcript_id, gene_id, seqnames, start, end, strand, type) %>%
  dplyr::mutate(type = "truncation")

bin_tbl <- tr_tbl %>%
  dplyr::filter(!transcript_id %in% unique(merge_table$tx_in)) %>%
  dplyr::mutate(bin_id = paste0("bin_", row_number())) %>%
  dplyr::relocate(bin_id) %>%
  dplyr::mutate(type = "bin")

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

map_tx_tr <- mapping_TxTr(tr_tbl = tr_tbl, overlap_table = overlap_table, n_threads = 64, output_path = map_tx_tr_path)

# Update transcript table with the reference to truncation table
tx_tbl <- tx_tbl %>%
  dplyr::left_join(map_tx_tr, by = "transcript_id") %>%
  dplyr::relocate(transcript_id, truncation_id, gene_id)


#----------------------------------------------------------------------------- #
## Map binned to transcripts and truncated ----
mapping_TrBin <- function(bin_tbl, merge_table, n_threads = 32, output_path = "", overwrite = F){
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

  ### 2. Identify terminal nodes and input nodes
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
  map_tr_bin <- map_tr_bin_s2 %>%
    dplyr::left_join(bin_tbl %>% dplyr::select(truncation_id, bin_id), by = c("output_nodes" = "truncation_id")) %>%
    dplyr::select(truncation_id = input_nodes, bin_id) %>%
    dplyr::bind_rows(map_tr_bin_s1) %>%
    dplyr::mutate(bin_id = factor(bin_id, levels = unique(bin_tbl$bin_id))) %>%
    dplyr::arrange(bin_id) %>%
    dplyr::mutate(bin_id = as.character(bin_id))

  if(output_path != "") map_tr_bin %>% saveRDS(output_path)
  return(map_tr_bin)
}

map_tr_bin <- mapping_TrBin(bin_tbl = bin_tbl, merge_table = merge_table, n_threads = 64, output_path = map_tr_bin_path)

# Update transcript table with the reference to truncation table
tr_tbl <- tr_tbl %>%
  dplyr::left_join(map_tr_bin, by = "truncation_id") %>%
  dplyr::relocate(truncation_id, transcript_id, bin_id, gene_id)
tx_tbl <- tx_tbl %>%
  dplyr::left_join(map_tr_bin, by = "truncation_id") %>%
  dplyr::relocate(transcript_id, truncation_id, bin_id, gene_id)


#----------------------------------------------------------------------------- #
## Generate Open Reading Frame similarity scores - Biostring ----

# Create the ORF table and update the transcript table to include the Open
# Reading Frame sequence and ORF ID
tx_tbl <- tx_tbl %>%
  dplyr::left_join(dplyr::select(sq3_class, isoform, orf_seq = ORF_seq), by = c("transcript_id" = "isoform")) %>%
  dplyr::group_by(orf_seq, gene_id) %>%
  dplyr::mutate(orf_id = paste0("orf_", dplyr::cur_group_id()), .after = gene_id) %>%
  dplyr::ungroup()

orf_tbl <- tx_tbl %>%
  dplyr::distinct(orf_id, orf_seq, gene_id)

#----------------------------------------------------------------------------- #
## Save current tables ----
tx_tbl %>% saveRDS(tx_tbl_path)
tr_tbl %>% saveRDS(tr_tbl_path)
bin_tbl %>% saveRDS(bin_tbl_path)
orf_tbl %>% saveRDS(orf_tbl_path)
