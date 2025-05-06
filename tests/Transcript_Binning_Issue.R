library(rtracklayer)
library(tidyverse)
library(igraph)
library(ggtranscript)
library(GenomicRanges)
library(GenomicFeatures)

lv <- vroom::vroom("/home/grocamora/RytenLab-Research/38-Endome_generation/results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.merge.tsv")
og_gtf <- rtracklayer::import("/home/grocamora/RytenLab-Research/38-Endome_generation/results/ORF_Filter/sq3.annotated_orf.filter.gtf")
tx_gtf <- rtracklayer::import("/home/grocamora/RytenLab-Research/38-Endome_generation/results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.gtf")

g <- igraph::graph_from_data_frame(lv %>% dplyr::select(tx_in, tx_out), directed = T)
g <- simplify(g)

all_nodes <- V(g)$name
is_terminal <- degree(g, mode = "out") == 0
terminal_nodes <- all_nodes[is_terminal]



################################################################################
################################################################################
# Find nodes with the issue
bkup_path <- "tests/Transcript_Binning_Issue.rds"
if(file_test("-f", bkup_path)){
  results <- readRDS(bkup_path)
}else{
  results <- BiocParallel::bplapply(terminal_nodes, BPPARAM = BiocParallel::SerialParam(progressbar = T), function(node){
    input_nodes <- subcomponent(g, node, mode = "in")$name

    if(length(input_nodes) > 2){
      input_txs <- tx_gtf %>% plyranges::filter(transcript_id %in% input_nodes, type == "transcript")
      # if(node == terminal_nodes[[92]]) start(input_txs[5]) <- 17138050

      if(sum(queryHits(findOverlaps(input_txs, drop.redundant = T)) == 1) < length(input_txs)){
        return(input_nodes[[1]])
      }
    }

    return("")
  })

  results %>% saveRDS(bkup_path)
}

which(results != "") %>% length
which(results != "") %>% length / terminal_nodes %>% length

################################################################################
################################################################################
# Example of Binned transcripts (no CDS)
terminal_node <- terminal_nodes[[4760]]
input_nodes <- subcomponent(g, terminal_node, mode = "in")$name

og_transcripts <- og_gtf %>% plyranges::filter(transcript_id %in% input_nodes, type == "exon") %>% tibble::as_tibble() %>% dplyr::mutate(status = "Original")
tx_transcripts <- gtf %>% plyranges::filter(transcript_id %in% input_nodes, type == "exon") %>% tibble::as_tibble() %>% dplyr::mutate(status = "Truncated")
merged_transcripts <- dplyr::bind_rows(og_transcripts, tx_transcripts)

merged_rescaled <- shorten_gaps(
  merged_transcripts,
  to_intron(merged_transcripts, c("transcript_id", "status")),
  group_var = c("transcript_id", "status")
) %>%
  dplyr::mutate(transcript_id = forcats::fct_rev(forcats::fct_inorder(transcript_id)))

tx_labels <- levels(merged_rescaled$transcript_id)
tx_alt_labels <- rev(LETTERS[1:length(tx_labels)])
merged_rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(xstart = start, xend = end, y = as.numeric(transcript_id), fill = status)) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_rescaled, type == "exon", status == "Original")) +
  geom_half_range(range.orientation = "bottom", data = dplyr::filter(merged_rescaled, type == "exon", status == "Truncated")) +
  scale_fill_manual(name = "Transcript Status", values = c("Original" = "gray", "Truncated" = "orange")) +
  geom_intron(data = dplyr::filter(merged_rescaled, type == "intron", status == "Original"), aes(strand = strand)) +
  scale_y_continuous(breaks = 1:length(tx_labels), labels = tx_labels, sec.axis = sec_axis(~., labels = tx_alt_labels, breaks = 1:length(tx_alt_labels))) +
  ggtitle(paste0("Binned transcripts to ", terminal_node)) +
  labs(y = "Transcript ID", x = "Rescaled Genomic Coordinates") +
  theme_bw()

merged_transcripts %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_id, fill = status)) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_transcripts, status == "Original")) +
  geom_half_range(range.orientation = "bottom", data = dplyr::filter(merged_transcripts, status == "Truncated")) +
  scale_fill_manual(name = "Transcript Status", values = c("Original" = "gray", "Truncated" = "orange")) +
  geom_intron(data = to_intron(merged_transcripts, c("transcript_id", "status")), aes(strand = strand)) +
  scale_y_discrete(limits = rev) +
  theme_bw()


################################################################################
################################################################################
# Example of Binned transcripts (with CDS)
terminal_node <- terminal_nodes[[4760]]
input_nodes <- subcomponent(g, terminal_node, mode = "in")$name

og_exons <- og_gtf %>% plyranges::filter(transcript_id %in% input_nodes, type == "exon") %>% tibble::as_tibble() %>% dplyr::mutate(status = "UTRs")
og_cds <- og_gtf %>% plyranges::filter(transcript_id %in% input_nodes, type == "CDS") %>% tibble::as_tibble() %>% dplyr::mutate(status = "CDSs")
tx_transcripts <- gtf %>% plyranges::filter(transcript_id %in% input_nodes, type == "exon") %>% tibble::as_tibble() %>% dplyr::mutate(status = "Truncated")
merged_types <- dplyr::bind_rows(og_exons, og_cds, tx_transcripts)

merged_rescaled <- shorten_gaps(
  merged_types,
  to_intron(merged_types, group_var = c("transcript_id", "status")),
  group_var = c("transcript_id", "status")
) %>%
  dplyr::mutate(transcript_id = forcats::fct_rev(forcats::fct_inorder(transcript_id)))

tx_labels <- levels(merged_rescaled$transcript_id)
tx_alt_labels <- rev(LETTERS[1:length(tx_labels)])
merged_rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(xstart = start, xend = end, y = as.numeric(transcript_id), fill = status)) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_rescaled, type == "exon", status == "UTRs"), height = 0.15) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_rescaled, type == "CDS", status == "CDSs")) +
  geom_half_range(range.orientation = "bottom", data = dplyr::filter(merged_rescaled, type == "exon", status == "Truncated"), height = 0.15) +
  scale_fill_manual(name = "Transcript Status", values = c("UTRs" = "white", "CDSs"= "#ac39ff", "Truncated" = "orange")) +
  geom_intron(data = dplyr::filter(merged_rescaled, type == "intron", status == "UTRs"), aes(strand = strand)) +
  scale_y_continuous(breaks = 1:length(tx_labels), labels = tx_labels, sec.axis = sec_axis(~., labels = tx_alt_labels, breaks = 1:length(tx_alt_labels))) +
  ggtitle(paste0("Binned transcripts to ", terminal_node, " - with CDS information")) +
  labs(y = "Transcript ID", x = "Rescaled Genomic Coordinates") +
  theme_bw()


merged_types %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_id, fill = status)) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_types, type == "exon", status == "UTRs"), height = 0.25) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_types, type == "CDS", status == "CDSs")) +
  geom_half_range(range.orientation = "bottom", data = dplyr::filter(merged_types, type == "exon", status == "Truncated"), height = 0.25) +
  scale_fill_manual(name = "Transcript Status", values = c("UTRs" = "gray", "CDSs"= "#ac39ff", "Truncated" = "orange")) +
  geom_intron(data = to_intron(dplyr::filter(merged_types, status == "UTRs"), "transcript_id"), aes(strand = strand)) +
  theme_bw()





################################################################################
################################################################################
# Test with CDS
sod1_annotation <- og_gtf %>% tibble::as_tibble() %>% dplyr::filter(transcript_id %in% input_nodes)
sod1_exons <- sod1_annotation %>% dplyr::filter(type == "exon")
sod1_cds <- sod1_annotation %>% dplyr::filter(type == "CDS")

sod1_exons %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
  geom_range(fill = "white", height = 0.25) +
  geom_range(data = sod1_cds) +
  geom_intron(data = to_intron(sod1_exons, "transcript_id"), aes(strand = strand), arrow.min.intron.length = 500)
