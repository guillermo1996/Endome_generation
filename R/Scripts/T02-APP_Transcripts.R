## _________________________________________________
##
## APP Transcript Visualization
##
## Aim:
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-05-20
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-05-20)
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
shhh(library(ggtranscript))
shhh(library(plyranges))
shhh(library(tidyverse))

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
truncation_direction = "3p"

#----------------------------------------------------------------------------- #
## 0.3 Script Paths ----
main_path <- "~/RytenLab-Research/38-Endome_generation"
script_path <- file.path(main_path, "R")

### Input Paths
gtf_path <- file.path(main_path, "results/ORF_Filter/sq3.annotated_orf.filter.gtf")
sq3_class_path<- file.path(main_path, "results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt")

gtf_path <- file.path(main_path, "results/ORF_Filter/sq3.annotated_orf.filter.gtf")
if(truncation_direction == "3p"){
  gtf_trunc_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.gtf")
  merge_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w500.3p.merge.tsv")
  tx_tbl_path <- file.path(script_path, "Results/Table_mapping_3p/tx_tbl.rds")

  APP_out_path <- file.path(script_path, "figures/APP_transcripts_3p.png")
}else if (truncation_direction == "5p"){
  gtf_trunc_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w300.5p.gtf")
  merge_table_path <- file.path(main_path, "results/gtxcutr/sq3.annotated.gtxcutr.w300.5p.merge.tsv")
  tx_tbl_path <- file.path(script_path, "Results/Table_mapping_5p/tx_tbl.rds")

  APP_out_path <- file.path(script_path, "figures/APP_transcripts_5p.png")
}else{
  stop("No valid truncation direction")
}

dir.create(dirname(APP_out_path), recursive = T, showWarnings = F)


############################################################################## #
# ---- 1. Load all the needed files ----
gtf <- rtracklayer::import(gtf_path)

gtf_trunc <- rtracklayer::import(gtf_trunc_path)
merge_table <- vroom::vroom(merge_table_path)
tx_tbl <- readRDS(tx_tbl_path)


############################################################################## #
# ---- 2. Plot APP transcripts ----

#----------------------------------------------------------------------------- #
## 2.1 Select APP transcripts ----
APP_id = "ENSG00000142192.22"

app_transcripts <- tx_tbl %>% dplyr::filter(gene_id == APP_id) %>% dplyr::pull(transcript_id)
app_transcripts <- c("ENST00000354192.7", "ENST00000348990.9", "ENST00000357903.7", "ENST00000359726.7",
                     # "ENST00000440126.7",
                     "ENST00000439274.6", "ENST00000358918.7", "ENST00000707132.1")

#----------------------------------------------------------------------------- #
## 2.2 Extract and merge Exon info ----
gtf_exons <- gtf %>%
  plyranges::filter(transcript_id %in% app_transcripts, type == "exon") %>%
  tibble::as_tibble() %>% dplyr::mutate(status = "UTRs") %>%
  dplyr::left_join(tx_tbl %>% dplyr::select(transcript_id, truncation_id, bin_id, orf_id), by = "transcript_id")
gtf_cds <- gtf %>%
  plyranges::filter(transcript_id %in% app_transcripts, type == "CDS") %>%
  tibble::as_tibble() %>% dplyr::mutate(status = "CDSs") %>%
  dplyr::left_join(tx_tbl %>% dplyr::select(transcript_id, truncation_id, bin_id, orf_id), by = "transcript_id")
tr_transcripts <- gtf_trunc %>%
  plyranges::filter(transcript_id %in% app_transcripts, type == "exon") %>%
  tibble::as_tibble() %>% dplyr::mutate(status = "Truncated") %>%
  dplyr::left_join(tx_tbl %>% dplyr::select(transcript_id, truncation_id, bin_id, orf_id), by = "transcript_id")

merged_transcripts <- dplyr::bind_rows(gtf_exons, gtf_cds, tr_transcripts)
merged_rescaled <- ggtranscript::shorten_gaps(
  merged_transcripts,
  ggtranscript::to_intron(merged_transcripts, c("transcript_id", "status")),
  group_var = c("transcript_id", "status")
) %>%
  dplyr::mutate(transcript_id = forcats::fct_rev(forcats::fct_inorder(transcript_id)))


#----------------------------------------------------------------------------- #
## 2.3 Graph Results ----
tx_labels <- levels(merged_rescaled$transcript_id)
tx_labels2 <- setNames(paste0(tx_labels, "\n(", rev(LETTERS[1:length(tx_labels)]), ")"), tx_labels)

merged_rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_id, fill = status)) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_rescaled, type == "exon", status == "UTRs"), height = 0.15) +
  geom_half_range(range.orientation = "top", data = dplyr::filter(merged_rescaled, type == "CDS", status == "CDSs")) +
  geom_half_range(range.orientation = "bottom", data = dplyr::filter(merged_rescaled, type == "exon", status == "Truncated"), height = 0.15) +
  scale_fill_manual(name = "Transcript Status", values = c("UTRs" = "gray", "CDSs"= "#ac39ff", "Truncated" = "orange")) +
  geom_intron(data = dplyr::filter(merged_rescaled, type == "intron", status == "UTRs"), aes(strand = strand)) +
  scale_y_discrete(labels = tx_labels2) +
  ggtitle(paste0("Binned transcripts of APP - with CDS information")) +
  labs(y = "Transcript ID", x = "Rescaled Genomic Coordinates") +
  ggh4x::facet_nested(bin_id + truncation_id ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0)
  )

ggsave(APP_out_path, dpi = 300, width = 13, height = 9, bg = "white")

## URGENT: There is an issue with transcript ENST00000440126.7, where its ending
## position does not match that of annotated APP. I was unable to find the
## sample where this transcript was found or any reason as to why did this
## particular transcript had that wrong exon in the end.
tx_issue = "ENST00000440126.7"

