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
## Date Created: 2025-03-21
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-03-21)
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
## 0.0 Define Snakemake interactive parameters ----
##
if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots=c(
      input='list',
      output='list',
      params='list',
      wildcards='list',
      threads='numeric'
    )
  )
  snakemake <- Snakemake(
    input=list(gtf="/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Rescue/sq3.annotated_rescued.gtf",
               sq3_class="/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt"),
    # output=list(gtf="/home/grocamora/RytenLab-Research/38-Endome_generation/results/test.gtf.gz",
    #             fa="/home/grocamora/RytenLab-Research/38-Endome_generation/results/test.fa.gz",
    #             tsv="/home/grocamora/RytenLab-Research/38-Endome_generation/results/test.tsv"),
    params=list(gtf_format = "ENSEMBL"),
    # wildcards=list(width="500"),
    threads=1
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(plyranges))
shhh(library(tidyverse))

## Package options - these are some examples and do not modify the script usage.
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)

#----------------------------------------------------------------------------- #
## 0.2 Script Paths ----
main_path <- "/home/grocamora/RytenLab-Research/38-Endome_generation" # here::here()

### Input Paths
input_gtf_path <- snakemake@input$gtf
input_sq3_class_path <- snakemake@input$sq3_class

### Output Paths

#### Create output directory
# dir.create(results_path, showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.3 Script Parameters ----
gtf_format <- snakemake@params$gtf_format %>% tolower()

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
# ---- 1. Load the GTF & Classification ----
input_gtf <- rtracklayer::import.gff(input_gtf_path)
input_sq3_class <- readr::read_delim(input_sq3_class_path)


#----------------------------------------------------------------------------- #
## Filter protein coding transcripts ----
pc_txs <- input_sq3_class %>% dplyr::filter(coding == "coding") %>% dplyr::pull(isoform)

#----------------------------------------------------------------------------- #
## Filter the input GTF for protein coding transcripts ----
input_df <- input_gtf %>% tibble::as_tibble()
pc_input_df <- input_df %>% dplyr::filter(transcript_id %in% pc_txs)

lv
system.time({
  r <- pc_input_df %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::group_split() %>%
    .[1:1000] %>%
    BiocParallel::bplapply(
      correct_transcript_gtf,
      BPPARAM = BiocParallel::SerialParam(progressbar = T)
      # BPPARAM = BiocParallel::MulticoreParam(workers = 32, progressbar = T)
    )
})

correct_transcript_gtf <- function(transcript_df){
  if(all(transcript_df$strand == "-")){
    transcript_df <- transcript_df %>%
      dplyr::mutate(start = -1*start,
                    end = -1*end) %>%
      dplyr::mutate(tmp = start, start = end, end = tmp) %>% dplyr::select(-tmp)
  }

  tx_start <- transcript_df %>% dplyr::filter(type == "transcript") %>% dplyr::pull(start)
  tx_end <- transcript_df %>% dplyr::filter(type == "transcript") %>% dplyr::pull(end)

  first_cds_start <- transcript_df %>% dplyr::filter(type == "CDS") %>% dplyr::pull(start) %>% min()
  last_cds_end <- transcript_df %>% dplyr::filter(type == "CDS") %>% dplyr::pull(end) %>% max()

  five_prime_end <- first_cds_start - 1
  five_prime_utr <- transcript_df %>%
    dplyr::filter(type == "exon" & start < five_prime_end) %>%
    dplyr::mutate(end = ifelse(five_prime_end < end, five_prime_end, end)) %>%
    dplyr::mutate(type = "five_prime_utr")

  three_prime_start <- last_cds_end + 1
  three_prime_utr <- transcript_df %>%
    dplyr::filter(type == "exon" & end > three_prime_start) %>%
    dplyr::mutate(start = ifelse(start < three_prime_start, three_prime_start, start)) %>%
    dplyr::mutate(type = "three_prime_utr")

  start_codon <- transcript_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::filter(start == min(start)) %>%
    dplyr::mutate(end = start + 2) %>%
    dplyr::mutate(type = "start_codon")

  stop_codon <- transcript_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::filter(end == max(end)) %>%
    dplyr::mutate(start = end - 2) %>%
    dplyr::mutate(type = "stop_codon")

  coding_sequences <- transcript_df %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::mutate(end = ifelse(end == max(end), end - 3, end))

  corrected_transcript <- dplyr::bind_rows(
    transcript_df %>% dplyr::filter(type %in% c("transcript", "exon")),
    coding_sequences,
    start_codon,
    stop_codon,
    five_prime_utr,
    three_prime_utr
  )

  if(all(transcript_df$strand == "-")){
    corrected_transcript <- corrected_transcript %>%
      dplyr::mutate(start = -1*start,
                    end = -1*end) %>%
      dplyr::mutate(tmp = start, start = end, end = tmp) %>% dplyr::select(-tmp)
  }

  return(corrected_transcript)
}

transcript_df <- pc_input_df %>% dplyr::filter(grepl("ENST00000000233", transcript_id))
lv <- ensembl_gtf %>% plyranges::filter(transcript_id == "ENST00000000233") %>% plyranges::select(type, gene_id)
