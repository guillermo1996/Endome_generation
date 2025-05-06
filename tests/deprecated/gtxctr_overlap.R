#!/usr/bin/env Rscript

################################################################################
## Mock `snakemake` preamble
################################################################################

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
    input=list(gtf="/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Rescue/sq3.annotated_rescued.gtf"),
    output=list(gtf="/home/grocamora/RytenLab-Research/38-Endome_generation/results/test.gtf.gz",
                fa="/home/grocamora/RytenLab-Research/38-Endome_generation/results/test.fa.gz",
                tsv="/home/grocamora/RytenLab-Research/38-Endome_generation/results/test.tsv"),
    params=list(mergeDist="200", genome="hg38"),
    wildcards=list(width="500"),
    threads=1
  )
}

################################################################################
## Libraries and Parameters
################################################################################

library(gtxcutr)
library(BSgenome)
library(GenomicFeatures)
library(BiocParallel)
library(magrittr)

## convert arguments
maxTxLength <- as.integer(snakemake@wildcards$width)
minDistance <- as.integer(snakemake@params$mergeDist)

## load genome
bsg <- getBSgenome(snakemake@params$genome)

## set cores
# BiocParallel::register(BiocParallel::MulticoreParam(snakemake@threads, progressbar = T))
BiocParallel::register(BiocParallel::SerialParam(progressbar = T))

################################################################################
## Load Data, Truncate, and Export
################################################################################

txdb <- makeTxDbFromGFF(file=snakemake@input$gtf, organism=organism(bsg))
txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
# seqlevels(txdb) <- "chr1"
seqlevelsStyle(txdb) <- "UCSC"

grlExons <- exonsBy(txdb, use.names=TRUE)
mcols(grlExons) <- NULL  # Remove metadata
dfTxGene <- AnnotationDbi::select(txdb, keys=names(grlExons), keytype="TXNAME", columns="GENEID")
mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)

message("Truncating transcripts...")
BPPARAM <- BiocParallel::bpparam()
clipped <- gtxcutr:::.clipTranscript_modded2(grlExons, maxTxLength = maxTxLength, BPPARAM = BPPARAM)

################################################################################
## Prepare the test data
################################################################################
test_txs <- levels(mcols(clipped)[["transcript_id"]])[1:100]

test_data <- clipped %>% plyranges::filter(transcript_id %in% test_txs)# %>% plyranges::mutate(transcript_id = as.character(transcript_id))
test_list <- S4Vectors::split(test_data, mcols(test_data)[["transcript_id"]])


################################################################################
## Default method

# It removes the first findings, leaving only the last one found to be different
overlaps <- findOverlaps(test_list, type = "equal",
                         ignore.strand=FALSE,
                         drop.self=TRUE, drop.redundant=TRUE)

## get duplicate indices
duplicates <- unique(queryHits(overlaps))
if (length(duplicates) > 0) {
  test_list_out <- test_list[-duplicates]
}

ref_out <- test_list_out
ref_df <- unlist(ref_out)
df_overlaps <- overlaps %>% as_tibble()

################################################################################
## My method 1
test_group <- test_data %>%
  tibble::as_tibble() %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::group_split()

tx_duplicated <- rev(test_group) %>% lapply(function(x) x %>% dplyr::select(-transcript_id)) %>% duplicated()
mod_df1 <-test_group[!rev(tx_duplicated)] %>% dplyr::bind_rows() %>% GRanges()
names(mod_df1) <- mcols(mod_df1)[["transcript_id"]]

################################################################################
## My method 2
keep_txs <- test_data %>%
  tibble::as_tibble() %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarise(signature = paste(paste0(seqnames, ":", start, "-", end, ":", strand), collapse = ";")) %>%
  dplyr::group_by(signature) %>%
  dplyr::slice_tail(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::pull(transcript_id) %>%
  sort() %>%
  as.character()

mod_df2 <- test_data %>% plyranges::filter(transcript_id %in% keep_txs)
names(mod_df2) <- mcols(mod_df2)[["transcript_id"]]

################################################################################
## Official benchmark
################################################################################
test_txs <- levels(mcols(clipped)[["transcript_id"]])[1:10000]
test_data <- clipped# %>% plyranges::filter(transcript_id %in% test_txs)
grlC_clipped <- S4Vectors::split(test_data, mcols(test_data)["transcript_id"])

ref_method <- function(grlC_clipped){
  overlaps <- findOverlaps(grlC_clipped, type = "equal",
                           ignore.strand=FALSE,
                           drop.self=TRUE, drop.redundant=TRUE)

  ## get duplicate indices
  duplicates <- unique(queryHits(overlaps))
  if (length(duplicates) > 0) {
    grlC_clipped <- grlC_clipped[-duplicates]
  }

  return(grlC_clipped)
}

mod_method1 <- function(test_data){
  test_group <- test_data %>%
    tibble::as_tibble() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::group_split()

  tx_duplicated <- rev(test_group) %>% lapply(function(x) x %>% dplyr::select(-transcript_id)) %>% duplicated()
  mod_df1 <-test_group[!rev(tx_duplicated)] %>% dplyr::bind_rows() %>% GRanges()

  return(mod_df1)
}

mod_method2 <- function(test_data){
  keep_txs <- test_data %>%
    tibble::as_tibble() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(signature = paste(paste0(seqnames, ":", start, "-", end, ":", strand), collapse = ";")) %>%
    dplyr::group_by(signature) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::pull(transcript_id) %>%
    sort() %>%
    as.character()

  mod_df2 <- test_data %>% plyranges::filter(transcript_id %in% keep_txs)
  return(mod_df2)
}

mod_method3 <- function(test_data){
  keep_txs <- test_data %>%
    tibble::as_tibble() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(signature = paste(paste0(seqnames, ":", start, "-", end, ":", strand), collapse = ";")) %>%
    dplyr::group_by(signature) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::pull(transcript_id)

  mod_df3 <- test_data %>% tibble::as_tibble() %>% dplyr::filter(transcript_id %in% keep_txs) %>% GRanges()
  return(mod_df3)
}

system.time({a <- ref_method(grlC_clipped) %>% unlist()})
system.time({b <- mod_method1(test_data) %>% `names<-`(mcols(.)[["transcript_id"]])})
system.time({c <- mod_method2(test_data) %>% `names<-`(mcols(.)[["transcript_id"]])})
system.time({d <- mod_method3(test_data) %>% `names<-`(mcols(.)[["transcript_id"]])})


library(microbenchmark)

benchmark_results <- microbenchmark(
  ref_method = ref_method(grlC_clipped),
  # mod_method1 = mod_method1(test_data),
  mod_method2 = mod_method2(test_data),
  # mod_method3 = mod_method3(test_data),
  times = 2
); benchmark_results
