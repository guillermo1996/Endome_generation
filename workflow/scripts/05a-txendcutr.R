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
      log="list",
      threads='numeric'
    )
  )
  snakemake <- Snakemake(
      input=list(gtf="data/test_data/test.pc.orf_filter.gtf"),
      output=list(gtf="data/test_data/txendcutr_test.gtf",
                  fa="data/test_data/txendcutr_test.fa",
                  transcript_overlap="data/test_data/txendcutr_test.overlaps.tsv",
                  merge_table="data/test_data/txendcutr_test.merge.tsv"),
      params=list(mergeDist="200", genome="hg38"),
      wildcards=list(width="500", txEnd="3p"),
      log = list("data/test_data/txendcutr_logs.log"),
      threads=8
  )
}

################################################################################
## Libraries and Parameters
################################################################################
suppressMessages({
  library(txendcutr)
  library(BSgenome)
  library(GenomicFeatures)
  library(GenomeInfoDb)
  library(GenomeInfoDbData)
  library(BiocParallel)
  library(magrittr)
  library(txdbmaker)
})


## convert arguments
maxTxLength <- as.integer(snakemake@wildcards$width)
txEnd <- snakemake@wildcards$txEnd
minDistance <- as.integer(snakemake@params$mergeDist)
overlap_path <- snakemake@output$transcript_overlap

## load genome
bsg <- getBSgenome(snakemake@params$genome)

## set cores
BiocParallel::register(BiocParallel::MulticoreParam(snakemake@threads, progressbar = F))
BPPARAM = bpparam()
# BiocParallel::register(BiocParallel::SerialParam(progressbar = T))

################################################################################
## Load Data, Truncate, and Export
################################################################################
txdb <- makeTxDbFromGFF(file=snakemake@input$gtf, organism=organism(bsg))

txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
seqlevelsStyle(txdb) <- "UCSC"

# overlap_path is optional
# txEnd is optional. Defaults to 3' truncation
txdb_result <- truncateTxome(txdb, maxTxLength = maxTxLength, overlapFile = overlap_path, txEnd = txEnd)

print("Export GTF")
exportGTF(txdb_result, snakemake@output$gtf)

print("Export FASTA")
exportFASTA(txdb_result, bsg, snakemake@output$fa)

print("Export MergeTable")
exportMergeTable(txdb_result, snakemake@output$merge_table, minDistance=minDistance)
