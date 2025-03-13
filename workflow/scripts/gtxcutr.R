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
seqlevelsStyle(txdb) <- "UCSC"

txdb_result <- truncateTxome(txdb, maxTxLength = maxTxLength)

print("Export GTF")
exportGTF(txdb_result, snakemake@output$gtf)

print("Export FASTA")
exportFASTA(txdb_result, bsg, snakemake@output$fa)

print("Export MergeTable")
exportMergeTable(txdb_result, snakemake@output$tsv, minDistance=minDistance)
