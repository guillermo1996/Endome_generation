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
      input=list(gtf="homo_sapiens/Homo_sapiens.GRCh38.104.mRNA_ends_found.gtf.gz"),
      output=list(gtf="/fscratch/fanslerm/ensembl.hg38.mRNA_ends_found.txcutr.w500.gtf",
                  fa="/fscratch/fanslerm/ensembl.hg38.mRNA_ends_found.txcutr.w500.fa",
                  tsv="/fscratch/fanslerm/ensembl.hg38.mRNA_ends_found.txcutr.w500.merge.tsv"),
      params=list(mergeDist="200", genome="hg38"),
      wildcards=list(width="500"),
      threads=1
  )
}

################################################################################
## Libraries and Parameters
################################################################################

library(txcutr)
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
BiocParallel::register(BiocParallel::MulticoreParam(snakemake@threads))

################################################################################
## Load Data, Truncate, and Export
################################################################################

txdb <- makeTxDbFromGFF(file=snakemake@input$gtf, organism=organism(bsg))
txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
seqlevelsStyle(txdb) <- "UCSC"

txdb_result <- truncateTxome(txdb, maxTxLength)

exportGTF(txdb_result, snakemake@output$gtf)

exportFASTA(txdb_result, bsg, snakemake@output$fa)

exportMergeTable(txdb_result, snakemake@output$tsv, minDistance=minDistance)
