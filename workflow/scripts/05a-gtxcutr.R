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
      input=list(gtf="~/RytenLab-Research/snakefile-refactor/results_k15/04-ORF_Identification/ORF_Filtration/pigeon.annotated_orf.filter.gtf"),
      output=list(gtf="~/RytenLab-Research/38-Endome_generation/results/gtxcutr_test.gtf",
                  fa="~/RytenLab-Research/38-Endome_generation/results/gtxcutr_test.fa",
                  transcript_overlap="~/RytenLab-Research/38-Endome_generation/results/gtxcutr_test.overlaps.tsv",
                  merge_table="~/RytenLab-Research/38-Endome_generation/results/gtxcutr_test.merge.tsv"),
      params=list(mergeDist="200", genome="hg38"),
      wildcards=list(width="500", txEnd="3p"),
      log = list("~/RytenLab-Research/38-Endome_generation/results/gtxcutr.log"),
      threads=8
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

require(txdbmaker)

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
txdb_result <- truncateTxome(txdb, maxTxLength = maxTxLength, overlap_path = overlap_path, txEnd = txEnd)

print("Export GTF")
exportGTF(txdb_result, snakemake@output$gtf)

print("Export FASTA")
exportFASTA(txdb_result, bsg, snakemake@output$fa)

print("Export MergeTable")
exportMergeTable(txdb_result, snakemake@output$merge_table, minDistance=minDistance)
