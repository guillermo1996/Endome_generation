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
  ## dataset.group.merge_method folder created by the `test_data` rule, and the
  ## orf_filter preset and truncation (width, txEnd) to test. Outputs use the
  ## "interactive" prefix so they never overwrite the linked results.
  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  test_orf_filter <- "ref_pc"
  test_width <- "500"
  test_txEnd <- "3p"
  test_out <- file.path(test_dir, paste0("interactive.", test_orf_filter, ".txendcutr.w", test_width, ".", test_txEnd))
  snakemake <- Snakemake(
      input=list(gtf=file.path(test_dir, paste0("test.", test_orf_filter, ".orf_filter.gtf"))),
      output=list(gtf=paste0(test_out, ".gtf"),
                  fa=paste0(test_out, ".fa.gz"),
                  transcript_overlap=paste0(test_out, ".overlaps.tsv"),
                  merge_table=paste0(test_out, ".merge.tsv")),
      params=list(mergeDist="200", genome="hg38"),
      wildcards=list(width=test_width, txEnd=test_txEnd),
      log = list(paste0(test_out, ".log")),
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
