#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcriptome truncation with txendcutr
##
## Aim: Truncate the ORF-filtered transcriptome to its 3' or 5' ends and export
## the ENDome GTF, its sequences and the merge table of its bins
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S05a-txendcutr
## Date Created: 2025-03-13
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.3 (2026-10-29)
##
## _________________________________________________
##
## Notes:
##    - Width "full" keeps the whole transcripts. The TxDb gets the same
##      metadata as txendcutr::truncateTxome and the overlap table is empty.
##    - Width "UTR" keeps the 3' or 5' UTR of each transcript (transcripts
##      without it are dropped). "UTR_exact" keeps the same features, but only
##      identical UTRs share a bin (no merge distance).
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.3 (2026-10-29): Widths "UTR" and "UTR_exact", from T05a-txendcutr_UTR.R.
##    - v1.2 (2026-10-28): Standardization of the script.
##    - v1.1 (2026-10-27): Width "full" (no truncation).
##    - v1.0 (2025-03-13): Initial version.
## _________________________________________________

############################################################################## #
# ---- 0. Setup ----

#----------------------------------------------------------------------------- #
## 0.0 Define Snakemake interactive parameters ----
if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots = c(
      input = 'list',
      output = 'list',
      params = 'list',
      wildcards = 'list',
      log = 'list',
      threads = 'numeric',
      scriptdir = 'character')
  )

  rs <- "results/Wood.control"
  prefix <- "Wood.control.st_ref"
  orf_filter <- "pc_novel"
  width <- "500"
  txEnd <- "3p"
  endome <- sprintf("%s.%s.w%s.%s", prefix, orf_filter, width, txEnd)
  dir.create(file.path(rs, "temp"), showWarnings = FALSE)

  snakemake <- Snakemake(
    input = list(
      gtf = file.path(rs, "04-ORF_Identification/ORF_Filtration", sprintf("%s.%s.orf_filter.gtf", prefix, orf_filter))
    ),
    output = list(
      gtf = file.path(rs, "temp", sprintf("%s.txendcutr.gtf", endome)),
      fa = file.path(rs, "temp", sprintf("%s.txendcutr.fa.gz", endome)),
      transcript_overlap = file.path(rs, "temp", sprintf("%s.txendcutr.overlaps.tsv", endome)),
      merge_table = file.path(rs, "temp", sprintf("%s.txendcutr.merge.tsv", endome))
    ),
    params = list(
      mergeDist = "200",
      genome = "hg38"
    ),
    wildcards = list(
      prefix = prefix,
      orf_filter = orf_filter,
      width = width,
      txEnd = txEnd
    ),
    threads = 8,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages
shhh({
  library(txendcutr)
  library(BSgenome)
  library(GenomicFeatures)
  library(GenomeInfoDb)
  library(BiocParallel)
  library(txdbmaker)
})

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
smk_inputs <- snakemake@input
smk_outputs <- snakemake@output
params <- snakemake@params
wc <- snakemake@wildcards

### Width can be any of the following: a number (last/first nt kept), "full"
### (whole transcripts) or the 3'/5' UTR, binned by the merge distance (UTR) or
### exactly (UTR_exact)
truncation_mode <- if (wc$width == "full") "full" else if (wc$width %in% c("UTR", "UTR_exact")) "utr" else "window"
maxTxLength <- if (truncation_mode == "window") as.integer(wc$width) else NA_integer_
if (truncation_mode == "window" && is.na(maxTxLength)) stop("Unknown width: ", wc$width, call. = FALSE)

txEnd <- wc$txEnd
merge_ends <- wc$width != "UTR_exact"
minDistance <- as.integer(params$mergeDist)

bsg <- getBSgenome(params$genome)
BiocParallel::register(BiocParallel::MulticoreParam(snakemake@threads, progressbar = FALSE))

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----

#' Transcriptome reduced to the 3' or 5' UTR of each transcript
#'
#' @param txdb TxDb with CDS (the S04b GTF has them).
#' @param txEnd "3p" or "5p".
#' @param overlapFile Path of the overlap table, same columns as txendcutr::truncateTxome.
#'
#' @returns TxDb with the txendcutr::truncateTxome metadata. Transcripts without
#'   that UTR are dropped, and identical UTRs of a gene are kept once.
utrTxome <- function(txdb, txEnd, overlapFile){
  grl <- switch(txEnd,
    "3p" = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE),
    "5p" = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE))
  txs <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"))
  tx_gene <- setNames(as.character(unlist(txs$gene_id)), txs$tx_name)
  message(sprintf("  %s of %s transcripts have a %s UTR",
    format(length(grl), big.mark = ","), format(length(txs), big.mark = ","), txEnd))

  ### Identical UTRs of the same gene, as txendcutr::truncateTxome
  hits <- as.data.frame(findOverlaps(grl, type = "equal", drop.self = TRUE, drop.redundant = TRUE))
  hits$queryTx <- names(grl)[hits$queryHits]
  hits$subjectTx <- names(grl)[hits$subjectHits]
  hits$queryGene <- unname(tx_gene[hits$queryTx])
  hits$subjectGene <- unname(tx_gene[hits$subjectTx])
  hits <- hits[hits$queryGene == hits$subjectGene, ]
  write.table(hits, overlapFile, sep = "\t", row.names = FALSE, quote = FALSE)
  if (nrow(hits) > 0) grl <- grl[-unique(hits$queryHits)]
  message(sprintf("  Removed %s identical UTRs", format(length(unique(hits$queryHits)), big.mark = ",")))

  ### One exon_id per exon line, as in the full-length branch
  exons <- unlist(grl, use.names = FALSE)
  exons <- sort(exons) 
  mcols(exons) <- S4Vectors::DataFrame(
    type = rep("exon", length(exons)),
    transcript_id = rep(names(grl), lengths(grl)),
    exon_id = paste0("exon_", seq_along(exons)))
  tx_ranges <- unlist(range(grl), use.names = FALSE)
  mcols(tx_ranges) <- S4Vectors::DataFrame(
    type = rep("transcript", length(grl)),
    transcript_id = names(grl),
    exon_id = rep(NA_character_, length(grl)))
  gr <- c(tx_ranges, exons)
  gr$gene_id <- unname(tx_gene[gr$transcript_id])

  df_metadata <- data.frame(
    name = c("Truncated by", "Maximum Transcript Length", "Truncation End"),
    value = c("txendcutr", "UTR", c("3p" = "3prime", "5p" = "5prime")[[txEnd]]))
  suppressWarnings(txdbmaker::makeTxDbFromGRanges(gr, taxonomyId = taxonomyId(txdb), metadata = df_metadata))
}

############################################################################## #
# ---- 1. Truncate the transcriptome ----
if(truncation_mode == "full"){
 message("Loading the full-length transcriptome (no truncation)...")

  ### One exon_id per exon line, as in txendcutr::truncateTxome. Otherwise
  ### exportGTF writes an exon shared by several transcripts several times.
  gtf <- rtracklayer::import(smk_inputs$gtf, format = "gtf")
  gtf <- keepStandardChromosomes(gtf, pruning.mode = "coarse")
  seqlevelsStyle(gtf) <- "UCSC"
  gtf$exon_id <- ifelse(gtf$type == "exon", paste0("exon_", seq_along(gtf)), NA_character_)

  ### Same metadata as txendcutr::truncateTxome, read by exportMergeTable
  df_metadata <- data.frame(
    name = c("Truncated by", "Maximum Transcript Length", "Truncation End"),
    value = c("none", "full", ifelse(txEnd == "5p", "5prime", "3prime"))
  )
  txdb_result <- suppressWarnings(makeTxDbFromGRanges(gtf, metadata = df_metadata))

  ### Empty overlap table with the columns of txendcutr::truncateTxome
  write.table(
    data.frame(queryHits = integer(), subjectHits = integer(), queryTx = character(),
               subjectTx = character(), queryGene = character(), subjectGene = character()),
    smk_outputs$transcript_overlap, sep = "\t", row.names = FALSE, quote = FALSE)
}else{
  txdb <- makeTxDbFromGFF(file = smk_inputs$gtf, organism = organism(bsg))
  txdb <- keepStandardChromosomes(txdb, pruning.mode = "coarse")
  seqlevelsStyle(txdb) <- "UCSC"

  if (truncation_mode == "window") {
    message("Truncating the transcriptome to ", maxTxLength, " nt (", txEnd, ")...")
    txdb_result <- truncateTxome(txdb, maxTxLength = maxTxLength, overlapFile = smk_outputs$transcript_overlap, txEnd = txEnd)
  } else {
    message("Reducing the transcriptome to its ", txEnd, " UTRs...")
    txdb_result <- utrTxome(txdb, txEnd, smk_outputs$transcript_overlap)
  }
}

############################################################################## #
# ---- 2. Write the outputs ----
message("Exporting the GTF...")
exportGTF(txdb_result, smk_outputs$gtf)

message("Exporting the FASTA...")
exportFASTA(txdb_result, bsg, smk_outputs$fa)

message("Exporting the merge table...")
if (merge_ends) {
  exportMergeTable(txdb_result, smk_outputs$merge_table, minDistance = minDistance)
} else {
  ### No merge distance: each transcript is its own bin, same columns as exportMergeTable
  txs <- GenomicFeatures::transcripts(txdb_result, columns = c("tx_name", "gene_id"))
  merge_df <- data.frame(tx_in = txs$tx_name, tx_out = txs$tx_name, gene_out = as.character(unlist(txs$gene_id)))
  write.table(merge_df[order(merge_df$tx_in), ], smk_outputs$merge_table, sep = "\t", row.names = FALSE, quote = FALSE)
}

message("Done.")
