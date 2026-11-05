#!/usr/bin/env Rscript

## _________________________________________________
##
## Import a pre-truncated GTF
##
## Aim: Use an ENDome built outside the pipeline in steps 06-07. Writes the
## same GTF, FASTA and merge table as S05a
##
## Author: Guillermo Rocamora Pérez
## Project: ENDome generation - S05e-Import_Truncated_GTF
## Date Created: 2026-10-26
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
## Email: guillermorocamora@gmail.com
##
## Latest Version: v1.2 (2026-10-28)
##
## _________________________________________________
##
## Notes:
##    - The GTF and FASTA are exported with txendcutr, so the output has the
##      same features and attributes as a txendcutr ENDome.
##    - Without a merge table (`pre-truncated_merge_tsv`), only identical
##      transcripts of a gene share a bin.
##    - Developed with the assistance of Claude Code (Opus 5.5).
##
## Changelog:
##    - v1.2 (2026-10-28): Standardization of the script.
##    - v1.1 (2026-10-26): GTF and FASTA exported with txendcutr.
##    - v1.0 (2026-10-26): Initial version.
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

  rs <- "results/gencode_3UTR"
  prefix <- "gencode_3UTR.none.ext"
  orf_filter <- "ext"
  width <- "UTR"
  txEnd <- "3p"
  endome <- sprintf("%s.%s.w%s.%s", prefix, orf_filter, width, txEnd)
  dir.create(file.path(rs, "temp"), showWarnings = FALSE, recursive = TRUE)

  snakemake <- Snakemake(
    input = list(
      gtf = "data/custom_inputs/gencode.v48.pc_3UTR.gtf",
      merge_tsv = character(0)
    ),
    output = list(
      gtf = file.path(rs, "temp", sprintf("%s.txendcutr.gtf", endome)),
      fa = file.path(rs, "temp", sprintf("%s.txendcutr.fa.gz", endome)),
      merge_table = file.path(rs, "temp", sprintf("%s.txendcutr.merge.tsv", endome))
    ),
    params = list(
      genome = "hg38"
    ),
    wildcards = list(
      prefix = prefix,
      orf_filter = orf_filter,
      width = width,
      txEnd = txEnd
    ),
    threads = 1,
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
  library(txdbmaker)
})

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
smk_inputs <- snakemake@input
smk_outputs <- snakemake@output
params <- snakemake@params
wc <- snakemake@wildcards

endome_name <- paste(wc$prefix, wc$orf_filter, paste0("w", wc$width), wc$txEnd, sep = ".")
tx_end <- c("3p" = "3prime", "5p" = "5prime")[[wc$txEnd]]

### Shorter transcripts have no k-mer in the kallisto index
kallisto_k <- 31

############################################################################## #
# ---- 1. Load the pre-truncated GTF ----
message("Loading the genome ", params$genome, "...")
bsg <- getBSgenome(params$genome)

message("Reading the pre-truncated GTF of ", endome_name, "...")

### One exon_id per exon line, as in txendcutr::truncateTxome
gtf <- rtracklayer::import(smk_inputs$gtf, format = "gtf")
gtf <- keepStandardChromosomes(gtf, pruning.mode = "coarse")
seqlevelsStyle(gtf) <- "UCSC"

### Same chromosomes as the assembled ENDomes (step 02 drops chrM)
gtf <- dropSeqlevels(gtf, intersect("chrM", seqlevels(gtf)), pruning.mode = "coarse")

gtf$exon_id <- ifelse(gtf$type == "exon", paste0("exon_", seq_along(gtf)), NA_character_)

### Same metadata as txendcutr::truncateTxome
df_metadata <- data.frame(
  name = c("Truncated by", "Maximum Transcript Length", "Truncation End"),
  value = c("external", wc$width, tx_end)
)
txdb <- suppressWarnings(makeTxDbFromGRanges(gtf, metadata = df_metadata))

message(sprintf("  %s transcripts | %s genes",
  format(length(transcripts(txdb)), big.mark = ","),
  format(length(genes(txdb)), big.mark = ",")))

############################################################################## #
# ---- 2. Write the outputs ----

#----------------------------------------------------------------------------- #
## 2.1 ENDome GTF ----
message("Exporting the GTF...")
txendcutr::exportGTF(txdb, smk_outputs$gtf)

#----------------------------------------------------------------------------- #
## 2.2 Transcript sequences ----
message("Exporting the FASTA...")
txendcutr::exportFASTA(txdb, bsg, smk_outputs$fa)

tx_lengths <- GenomicFeatures::transcriptLengths(txdb)$tx_len
n_short <- sum(tx_lengths < kallisto_k)
if (n_short > 0) {
  message(sprintf("  %s transcripts are shorter than %s nt and cannot be quantified by kallisto",
    format(n_short, big.mark = ","), kallisto_k))
}

#----------------------------------------------------------------------------- #
## 2.3 Merge table ----
if (length(smk_inputs$merge_tsv) > 0) {
  message("Copying the merge table ", smk_inputs$merge_tsv, " ...")
  file.copy(smk_inputs$merge_tsv, smk_outputs$merge_table, overwrite = TRUE)
} else {
  message("Building the merge table: identical transcripts of a gene share a bin...")
  exons <- unlist(GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE))
  exon_keys <- paste0(seqnames(exons), ":", start(exons), "-", end(exons), ":", strand(exons))
  tx_chain <- tapply(exon_keys, names(exons), function(keys) paste(sort(keys), collapse = ","))

  txs <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"))
  tx_gene <- setNames(as.character(unlist(txs$gene_id)), txs$tx_name)

  ### Same gene and exon chain: merged into the first transcript id
  merge_df <- data.frame(tx_in = names(tx_chain), chain = unname(tx_chain))
  merge_df$gene_out <- unname(tx_gene[merge_df$tx_in])
  merge_df$tx_out <- ave(merge_df$tx_in, merge_df$gene_out, merge_df$chain, FUN = min)

  ### Same format as txendcutr::exportMergeTable
  utils::write.table(merge_df[, c("tx_in", "tx_out", "gene_out")], smk_outputs$merge_table,
    sep = "\t", row.names = FALSE, quote = FALSE)
}

merge_df <- utils::read.delim(smk_outputs$merge_table)
message(sprintf("  %s transcripts -> %s bins",
  format(nrow(merge_df), big.mark = ","),
  format(length(unique(merge_df$tx_out)), big.mark = ",")))

message("Done.")
