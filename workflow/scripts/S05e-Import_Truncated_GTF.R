## _________________________________________________
##
## Import a pre-truncated GTF
##
## Aim: Use an ENDome built outside the pipeline (e.g. a UTR-only transcriptome)
## in steps 06-07 without running steps 03-05. Writes the same three files as
## S05a (txendcutr): the GTF, the transcript sequences and the merge table
##
## Project: ENDome generation - Custom entry points
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-26
##
## Latest Version: v1.1 (2026-10-26)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##   - The GTF is loaded and exported exactly as in S05a, with txendcutr's
##     exportGTF / exportFASTA. The output has the same
##     gene, transcript and exon lines and attributes as a txendcutr ENDome.
##     Transcript and gene lines are built from the exons when the input has
##     none.
##   - The TxDb gets the same metadata as txendcutr::truncateTxome.
##   - Without a merge table (`pre-truncated_merge_tsv`), only identical
##     transcripts are merged: those of a gene with the same exon chain share a
##     bin. scUTRquant counts each `tx_out` as one feature, so identical
##     sequences must share a bin.
##
## Changelog:
##    - v1.1 (2026-10-26): Export the GTF and FASTA with txendcutr, so the
##      outputs match S05a.
##    - v1.0 (2026-10-26): Initial version.
##
## Contact: guillermorocamora@gmail.com
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
  test_dir <- "data/test_data/gencode_3UTR.none.ext"
  test_endome <- "ext.wUTR.3p"
  dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)

  snakemake <- Snakemake(
    input = list(
      gtf = "data/custom_inputs/gencode.v48.pc_3UTR.gtf",
      merge_tsv = character(0)
    ),
    output = list(
      gtf = file.path(test_dir, paste0("interactive.", test_endome, ".txendcutr.gtf")),
      fa = file.path(test_dir, paste0("interactive.", test_endome, ".txendcutr.fa.gz")),
      merge_table = file.path(test_dir, paste0("interactive.", test_endome, ".txendcutr.merge.tsv"))
    ),
    params = list(
      genome = "hg38"
    ),
    wildcards = list(
      prefix = basename(test_dir),
      orf_filter = "ext",
      width = "UTR",
      txEnd = "3p"
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
smk_params <- snakemake@params
wc <- snakemake@wildcards

endome_name <- paste(wc$prefix, wc$orf_filter, paste0("w", wc$width), wc$txEnd, sep = ".")

### Truncation end as txendcutr names it ("3prime" / "5prime")
tx_end <- c("3p" = "3prime", "5p" = "5prime")[[wc$txEnd]]

### Shorter transcripts have no k-mer in the kallisto index (k = 31)
kallisto_k <- 31

############################################################################## #
# ---- 1. Load the pre-truncated GTF ----
message("Loading the genome ", smk_params$genome, "...")
bsg <- getBSgenome(smk_params$genome)

message("Reading the pre-truncated GTF of ", endome_name, "...")
## Same metadata as txendcutr::truncateTxome
df_metadata <- data.frame(
  name = c("Truncated by", "Maximum Transcript Length", "Truncation End"),
  value = c("external", wc$width, tx_end)
)

## One exon_id per exon line, so an exon shared by several transcripts is kept
## once per transcript (as in txendcutr::truncateTxome). Otherwise exportGTF
## writes it several times in each transcript.
gtf <- rtracklayer::import(smk_inputs$gtf, format = "gtf")
gtf <- keepStandardChromosomes(gtf, pruning.mode = "coarse")
seqlevelsStyle(gtf) <- "UCSC"
gtf$exon_id <- ifelse(gtf$type == "exon", paste0("exon_", seq_along(gtf)), NA_character_)

txdb <- suppressWarnings(makeTxDbFromGRanges(gtf, metadata = df_metadata))

message(sprintf("  %s transcripts | %s genes",
  format(length(transcripts(txdb)), big.mark = ","),
  format(length(genes(txdb)), big.mark = ",")))

############################################################################## #
# ---- 2. Write the outputs ----
for (out_file in unlist(smk_outputs)) {
  if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
}

#----------------------------------------------------------------------------- #
## 2.1 The ENDome GTF ----
message("Exporting the GTF...")
txendcutr::exportGTF(txdb, smk_outputs$gtf)

#----------------------------------------------------------------------------- #
## 2.2 The transcript sequences for the kallisto index ----
message("Exporting the FASTA...")
txendcutr::exportFASTA(txdb, bsg, smk_outputs$fa)

tx_lengths <- GenomicFeatures::transcriptLengths(txdb)$tx_len
n_short <- sum(tx_lengths < kallisto_k)
if (n_short > 0) {
  message(sprintf("  %s transcripts are shorter than %s nt and cannot be quantified by kallisto",
    format(n_short, big.mark = ","), kallisto_k))
}

#----------------------------------------------------------------------------- #
## 2.3 The merge table ----
if (length(smk_inputs$merge_tsv) > 0) {
  message("Copying the merge table ", smk_inputs$merge_tsv, " ...")
  file.copy(smk_inputs$merge_tsv, smk_outputs$merge_table, overwrite = TRUE)
} else {
  message("Building the merge table: identical transcripts of a gene share a bin...")
  ## Exon chain of each transcript: its exons as "chr:start-end:strand"
  exons <- unlist(GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE))
  exon_keys <- paste0(seqnames(exons), ":", start(exons), "-", end(exons), ":", strand(exons))
  tx_chain <- tapply(exon_keys, names(exons), function(keys) paste(sort(keys), collapse = ","))

  ## Gene of each transcript
  txs <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"))
  tx_gene <- setNames(as.character(unlist(txs$gene_id)), txs$tx_name)

  ## Same gene and exon chain: merged into the first transcript id, as txendcutr
  ## does for transcripts with equal ends
  merge_df <- data.frame(tx_in = names(tx_chain), chain = unname(tx_chain))
  merge_df$gene_out <- unname(tx_gene[merge_df$tx_in])
  merge_df$tx_out <- ave(merge_df$tx_in, merge_df$gene_out, merge_df$chain, FUN = min)

  ## Same format as txendcutr::exportMergeTable
  utils::write.table(merge_df[, c("tx_in", "tx_out", "gene_out")], smk_outputs$merge_table,
    sep = "\t", row.names = FALSE, quote = FALSE)
}

merge_df <- utils::read.delim(smk_outputs$merge_table)
message(sprintf("  %s transcripts -> %s bins",
  format(nrow(merge_df), big.mark = ","),
  format(length(unique(merge_df$tx_out)), big.mark = ",")))

message("Done.")
