## _________________________________________________
##
## txcutr vs gtxcutr (modified implementation)
##
## Aim: To test whether my modifications on the txcutr package change the output
## of Sonia's generated UTRome with the Ebbert data. The
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-04-08
##
## Copyright (c) Guillermo Rocamora Pérez, 2025
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-04-08)
##
## _________________________________________________
##
## - Notes:
##
## - Changelog:
##
## - Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(plyranges)
library(waldo)  # Testing this new library to replace "all.equal()" with a more interpretable output (https://waldo.r-lib.org/index.html)

# Import the GTFs
txcutr_gtf <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/data/Sonia_Ebbert_UTRome/ebbert_LR.all.txcutr_w500.gtf")
gtxcutr_gtf <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/gtxcutr/sonia_ebbert.gtxcutr.w500.3p.gtf.gz")

# Fix different seqlevelstyle
seqlevelsStyle(txcutr_gtf) <- "UCSC"

# Extract each feature type from the GTFs
txcutr_genes <- plyranges::filter(txcutr_gtf, type == "gene")
txcutr_transcripts <- plyranges::filter(txcutr_gtf, type == "transcript")
txcutr_exons <- plyranges::filter(txcutr_gtf, type == "exon")

gtxcutr_genes <- plyranges::filter(gtxcutr_gtf, type == "gene")
gtxcutr_transcripts <- plyranges::filter(gtxcutr_gtf, type == "transcript")
gtxcutr_exons <- plyranges::filter(gtxcutr_gtf, type == "exon")

# Ensure that number of elements is the same
waldo::compare(length(txcutr_genes), length(gtxcutr_genes))
waldo::compare(length(txcutr_transcripts), length(gtxcutr_transcripts))
waldo::compare(length(txcutr_exons), length(gtxcutr_exons))

# Ensure that the ranges of the elements are the same
waldo::compare(ranges(txcutr_genes), ranges(gtxcutr_genes))
waldo::compare(ranges(txcutr_transcripts), ranges(gtxcutr_transcripts))
waldo::compare(ranges(txcutr_exons), ranges(gtxcutr_exons))

# Ensure that the metadata columns are the same
waldo::compare(mcols(txcutr_genes), mcols(gtxcutr_genes))
waldo::compare(mcols(txcutr_transcripts), mcols(gtxcutr_transcripts))
waldo::compare(mcols(txcutr_exons), mcols(gtxcutr_exons))

# Ensure that both objects are exactly the same
waldo::compare(txcutr_gtf, gtxcutr_gtf)



# Same tests performed with stopifnot & all.equal to ensure that it is all working as expected
stopifnot(length(txcutr_gtf) == length(gtxcutr_gtf))
stopifnot(ranges(txcutr_gtf) == ranges(gtxcutr_gtf))
stopifnot(isTRUE(all.equal(mcols(txcutr_gtf), mcols(gtxcutr_gtf))))

# stopifnot(length(txcutr_genes) == length(gtxcutr_genes))
# stopifnot(length(txcutr_transcripts) == length(gtxcutr_transcripts))
# stopifnot(length(txcutr_exons) == length(gtxcutr_exons))
#
# stopifnot(ranges(txcutr_genes) == ranges(gtxcutr_genes))
# stopifnot(ranges(txcutr_transcripts) == ranges(gtxcutr_transcripts))
# stopifnot(ranges(txcutr_exons) == ranges(gtxcutr_exons))
#
# stopifnot(isTRUE(all.equal(mcols(txcutr_genes), mcols(gtxcutr_genes))))
# stopifnot(isTRUE(all.equal(mcols(txcutr_transcripts), mcols(gtxcutr_transcripts))))
# stopifnot(isTRUE(all.equal(mcols(txcutr_exons), mcols(gtxcutr_exons))))



# Test as TxDb
txcutr_txdb <- makeTxDbFromGFF("/home/grocamora/RytenLab-Research/38-Endome_generation/data/Sonia_Ebbert_UTRome/ebbert_LR.all.txcutr_w500.gtf")
gtxcutr_txdb <- makeTxDbFromGFF("/home/grocamora/RytenLab-Research/38-Endome_generation/results/gtxcutr/sonia_ebbert.gtxcutr.w500.3p.gtf.gz")

SseqlevelsStyle(txcutr_txdb) <- "UCSC"

waldo::compare(metadata(txcutr_txdb), metadata(gtxcutr_txdb))
waldo::compare(genes(txcutr_txdb), genes(gtxcutr_txdb))
waldo::compare(transcripts(txcutr_txdb), transcripts(gtxcutr_txdb))
waldo::compare(exons(txcutr_txdb), exons(gtxcutr_txdb))
waldo::compare(cds(txcutr_txdb), cds(gtxcutr_txdb))
