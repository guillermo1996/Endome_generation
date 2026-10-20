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
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomeInfoDb)
  library(GenomeInfoDbData)
  library(BiocParallel)
  library(magrittr)
  library(txdbmaker)
  library(tidyverse)
})

### Conflicts - declare function preferences
shhh <- suppressPackageStartupMessages # Shortcut to hide package start-up messages
shhh(library(conflicted))
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflict_prefer_all("tidyr", quiet = TRUE)

## convert arguments
maxTxLength <- snakemake@wildcards$width
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
## New Truncation Function
################################################################################
truncateTxome_mod <- function(txdb,
                              txEnd = "3prime",
                              overlapFile = NULL){
  ############################################################################
  # Ensure correct values of `txEnd`
  valid_3prime <- c("3", "3'", "3p", "3prime", "3_prime")
  valid_5prime <- c("5", "5'", "5p", "5prime", "5_prime")
  
  if (txEnd %in% valid_3prime) txEnd <- "3prime"
  if (txEnd %in% valid_5prime) txEnd <- "5prime"
  
  if (!txEnd %in% c("3prime", "5prime")) stop("txEnd parameter not valid - only '3prime' or '5prime' parameters are accepted.")
  
  ############################################################################
  # Split exons by transcripts and create a mapping dictionary from
  # transcript_id to gene_id
  grlExons <- exonsBy(txdb, use.names = TRUE)
  dfTxGene <- suppressMessages(AnnotationDbi::select(txdb, keys = names(grlExons), keytype = "TXNAME", columns = "GENEID"))
  mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)
  
  ############################################################################
  # Hijack truncation logic
  grl <- switch(txEnd,
    "3prime" = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE),
    "5prime" = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE),
    stop("Unknown txEnd: ", txEnd)
  )
  
  clipped <- slot(txendcutr:::.mutateEach(grl, transcript_id = names(grl)), "unlistData")
  clipped$transcript_id <- factor(clipped$transcript_id, levels = names(grl))
  
  ############################################################################
  # Removal of duplicated transcripts
  grlC_clipped <- S4Vectors::split(clipped, mcols(clipped)["transcript_id"])
  grlC_clipped_names <- names(grlC_clipped)
  
  ## Generate the overlap: look for exact matches between transcripts
  message("Checking for duplicate transcripts...")
  overlaps <- findOverlaps(grlC_clipped, type="equal",
                           ignore.strand=FALSE,
                           drop.self=TRUE, drop.redundant=TRUE)
  
  ##  Ensure that overlaps are from the same gene
  matched_overlaps <- tibble::as_tibble(overlaps) %>% 
    dplyr::mutate(queryTx = grlC_clipped_names[queryHits],
                  subjectTx = grlC_clipped_names[subjectHits]) %>% 
    dplyr::mutate(queryGene = mapTxToGene[queryTx],
                  subjectGene = mapTxToGene[subjectTx]) %>% 
    dplyr::filter(queryGene == subjectGene)
  
  ## Export overlap data.frame
  if (!is.null(overlapFile) && overlapFile != ""){
    output_dir <- dirname(overlapFile)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive=TRUE)
    }
    
    ## Store output in disk
    write.table(matched_overlaps, overlapFile, sep="\t", row.names=FALSE, quote=FALSE)
    message(sprintf("Post-truncation transcript overlaps exported to: %s", overlapFile))
  }
  
  ## Remove overlaps
  if (nrow(matched_overlaps) > 0) {
    tx_to_remove <- unique(matched_overlaps$queryHits)
    grlC_clipped <- grlC_clipped[-tx_to_remove]
    message(sprintf("Removed %d duplicates.", length(tx_to_remove)))
  } else {
    message("No duplicated transcripts found.")
  }
  
  ############################################################################
  # Create the final exon ranges
  message("Creating exon ranges...")
  
  ## flatten with tx_id in metadata
  grExons <- slot(grlC_clipped, "unlistData")
  mcols(grExons)["type"] <- "exon"
  
  ## add gene id
  mcols(grExons)["gene_id"] <- mapTxToGene[as.character(mcols(grExons)$transcript_id)]
  
  ## reindex exon info
  grExons <- sort(grExons)
  mcols(grExons)["exon_id"] <- seq_along(grExons)
  mcols(grExons)["exon_name"] <- NULL
  message("Done.")
  
  ############################################################################
  # Create the final transcript ranges
  message("Creating tx ranges...")
  grTxs <- unlist(range(grlC_clipped))
  mcols(grTxs)["transcript_id"] <- factor(names(grlC_clipped), levels = levels(grExons$transcript_id))
  mcols(grTxs)["type"] <- "transcript"
  
  ## add gene id
  mcols(grTxs)["gene_id"] <- mapTxToGene[as.character(grTxs$transcript_id)]
  
  grTxs <- grTxs[order(grTxs$transcript_id)]
  message("Done.")
  
  ############################################################################
  # Create the final gene ranges
  message("Creating gene ranges...")
  grGenes <- unlist(range(S4Vectors::split(grTxs, mcols(grTxs)$gene_id)))
  mcols(grGenes)["gene_id"] <- names(grGenes)
  mcols(grGenes)["type"] <- "gene"
  message("Done.")
  
  ############################################################################
  # Generate the final TxDb object
  dfMetadata <- data.frame(
    name=c("Truncated by", "Maximum Transcript Length", "Truncation End"),
    value=c("txendcutr", "UTR", txEnd)
  )
  
  txendcutr:::.suppressTxDbGenomeWarning(
    txdbmaker::makeTxDbFromGRanges(c(grGenes, grTxs, grExons),
                        taxonomyId = taxonomyId(txdb),
                        metadata = dfMetadata
    )
  )
}

generateMergeTable_mod <- function(txdb, type = NULL){
  ############################################################################
  # Extract txEnd from txdb metadata
  if (!"Truncation End" %in% metadata(txdb)$name) {
    warning("'Truncation End' parameter not found in the TxDb metadata. Defaults to '3prime'. Was this object created with txendcutr::truncateTxome?")
    txEnd <- "3prime"
  } else {
    txEnd <- metadata(txdb)[metadata(txdb)$name == "Truncation End", "value"]
  }
  if (!txEnd %in% c("3prime", "5prime")) stop("'Truncation End' parameter is not valid - only '3prime' or '5prime' are accepted.")
  
  ############################################################################
  # Merge pipeline
  grTxs <- transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"))

  ## Invert strand if 5' truncation
  if (txEnd == "5prime") grTxs <- invertStrand(grTxs)
  
  overlaps <- tibble::as_tibble(findOverlaps(grTxs, type = type, ignore.strand = FALSE, drop.self = TRUE))
  overlaps["tx_in"] <- unlist(grTxs$tx_name[overlaps$queryHits])
  overlaps["tx_out"] <- unlist(grTxs$tx_name[overlaps$subjectHits])
  overlaps["gene_in"] <- unlist(grTxs$gene_id[overlaps$queryHits])
  overlaps["gene_out"] <- unlist(grTxs$gene_id[overlaps$subjectHits])
  
  ## filter unmatched genes
  overlaps <- overlaps[overlaps$gene_in == overlaps$gene_out, ]
  
  ## keep downstream
  overlaps["strand"] <- as.character(strand(grTxs))[overlaps$queryHits]
  overlaps["end_in"] <- ifelse(overlaps$strand == "+", end(grTxs[overlaps$queryHits]), -start(grTxs[overlaps$queryHits]))
  overlaps["end_out"] <- ifelse(overlaps$strand == "+", end(grTxs[overlaps$subjectHits]), -start(grTxs[overlaps$subjectHits]))
  overlaps <- overlaps[overlaps$end_in <= overlaps$end_out, ]
  
  idxEqual <- which(overlaps$end_in == overlaps$end_out)
  if (length(idxEqual) > 0) {
    isLaterTx <- overlaps[idxEqual, "tx_in"] < overlaps[idxEqual, "tx_out"]
    overlaps <- overlaps[-idxEqual[isLaterTx], ]
  }
  
  ## pick unique out
  dfMerge <- overlaps %>% 
    dplyr::group_by(queryHits) %>% 
    dplyr::arrange(-end_out, tx_out, .by_group = T) %>% 
    dplyr::slice_head(n = 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(tx_in, tx_out)
  dfMerge <- txendcutr:::.propagateMap(dfMerge)
  
  ## include self-maps
  unmergedTxs <- grTxs$tx_name[!(grTxs$tx_name %in% dfMerge$tx_in)]
  dfMerge <- rbind(dfMerge, data.frame(tx_in = unmergedTxs, tx_out = unmergedTxs))
  
  ## append gene
  txToGeneMap <- setNames(object = unlist(grTxs$gene_id), nm = grTxs$tx_name)
  dfMerge["gene_out"] <- txToGeneMap[as.character(dfMerge$tx_out)]
  
  ## reorder and drop rownames
  dfMerge <- dfMerge[with(dfMerge, order(tx_in)), ]
  rownames(dfMerge) <- NULL
  
  dfMerge
}

exportMergeTable_mod <- function(txdb, file, type = NULL){
  if(is.null(type) || type == "default"){
    df <- generateMergeTable_mod(txdb)
  }else if(type == "equal"){
    df <- generateMergeTable_mod(txdb, type = "equal")
  }else if(type == "bin"){
    df <- txendcutr::generateMergeTable(txdb)
  }else if(type == "empty"){
    grTxs <- transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"))
    unmergedTxs <- grTxs$tx_name
    dfMerge <- data.frame(tx_in = unmergedTxs, tx_out = unmergedTxs)

    txToGeneMap <- setNames(object = unlist(grTxs$gene_id), nm = grTxs$tx_name)
    dfMerge["gene_out"] <- txToGeneMap[as.character(dfMerge$tx_out)]
    dfMerge <- dfMerge[with(dfMerge, order(tx_in)), ]
    rownames(dfMerge) <- NULL
    df <- dfMerge
  }else{
    stop()
  }

  if (!inherits(file, "connection")) {
    if (grepl(".gz$", file)) {
      file <- gzfile(file, "wb")
      on.exit(close(file))
    }
  }
  write.table(df, file, sep="\t", row.names=FALSE, quote=FALSE)
  invisible(txdb)
}

################################################################################
## Load Data, Truncate, and Export
################################################################################
txdb <- txdbmaker::makeTxDbFromGFF(file=snakemake@input$gtf, organism=organism(bsg))

txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
seqlevelsStyle(txdb) <- "UCSC"

# overlap_path is optional
# txEnd is optional. Defaults to 3' truncation

txdb_result <- truncateTxome_mod(txdb, overlapFile = overlap_path, txEnd = txEnd)

print("Export GTF")
exportGTF(txdb_result, snakemake@output$gtf)

print("Export FASTA")
exportFASTA(txdb_result, bsg, snakemake@output$fa)

print("Export MergeTable")

if(snakemake@wildcards$width == "UTR"){
  type = "default"
}else if(snakemake@wildcards$width == "UTR_empty"){
  type = "empty"
}else if(snakemake@wildcards$width == "UTR_equal"){
  type = "equal"
}else if(snakemake@wildcards$width == "UTR_bin"){
  type = "bin"
}else{
  stop("Not recognized ", snakemake@wildcards$width)
}
exportMergeTable_mod(txdb_result, snakemake@output$merge_table, type = type)
