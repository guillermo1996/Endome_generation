#!/usr/bin/env Rscript

################################################################################
## 06a - satuRn: Differential Transcript Usage on the scUTRquant ENDome
################################################################################
## Consumes the scUTRquant SingleCellExperiment (`*.txs.se.rds`) and runs DTU.
##
## scUTRquant SCE -> satuRn mapping:
##   assays(sce)$counts        -> count matrix (transcripts x samples)
##   rowData(sce)$transcript_id -> txInfo$isoform_id
##   rowData(sce)$gene_id       -> txInfo$gene_id
##   colData(sce)               -> experimental design (condition / donor / cell type)
##
## DTU "unit" is controlled by params$mode:
##   "pseudobulk" : sum counts per donor (x cell type); donor = biological replicate
##   "cells"      : each cell is a column, grouped by condition (x cell type)
##   "both"       : run both, write both result tables
##
## NOTE: requires a condition column in colData. Until the Hardy sample sheet /
## cell_annots carry case/control + cell-type labels this script will stop early
## with an informative error (by design — design is "not ready yet").
################################################################################

################################################################################
## Mock `snakemake` preamble
################################################################################

if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots = c(input = 'list', output = 'list', params = 'list',
              wildcards = 'list', log = "list", threads = 'numeric')
  )
  snakemake <- Snakemake(
    input  = list(sce = "debug_results/.../scUTRquant/data/sce/endome/hardy_snRNAseq.txs.se.rds"),
    output = list(pseudobulk = "results/saturn/endome.pseudobulk.DTU.tsv",
                  cells      = "results/saturn/endome.cells.DTU.tsv",
                  rds        = "results/saturn/endome.saturn.rds"),
    params = list(mode = "both",
                  condition_col = "group",
                  sample_col    = "sample_id.sq",
                  celltype_col  = NULL,         # or e.g. "cell_type"
                  contrasts     = c("case-control"),
                  min_count = 10, min_total_count = 30),
    log     = list("results/saturn/saturn.log"),
    threads = 8
  )
}

################################################################################
## Libraries and Parameters
################################################################################

suppressPackageStartupMessages({
  library(satuRn)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(BiocParallel)
  library(Matrix)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

log_con <- file(snakemake@log[[1]], open = "wt"); sink(log_con); sink(log_con, type = "message")

mode          <- snakemake@params$mode
condition_col <- snakemake@params$condition_col
sample_col    <- snakemake@params$sample_col
celltype_col  <- snakemake@params$celltype_col
contrasts     <- snakemake@params$contrasts
min_count       <- snakemake@params$min_count       %||% 10
min_total_count <- snakemake@params$min_total_count %||% 30

################################################################################
## Load scUTRquant SCE and validate design
################################################################################

sce <- readRDS(snakemake@input$sce)
cd  <- as.data.frame(colData(sce))

if (!condition_col %in% colnames(cd)) {
  stop(sprintf(
    "Condition column '%s' not found in colData(sce). Columns present: %s\n",
    condition_col, paste(colnames(cd), collapse = ", ")),
    "Add a case/control label to the scUTRquant sample sheet / cell_annots first.")
}

## transcript -> gene map from the SCE row annotations
txInfo_full <- data.frame(
  isoform_id = rowData(sce)$transcript_id,
  gene_id    = rowData(sce)$gene_id,
  stringsAsFactors = FALSE
)
rownames(txInfo_full) <- txInfo_full$isoform_id
txInfo_full <- txInfo_full[!is.na(txInfo_full$gene_id), ]

################################################################################
## Core satuRn runner (shared by both modes)
################################################################################

run_saturn <- function(counts, coldata) {
  ## 1. align txInfo to the matrix
  txInfo <- txInfo_full[rownames(counts), ]

  ## 2. feature filtering, then re-enforce >=2 isoforms / gene (satuRn requirement)
  design0 <- model.matrix(as.formula(paste0("~ ", condition_col)), data = coldata)
  keep    <- edgeR::filterByExpr(counts, design = design0,
                                 min.count = min_count, min.total.count = min_total_count)
  counts  <- counts[keep, , drop = FALSE]
  txInfo  <- txInfo[rownames(counts), ]
  multi   <- names(which(table(txInfo$gene_id) >= 2))
  txInfo  <- txInfo[txInfo$gene_id %in% multi, ]
  counts  <- counts[rownames(txInfo), , drop = FALSE]

  ## 3. assemble SummarizedExperiment + design formula
  se <- SummarizedExperiment(assays = list(counts = as.matrix(counts)),
                             colData = coldata, rowData = txInfo)
  metadata(se)$formula <- as.formula(paste0("~ 0 + ", condition_col))

  ## 4. fit + test
  se <- satuRn::fitDTU(object = se, formula = metadata(se)$formula,
                       parallel = TRUE, BPPARAM = MulticoreParam(snakemake@threads),
                       verbose = TRUE)

  design <- model.matrix(metadata(se)$formula, data = colData(se))
  colnames(design) <- levels(as.factor(colData(se)[[condition_col]]))
  L  <- limma::makeContrasts(contrasts = contrasts, levels = design)
  se <- satuRn::testDTU(object = se, contrasts = L, diagplot1 = FALSE, diagplot2 = FALSE)
  se
}

## flatten satuRn's per-contrast results (rowData(se)$fitDTUResult_<contrast>) to one tidy table
flatten_results <- function(se, mode_label) {
  res_cols <- grep("^fitDTUResult_", colnames(rowData(se)), value = TRUE)
  do.call(rbind, lapply(res_cols, function(rc) {
    df <- as.data.frame(rowData(se)[[rc]])
    df$isoform_id <- rownames(rowData(se))
    df$gene_id    <- rowData(se)$gene_id
    df$contrast   <- sub("^fitDTUResult_", "", rc)
    df$mode       <- mode_label
    df
  }))
}

################################################################################
## Build inputs per mode and run
################################################################################

results <- list(); se_objs <- list()

if (mode %in% c("pseudobulk", "both")) {
  key  <- interaction(cd[[sample_col]],
                      if (!is.null(celltype_col)) cd[[celltype_col]] else "all",
                      drop = TRUE)
  agg  <- t(rowsum(t(as.matrix(counts(sce))), group = key))          # tx x pseudobulk-sample
  pcol <- cd[match(levels(key), key), , drop = FALSE]
  rownames(pcol) <- levels(key)
  se_objs$pseudobulk <- run_saturn(agg, DataFrame(pcol))
  results$pseudobulk <- flatten_results(se_objs$pseudobulk, "pseudobulk")
  if (!is.null(snakemake@output$pseudobulk))
    write.table(results$pseudobulk, snakemake@output$pseudobulk, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (mode %in% c("cells", "both")) {
  se_objs$cells <- run_saturn(counts(sce), colData(sce))
  results$cells <- flatten_results(se_objs$cells, "cells")
  if (!is.null(snakemake@output$cells))
    write.table(results$cells, snakemake@output$cells, sep = "\t", quote = FALSE, row.names = FALSE)
}

saveRDS(se_objs, snakemake@output$rds)

cat("satuRn DTU complete. Modes run:", paste(names(se_objs), collapse = ", "), "\n")
sink(type = "message"); sink(); close(log_con)
