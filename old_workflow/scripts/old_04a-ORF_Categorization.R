#!/usr/bin/env Rscript

## _________________________________________________
##
## Transcript classification by ORF status
##
## Aim: Complete the GTF with ORF information
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-03-26
##
## Copyright (c) Guillermo Rocamora Pérez, 2025
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-03-26)
##
## _________________________________________________
##
## - Notes:
##
## - Changelog:
##
## - Please contact guillermorocamora@gmail.com for further assistance.
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
      threads = 'numeric'
    )
  )
  snakemake <- Snakemake(
    input = list(gtf = "~/RytenLab-Research/38-Endome_generation/results_k15/04-ORF_Filtration/gffread/sq3.annotated.gtf",
                 ref_annotation = "~/RytenLab-Research/Resources/GENCODE/gencode.v48.annotation.gtf",
                 sq3_class = "~/RytenLab-Research/38-Endome_generation/results/03-Sqanti3/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt"),
    output = list(gtf = "~/RytenLab-Research/38-Endome_generation/results/04-ORF_Filtration/ORF_Category/sq3.annotated_orf.gtf"),
    threads = 1
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

## Package options - these are some examples and do not modify the script usage.
options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)

#----------------------------------------------------------------------------- #
## 0.2 Script Paths ----

### Input Paths
input_gtf_path <- snakemake@input$gtf
input_reference_gtf_path <- snakemake@input$ref_annotation
input_sq3_class_path <- snakemake@input$sq3_class

### Output Paths
output_gtf_path <- snakemake@output$gtf

#### Create output directory
dir.create(dirname(output_gtf_path), showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.3 Script Parameters ----
#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----
#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----

############################################################################## #
# ---- 1. Load transcriptome and SQ3 Classification ----
gtf <- rtracklayer::import.gff(input_gtf_path)
annotation_gtf <- rtracklayer::readGFF(input_reference_gtf_path)
sq3_class <- readr::read_delim(input_sq3_class_path)


############################################################################## #
# ---- 2. Generate set of categories ----

#----------------------------------------------------------------------------- #
## 2.1 Assign categories to the transcripts ----

### GR: Categorization of the transcript is made based on the structural
### category assignation from Sqanti3 and the ORF presence. For more information
### regarding the possible categories from Sqanti3, please refer to their
### documentation:
### https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories
sq3_isoforms <- sq3_class %>%
  dplyr::mutate(sq3_type = dplyr::case_when(
    structural_category == "full-splice_match" & coding == "coding" ~ "known-ORF",
    structural_category != "full-splice_match" & coding == "coding" ~ "novel-ORF",
    coding != "coding" ~ "not-ORF",
    .default = "Unknown"
  )) %>%
  dplyr::mutate(hasCDS = !is.na(ORF_length),
                in_ref = !is.na(ref_length)) %>%
  dplyr::select(isoform, sq3_type, ORF_length, in_ref, hasCDS, associated_gene)

### GR: Load GENCODE annotations for gene and transcript types
annotation_genes <- annotation_gtf %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(gene_id, gene_type) %>%
  dplyr::mutate(gene_id = gsub("\\..*", "", gene_id))

annotation_transcripts <- annotation_gtf %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, transcript_type) %>%
  dplyr::mutate(transcript_id = gsub("\\..*", "", transcript_id))

isoform_types <- sq3_isoforms %>%
  dplyr::mutate(clean_iso = gsub("\\..*", "", isoform),
                clean_gene = gsub("\\..*", "", associated_gene)) %>%
  dplyr::left_join(annotation_genes, by = c("clean_gene" = "gene_id")) %>%
  dplyr::left_join(annotation_transcripts, by = c("clean_iso" = "transcript_id")) %>%
  dplyr::select(isoform, associated_gene, in_ref, hasCDS, ORF_length, sq3_type, ref_tx_type = transcript_type, ref_gene_type = gene_type)

#----------------------------------------------------------------------------- #
## 2.2 Map transcript ID to ORF Length ----

### GR: Set-up a dictionary to convert transcript id to the different other
### columns of interest
isoform_length <- isoform_types %>% dplyr::select(isoform, ORF_length) %>% tibble::deframe()

isoform_in_ref <- isoform_types %>% dplyr::select(isoform, in_ref) %>% tibble::deframe()
isoform_hasCDS <- isoform_types %>% dplyr::select(isoform, hasCDS) %>% tibble::deframe()
isoform_sq3_type <- isoform_types %>% dplyr::select(isoform, sq3_type) %>% tibble::deframe()
isoform_ref_tx_type <- isoform_types %>% dplyr::select(isoform, ref_tx_type) %>% tibble::deframe()
isoform_ref_gene_type <- isoform_types %>% dplyr::select(isoform, ref_gene_type) %>% tibble::deframe()


############################################################################## #
# ---- 3. Add information to the GTF ----

#----------------------------------------------------------------------------- #
## 3.1 Add ORF Category ----
mcols(gtf)[gtf$type == "transcript", "in_ref"] <- isoform_in_ref[mcols(gtf)[gtf$type == "transcript", "transcript_id"]]
mcols(gtf)[gtf$type == "transcript", "hasCDS"] <- isoform_hasCDS[mcols(gtf)[gtf$type == "transcript", "transcript_id"]]

#----------------------------------------------------------------------------- #
## 3.2 Add ORF Length ----
mcols(gtf)[gtf$type == "transcript", "ORF_length"] <- isoform_length[mcols(gtf)[gtf$type == "transcript", "transcript_id"]]

#----------------------------------------------------------------------------- #
## 3.3 Add Transcript Categories Category ----
mcols(gtf)[gtf$type == "transcript", "sq3_type"] <- isoform_sq3_type[mcols(gtf)[gtf$type == "transcript", "transcript_id"]]
mcols(gtf)[gtf$type == "transcript", "ref_tx_type"] <- isoform_ref_tx_type[mcols(gtf)[gtf$type == "transcript", "transcript_id"]]
mcols(gtf)[gtf$type == "transcript", "ref_gene_type"] <- isoform_ref_gene_type[mcols(gtf)[gtf$type == "transcript", "transcript_id"]]

############################################################################## #
# ---- 4. Output the GTF ----
gtf %>% rtracklayer::export(output_gtf_path)
