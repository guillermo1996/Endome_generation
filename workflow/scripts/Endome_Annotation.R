## _________________________________________________
##
## Script title
##
## Aim:
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-03-27
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-03-27)
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
    input = list(gtf = "/home/grocamora/RytenLab-Research/38-Endome_generation/results/gtxcutr/sq3.annotated.gtxcutr.w500.gtf.gz"),
    output = list(gtf = "results/Endome_Annotation/sq3.annotated.w500.annotated.gtf"),
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
main_path <- "/home/grocamora/RytenLab-Research/38-Endome_generation" # here::here()

### Input Paths
input_gtf_path <- snakemake@input$gtf
input_sq3_class_path <- snakemake@input$sq3_class

### Output Paths
output_gtf_path <- snakemake@output$gtf

#### Create output directory
# dir.create(dirname(output_gtf_path), showWarnings = F, recursive = T)

#----------------------------------------------------------------------------- #
## 0.3 Script Parameters ----

#----------------------------------------------------------------------------- #
## 0.4 User Input arguments ----

#----------------------------------------------------------------------------- #
## 0.5 Helper Functions ----

############################################################################## #
# ---- 1. Load transcriptome and SQ3 Classification ----
gtf <- rtracklayer::import.gff(input_gtf_path)
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
  dplyr::mutate(ORF_cat = dplyr::case_when(
    structural_category == "full-splice_match" & coding == "coding" ~ "known-ORF",
    structural_category != "full-splice_match" & coding == "coding" ~ "novel-ORF",
    coding != "coding" ~ "not-ORF",
    .default = "Unknown"
  )) %>%
  dplyr::select(isoform, ORF_cat, ORF_length)

#----------------------------------------------------------------------------- #
## 2.2 Define the transcript lists ----
tx_known <- sq3_isoforms$isoform[sq3_isoforms$ORF_cat == "known-ORF"]
tx_novel <- sq3_isoforms$isoform[sq3_isoforms$ORF_cat == "novel-ORF"]
tx_not <- sq3_isoforms$isoform[sq3_isoforms$ORF_cat == "not-ORF"]

#----------------------------------------------------------------------------- #
## 2.3 Map transcript ID to ORF Length ----

### GR: Set-up a dictionary to convert transcript id to ORF length. This metric
### might be useful in future modules.
isoform_length <- sq3_isoforms %>% dplyr::select(isoform, ORF_length) %>% tibble::deframe()


############################################################################## #
# ---- 3. Add information to the GTF ----

#----------------------------------------------------------------------------- #
## 3.1 Add ORF Category ----
mcols(gtf)$ORF_cat <- NA

mcols(gtf)[gtf$type == "transcript" & gtf$transcript_id %in% tx_known, "ORF_cat"] <- "known-ORF"
mcols(gtf)[gtf$type == "transcript" & gtf$transcript_id %in% tx_novel, "ORF_cat"] <- "novel-ORF"
mcols(gtf)[gtf$type == "transcript" & gtf$transcript_id %in% tx_not, "ORF_cat"] <- "not-ORF"


#----------------------------------------------------------------------------- #
## 3.2 Add ORF Length ----
mcols(gtf)$ORF_length <- NA
mcols(gtf)[gtf$ORF_cat %in% c("known-ORF", "novel-ORF"), "ORF_length"] <- isoform_length[mcols(gtf)[gtf$ORF_cat %in% c("known-ORF", "novel-ORF"), "transcript_id"]]


############################################################################## #
# ---- 4. Output the GTF ----
gtf %>% rtracklayer::export(output_gtf_path)
