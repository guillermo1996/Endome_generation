## _________________________________________________
##
## Hardy snRNA-seq: cell annotation and sample tables for step 07
##
## Aim: Convert the Hardy snRNA-seq annotations and metadata into the two
##      tables that step 07 reads.
##
## Project: ENDome generation - UTR Usage
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-10-15
##
## Latest Version: v2.0 (2026-10-16)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
####
## Notes:
##    - Inputs:
##        - annotations: Melissa's nucleus annotations. The index is
##          "<barcode>-<annotation_label>"; the label is the sample
##          (donor_region), plus "B1"/"B2" for the samples sequenced twice.
##        - sample_metadata: the CRN SAMPLE.csv. Gives the S-number of each
##          library (in file_name), its replicate (rep1/rep2) and its region.
##        - sample_covariates: the covariates used in the manuscript.
##
## Changelog:
##    - v2.0 (2026-10-16): Standalone script (paths below, no Snakemake object).
##      Fewer checks.
##    - v1.0 (2026-10-15): Initial version.
##
## Contact: guillermorocamora@gmail.com
## _________________________________________________

############################################################################## #
# ---- 0. Setup ----

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages

shhh({
  library(conflicted)
  library(tidyverse)
})

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflict_prefer_all("tidyr", quiet = TRUE)

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
input_dir <- "data/Hardy_snRNAseq/annotations"

annotations_path <- file.path(input_dir, "annotation_all_levels_all_nuclei_updated.csv")
sample_metadata_path <- file.path(input_dir, "SAMPLE.csv")
sample_covariates_path <- file.path(input_dir, "Hardy_final_covariates_weighted_and_combined_batches.csv")

output_dir <- "data/Hardy_snRNAseq/prepared"
cell_annotations_out <- file.path(output_dir, "Hardy.cell_annotations.tsv.gz")
sample_table_out <- file.path(output_dir, "Hardy.sample_table.tsv")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

annotation_level <- "annotation_level_2"

############################################################################## #
# ---- 1. Libraries and their annotation labels ----
## One row per sequencing library (S-number in the FASTQ name)
libraries <- readr::read_csv(sample_metadata_path, col_types = readr::cols(.default = readr::col_character())) %>%
  dplyr::transmute(sample_id, replicate, region = brain_region, donor_id = subject_id,
                   s_number = stringr::str_match(file_name, "_(S\\d+)_")[, 2]) %>%
  dplyr::filter(!is.na(s_number)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(library_id = paste0(sample_id, "_", s_number))

## Annotation label of each library: the sample, plus B1/B2 for replicated samples
replicate_n <- as.integer(sub("^rep", "", libraries$replicate))
label_n <- dplyr::if_else(libraries$sample_id == "zupam_IPL", 3 - replicate_n, replicate_n)

libraries <- libraries %>%
  dplyr::mutate(
    sample_id_batch = dplyr::if_else(is.na(replicate), sample_id, paste0(sample_id, "B", replicate_n)),
    annotation_label = dplyr::if_else(is.na(replicate), sample_id, paste0(sample_id, "B", label_n))
  )

## borah_IPL and nobin_ACG: replicated, but annotated without suffix -> rep2 library
unsuffixed <- libraries %>%
  dplyr::filter(replicate == "rep2", sample_id %in% c("borah_IPL", "nobin_ACG")) %>%
  dplyr::mutate(annotation_label = sample_id)
label_map <- dplyr::bind_rows(libraries, unsuffixed) %>%
  dplyr::select(annotation_label, library_id, sample_id, sample_id_batch)
message(nrow(libraries), " libraries in ", dplyr::n_distinct(libraries$sample_id), " samples")

############################################################################## #
# ---- 2. Cell annotation table ----
annots <- readr::read_csv(annotations_path, col_types = readr::cols(.default = readr::col_character()))
annotation_cols <- grep("^annotation_level_", colnames(annots), value = TRUE)

cell_annotations <- annots %>%
  tidyr::separate_wider_regex(index, c(barcode = "[ACGTN]+", "-", annotation_label = ".+")) %>%
  dplyr::filter(annotation_level_0 != "drop", dplyr::if_all(dplyr::all_of(annotation_cols), ~ !is.na(.x)))

## Labels of samples that are not in SAMPLE.csv are dropped
unmatched <- setdiff(cell_annotations$annotation_label, label_map$annotation_label)
if (length(unmatched) > 0) message("Annotation labels without a library (removed): ", paste(unmatched, collapse = ", "))

cell_annotations <- cell_annotations %>%
  dplyr::inner_join(label_map, by = "annotation_label", relationship = "many-to-one") %>%
  dplyr::select(library_id, barcode, sample_id, sample_id_batch, annotation_label, celltype = .data[[annotation_level]]) %>% # dplyr::all_of(annotation_cols)) %>%
  dplyr::arrange(library_id, barcode)

message(nrow(cell_annotations), " nuclei in ", dplyr::n_distinct(cell_annotations$library_id), " libraries (",
        dplyr::n_distinct(cell_annotations$sample_id), " samples)")

############################################################################## #
# ---- 3. Sample table ----
## Covariates of the manuscript, plus region and donor from SAMPLE.csv
sample_table <- readr::read_csv(sample_covariates_path, col_types = readr::cols(sample = readr::col_character(), group = readr::col_character())) %>%
  dplyr::rename(sample_id = sample) %>%
  dplyr::left_join(libraries %>% dplyr::distinct(sample_id, region, donor_id), by = "sample_id", relationship = "one-to-one") %>%
  dplyr::select(sample_id, group, region, donor_id, dplyr::everything()) 

############################################################################## #
# ---- 4. Outputs ----
readr::write_tsv(cell_annotations, cell_annotations_out)
readr::write_tsv(sample_table, sample_table_out)
message("Written: ", cell_annotations_out, ", ", sample_table_out)
