## _________________________________________________
##
## Build the ENDome DuckDB
##
## Aim: Build one DuckDB per ENDome GTF from the outputs of the truncation
## pipeline and the upstream ORF categorization
##
## Project: ENDome generation - Bin Information Content
##
## Author: Guillermo Rocamora Pérez (guillermorocamora@gmail.com)
##
## Date Created: 2026-08-19
##
## Latest Version: v1.1 (2026-09-30)
##
## Copyright (c) Guillermo Rocamora Pérez, 2026
##
## _________________________________________________
##
## Notes:
##    - One DuckDB per (prefix, orf_filter, width, txEnd), named to mirror the
##      ENDome GTF it describes.
##    - The `transcripts` table is built from the PRE-truncation ORF-filtered
##      GTF. `bin_id` / `superseded_by` are the only truncation-dependent columns
##    - Surrogate keys (aa_id, cds_id, utr5_id, utr3_id) are content hashes.
##
## Changelog:
##    - v1.1 (2026-09-30): transcripts table follows 04a v1.2: in_ref replaced
##      by ref_isoform; added gene_source, merge_gene_id and ref_transcript_id.
##    - v1.0 (2026-08-19): Initial version
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
  ## dataset.group.merge_method folder created by the `test_data` rule (its
  ## name is the prefix wildcard), and the orf_filter preset and truncation
  ## (width, txEnd) to test. The output uses the "interactive" prefix so it
  ## never overwrites the linked results.
  test_dir <- "data/test_data/Ebbert.control.iso_ref"
  test_orf_filter <- "ref_pc"
  test_width <- "500"
  test_txEnd <- "3p"
  test_txendcutr <- file.path(test_dir, paste0("test.", test_orf_filter, ".txendcutr.w", test_width, ".", test_txEnd))
  snakemake <- Snakemake(
    input = list(
      gtf_filter = file.path(test_dir, paste0("test.", test_orf_filter, ".orf_filter.gtf")),
      isoform_summary = file.path(test_dir, "test.isoform_summary.tsv"),
      protein_fa = file.path(test_dir, "protein.fa"),
      utr5_fa = file.path(test_dir, "utr5.fa"),
      utr3_fa = file.path(test_dir, "utr3.fa"),
      merge_table = paste0(test_txendcutr, ".merge.tsv"),
      overlap_table = paste0(test_txendcutr, ".overlaps.tsv"),
      ref_genome = "/home/MinaRyten/Guillermo/Resources/Genome/GRCh38.primary_assembly.genome.fa"
    ),
    output = list(
      duckdb = file.path(test_dir, paste0("interactive.", test_orf_filter, ".w", test_width, ".", test_txEnd, ".duckdb"))
    ),
    params = list(),
    wildcards = list(
      prefix = basename(test_dir),
      orf_filter = test_orf_filter,
      width = test_width,
      txEnd = test_txEnd
    ),
    threads = 4,
    scriptdir = "workflow/scripts"
  )
}

#----------------------------------------------------------------------------- #
## 0.1 Required Libraries ----
shhh <- suppressPackageStartupMessages

shhh({
  library(GenomicFeatures)
  library(txdbmaker)
  library(Biostrings)
  library(digest)
  library(DBI)
  library(duckdb)
  library(igraph)
  library(conflicted)
  library(plyranges)
  library(tidyverse)
})

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflict_prefer_all("tidyr", quiet = TRUE)

#----------------------------------------------------------------------------- #
## 0.2 Script Parameters ----
smk_inputs <- snakemake@input
smk_outputs <- snakemake@output
params <- snakemake@params
wc <- snakemake@wildcards

run_id <- paste(wc$prefix, wc$orf_filter, paste0("w", wc$width), wc$txEnd, sep = ".")

### Load the reference genome
fa <- Rsamtools::FaFile(smk_inputs$ref_genome)
if (!file.exists(Rsamtools::index(fa))) Rsamtools::indexFa(fa)
open(fa)

#----------------------------------------------------------------------------- #
## 0.3 Helper Functions ----
lib_dir <- file.path(snakemake@scriptdir, "../lib")
source(file.path(lib_dir, "Phf-DuckDB.R"))
source(file.path(lib_dir, "Phf-Bin_Information_Score.R"))

#' Start/stop codon sequence and canonicality per transcript
#'
#' @param grl_cds Named GRangesList of CDS exons per transcript (e.g. `cdsBy()`).
#' @param genome BSgenome/FaFile to pull sequence from.
#'
#' @returns Tibble with one row per transcript with `transcript_id`,
#'   `start_codon_sequence`, `start_codon_canonical`, `stop_codon_sequence`,
#'   `stop_codon_canonical` (ATG for start; TAA/TAG/TGA for stop).
extractStartStopCodonTbl <- function(grl_cds, genome = fa){
  stopifnot(!is.null(names(grl_cds)))

  cds_widths <- sum(width(grl_cds))
  codon_tx <- list(
    start = GRanges(seqnames = names(grl_cds), ranges = IRanges(1, 3)),
    stop  = GRanges(seqnames = names(grl_cds), ranges = IRanges(cds_widths - 2, cds_widths))
  )
  canonical_codons <- list(start = "ATG", stop = c("TAA", "TAG", "TGA"))

  codon_tbl <- function(type){
    codons_gr <- unlist(pmapFromTranscripts(codon_tx[[type]], grl_cds))
    codons_gr$transcript_id <- names(codons_gr)
    codons_gr <- codons_gr %>% dplyr::filter(hit)

    seqs <- extractTranscriptSeqs(genome, GenomicRanges::split(codons_gr, codons_gr$transcript_id))

    tibble::tibble(
      transcript_id = names(seqs),
      sequence = as.character(seqs)
    ) %>%
      dplyr::mutate(canonical = sequence %in% canonical_codons[[type]]) %>%
      dplyr::rename_with(~ paste0(type, "_codon_", .x), c(sequence, canonical))
  }

  dplyr::full_join(codon_tbl("start"), codon_tbl("stop"), by = "transcript_id")
}

############################################################################## #
# ---- 1. Load the pipeline outputs ----

#----------------------------------------------------------------------------- #
## 1.1 Tables ----
message("Loading isoform summary / merge table / overlaps table...")

isoform_summary <- readr::read_tsv(smk_inputs$isoform_summary)
merge_table <- readr::read_tsv(smk_inputs$merge_table, col_types = "ccc")
overlaps_table <- readr::read_tsv(smk_inputs$overlap_table)

#----------------------------------------------------------------------------- #
## 1.2 Transcriptome ----
message("Building TxDb from the full-length ORF-filtered GTF...")

txdb <- txdbmaker::makeTxDbFromGFF(smk_inputs$gtf_filter)
txdb <- GenomeInfoDb::keepStandardChromosomes(txdb, pruning.mode = "coarse")
GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"

############################################################################## #
# ---- 2. Build the DuckDB tables ----

#----------------------------------------------------------------------------- #
## 2.1 Proteins ----
message("Extracting the protein sequences from ORFannotate...")
protein_fa_df <- StringSetToTibble(Biostrings::readAAStringSet(smk_inputs$protein_fa))

### Remove the trailing stop character
protein_map <- protein_fa_df %>%
  dplyr::filter(!is.na(sequence), nzchar(sequence)) %>%
  dplyr::transmute(
    transcript_id,
    aa_seq = sub("\\*+$", "", sequence),
    aa_len = nchar(aa_seq)
  ) %>%
  dplyr::filter(aa_len > 0) %>%
  dplyr::mutate(aa_id = assignId(., "aa_seq", prefix = "aa"))

protein_tbl <- protein_map %>%
  dplyr::distinct(aa_id, aa_len, aa_seq) %>%
  dplyr::mutate(run_id = run_id, .before = 0)

transcript_aa_map <- protein_map %>% dplyr::select(transcript_id, aa_id)

#----------------------------------------------------------------------------- #
## 2.2 CDS / 5'UTR / 3'UTR ----
message("Extracting the CDS / UTR structure and sequence...")

grl_cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
grl_utr5 <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
grl_utr3 <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)

utr5_fa_df <- StringSetToTibble(Biostrings::readDNAStringSet(smk_inputs$utr5_fa))
utr3_fa_df <- StringSetToTibble(Biostrings::readDNAStringSet(smk_inputs$utr3_fa))

### ORFannotate emits the UTR FASTAs but not a CDS FASTA, so the CDS sequence is
### spliced out of the reference genome here.
start_stop_codons_df <- extractStartStopCodonTbl(grl_cds, genome = fa)
seq_cds_df <- StringSetToTibble(GenomicFeatures::extractTranscriptSeqs(fa, grl_cds))

close(fa) # Last read from the reference genome

### 2.2.1 CDS ----
cds_map <- tibble::as_tibble(grl_cds) %>% 
  dplyr::transmute(transcript_id = group_name, seqnames, start, end, strand) %>% 
  dplyr::group_by(transcript_id, seqnames, strand) %>% 
  dplyr::summarise(n_exons = dplyr::n(), n_junctions = n_exons - 1, start = min(start), end = max(end), .groups = "drop") %>% 
  dplyr::mutate(locus = sprintf("%s:%i-%i:%s", seqnames, start, end, strand)) %>% 
  dplyr::left_join(seq_cds_df %>% dplyr::select(transcript_id, len = width, seq = sequence), by = "transcript_id") %>% 
  dplyr::left_join(start_stop_codons_df, by = "transcript_id") %>%
  dplyr::left_join(transcript_aa_map, by = "transcript_id") %>% 
  dplyr::select(-seqnames, -strand, -start, -end) %>% 
  dplyr::mutate(cds_id = assignId(., c("locus", "seq", "n_exons"), prefix = "cds"), .before = 0) %>% 
  dplyr::relocate(transcript_id, cds_id, aa_id, locus, n_exons, n_junctions) %>% 
  dplyr::rename(cds_locus = locus, cds_len = len, cds_seq = seq)

cds_tbl <- cds_map %>% 
  dplyr::select(-transcript_id) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(run_id = run_id, .before = 0)

transcript_cds_map <- cds_map %>% dplyr::select(transcript_id, cds_id)

### 2.2.2 5' UTR ----
utr5_map <- tibble::as_tibble(grl_utr5) %>%
  dplyr::transmute(transcript_id = group_name, seqnames, start, end, strand) %>%
  dplyr::group_by(transcript_id, seqnames, strand) %>% 
  dplyr::summarise(n_exons = dplyr::n(), n_junctions = n_exons - 1, start = min(start), end = max(end), .groups = "drop") %>% 
  dplyr::mutate(locus = sprintf("%s:%i-%i:%s", seqnames, start, end, strand)) %>% 
  dplyr::left_join(utr5_fa_df %>% dplyr::select(transcript_id, len = width, seq = sequence), by = "transcript_id") %>% 
  dplyr::select(-seqnames, -strand, -start, -end) %>% 
  dplyr::mutate(utr5_id = assignId(., c("locus", "seq", "n_exons"), prefix = "utr5"), .before = 0) %>%
  dplyr::relocate(transcript_id, utr5_id, locus, n_exons, n_junctions) %>% 
  dplyr::rename(utr5_locus = locus, utr5_len = len, utr5_seq = seq)

utr5_tbl <- utr5_map %>% 
  dplyr::select(-transcript_id) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(run_id = run_id, .before = 0)

transcript_utr5_map <- utr5_map %>% dplyr::select(transcript_id, utr5_id)

### 2.2.3 3' UTR ----
utr3_map <- tibble::as_tibble(grl_utr3) %>%
  dplyr::transmute(transcript_id = group_name, seqnames, start, end, strand) %>%
  dplyr::group_by(transcript_id, seqnames, strand) %>% 
  dplyr::summarise(n_exons = dplyr::n(), n_junctions = n_exons - 1, start = min(start), end = max(end), .groups = "drop") %>% 
  dplyr::mutate(locus = sprintf("%s:%i-%i:%s", seqnames, start, end, strand)) %>% 
  dplyr::left_join(utr3_fa_df %>% dplyr::select(transcript_id, len = width, seq = sequence), by = "transcript_id") %>% 
  dplyr::select(-seqnames, -strand, -start, -end) %>% 
  dplyr::mutate(utr3_id = assignId(., c("locus", "seq", "n_exons"), prefix = "utr3"), .before = 0) %>%
  dplyr::relocate(transcript_id, utr3_id, locus, n_exons, n_junctions) %>% 
  dplyr::rename(utr3_locus = locus, utr3_len = len, utr3_seq = seq)

utr3_tbl <- utr3_map %>% 
  dplyr::select(-transcript_id) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(run_id = run_id, .before = 0)

transcript_utr3_map <- utr3_map %>% dplyr::select(transcript_id, utr3_id)

#----------------------------------------------------------------------------- #
## 2.3 Bins and their members ----
direct_members <- merge_table %>% 
  dplyr::transmute(bin_id = tx_out, transcript_id = tx_in, is_superseded = FALSE, superseded_by = NA_character_)

overlap_graph <- igraph::graph_from_data_frame(overlaps_table %>% dplyr::select(queryTx, subjectTx))
superseded_df <- resolveTerminalNodes(overlap_graph, report_self = FALSE)

superseded_members <- superseded_df %>% 
  dplyr::rename(transcript_id = node, superseded_by = terminal_node) %>% 
  dplyr::filter(!transcript_id %in% direct_members$transcript_id) %>% 
  dplyr::inner_join(direct_members %>% dplyr::select(bin_id, superseded_by = transcript_id), by = "superseded_by") %>% 
  dplyr::mutate(is_superseded = TRUE)

bin_members <- dplyr::bind_rows(direct_members, superseded_members)

bins_tbl <- bin_members %>% 
  dplyr::left_join(merge_table %>% dplyr::distinct(bin_id = tx_out, gene_id = gene_out), by = "bin_id") %>% 
  dplyr::group_by(bin_id, gene_id) %>% 
  dplyr::summarise(n_members = dplyr::n(), n_superseded = sum(is_superseded), .groups = "drop") %>%
  dplyr::mutate(run_id = run_id, .before = 0)

#----------------------------------------------------------------------------- #
## 2.4 Transcripts ----
message("Aggregating the transcript information...")

transcript_tbl <- transcripts(txdb) %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(isoform_summary, by = c("tx_name" = "transcript_id")) %>% 
    dplyr::select(
    transcript_id = tx_name, gene_id, gene_name, gene_source, merge_gene_id,
    seqnames, start, end, strand,
    structural_category, subcategory, ref_isoform, ref_transcript_id,
    coding_prob, ref_transcript_type, ref_gene_type,
    has_orf, total_junctions, ref_exons, NMD_sensitive
  ) %>%
  dplyr::mutate(run_id = run_id, .before = 0) %>% 
  dplyr::left_join(transcript_cds_map, by = "transcript_id") %>% 
  dplyr::left_join(transcript_utr5_map, by = "transcript_id") %>% 
  dplyr::left_join(transcript_utr3_map, by = "transcript_id") %>% 
  dplyr::left_join(transcript_aa_map, by = "transcript_id") %>% 
  dplyr::left_join(bin_members %>% dplyr::select(transcript_id, bin_id, superseded_by), by = "transcript_id") %>% 
  dplyr::relocate(run_id, transcript_id, bin_id, gene_id, aa_id, cds_id, utr5_id, utr3_id, superseded_by)

############################################################################## #
# ---- 3. Sanity checks ----

### Every bin member must be a transcript of this transcriptome, or the foreign
### keys below fail with a much less readable error.
orphan_members <- setdiff(bin_members$transcript_id, transcript_tbl$transcript_id)
if (length(orphan_members) > 0) {
  stop(sprintf("%i bin member(s) absent from the ORF-filtered GTF: %s",
    length(orphan_members), paste(utils::head(orphan_members, 10), collapse = ", ")), call. = FALSE)
}

### One row per transcript
stopifnot(!any(duplicated(transcript_tbl$transcript_id)))
stopifnot(!any(duplicated(bin_members[c("bin_id", "transcript_id")])))

############################################################################## #
# ---- 4. Write the DuckDB ----
message("Writing ", smk_outputs$duckdb, " ...")
if(file_test("-f", smk_outputs$duckdb)) file.remove(smk_outputs$duckdb)
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = smk_outputs$duckdb)

duckdb_schema <- list(
  proteins = list(df = protein_tbl, pk = "aa_id"),
  bins = list(df = bins_tbl, pk = "bin_id"),
  utr5 = list(df = utr5_tbl, pk = "utr5_id"),
  utr3 =  list(df = utr3_tbl, pk = "utr3_id"),
  cds = list(df = cds_tbl, pk = "cds_id", fk = list(aa_id = "proteins.aa_id")),
  transcripts = list(
    df = transcript_tbl, pk = "transcript_id",
    fk = list(
      bin_id = "bins.bin_id",
      utr5_id = "utr5.utr5_id",
      utr3_id = "utr3.utr3_id",
      cds_id = "cds.cds_id",
      superseded_by = "transcripts.transcript_id"
    )
  )
)

for(tbl_name in names(duckdb_schema)){
  n <- writeDuckTable(con, tbl_name, duckdb_schema[[tbl_name]])
  message(sprintf("  %-12s %8s rows", tbl_name, format(n, big.mark = ",")))
}

DBI::dbDisconnect(con, shutdown = TRUE)
message("Done.")

############################################################################## #
# ---- 5. Interactive only steps ----
if(interactive()){
  shhh(library(dm))

  dm_draw_svg = function(dm, ...){
    if (!requireNamespace("DiagrammeRsvg", quietly = TRUE)) {
      stop(
        "Package \"DiagrammeRsvg\" must be installed to use this function.",
        call. = FALSE
      )
    }
    
    dm::dm_draw(dm = dm, ...) %>%
      DiagrammeRsvg::export_svg() %>%
      htmltools::HTML() %>%
      htmltools::html_print()
  }

  utrome <- dm(bins_tbl, utr5_tbl, utr3_tbl, cds_tbl, protein_tbl, transcript_tbl)
  dm_examine_constraints(utrome) 

  utrome <- utrome %>% 
    dm_add_pk(bins_tbl, bin_id) %>% 
    dm_add_pk(utr5_tbl, utr5_id) %>% 
    dm_add_pk(utr3_tbl, utr3_id) %>% 
    dm_add_pk(cds_tbl, cds_id) %>% 
    dm_add_pk(protein_tbl, aa_id) %>% 
    dm_add_pk(transcript_tbl, transcript_id)

  utrome <- utrome %>% 
    dm_add_fk(transcript_tbl, bin_id, bins_tbl) %>% 
    dm_add_fk(transcript_tbl, utr5_id, utr5_tbl) %>% 
    dm_add_fk(transcript_tbl, utr3_id, utr3_tbl) %>% 
    dm_add_fk(transcript_tbl, cds_id, cds_tbl) %>% 
    dm_add_fk(transcript_tbl, aa_id, protein_tbl) %>% 
    dm_add_fk(cds_tbl, aa_id, protein_tbl)

  dm_examine_constraints(utrome) 

  dm_draw_svg(utrome, rankdir = "TB", view_type = "all")
  # dm_draw(utrome)
  # dm_get_all_pks(utrome)
  # dm_get_all_fks(utrome)
  # dm_enum_pk_candidates(utrome, table = bins_tbl)
  # dm_enum_fk_candidates(utrome, table = transcript_tbl, ref_table = bins_tbl)
}