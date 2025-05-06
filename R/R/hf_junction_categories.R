# Variables needed for the junction categories
category_names <- c(
  "annotated" = "Annotated",
  "novel_acceptor" = "Novel Acceptor",
  "novel_donor" = "Novel Donor",
  "diff_acceptor_donor" = "Diff. Acceptor/Donor",
  "novel_exon_skip" = "Novel Exon Skip",
  "novel_combo" = "Novel Combo",
  "ambig_gene" = "Ambig. genes",
  "unannotated" = "Unannotated"
)

category_names_red <- c(
  "annotated" = "Annotated",
  "novel_acceptor" = "Partially Annotated",
  "novel_donor" = "Partially Annotated",
  "novel_exon_skip" = "Partially Annotated",
  "novel_combo" = "Partially Annotated",
  "ambig_gene" = "Ambig. genes",
  "unannotated" = "Unannotated"
)

category_protein_coding <- c("PC", "Other", "Non-PC")


#' Load splice junctions from first pass STAR alignment
#'
#' Source:
#' https://github.com/RHReynolds/LBD-seq-bulk-analyses/blob/main/R/load_sj_df.R
#'
#' @param metadata dataframe, contains the relevant information of all samples.
#'   It is required to have a "sample_id" field with a name that can be matched
#'   with the SJ.out.tab file.
#' @param sj_path character vector, path to the directory containing the splice
#'   junctions files (SJ.out.tab).
#' @param alternative_sample_id, string, if a different sample ID naming is
#'   employed between the metadata dataframe and the SJ.out.tab files, you can
#'   employ this argument to use that field as the Sj.out.tab file identifier.
#' @param overwrite boolean, whether to overwrite previously generated results
#'   from the function. If set to FALSE, the function looks for the files in
#'   disk and loads them if possible. By default, FALSE.
#' @param output_path string, path to where to store the results. Defaults to
#'   empty string so that no result is stored in disk.
#'
#' @return dataframe with all the junctions from the SJ.out.tab files.
#' @export
loadSJ <- function(metadata,
                   sj_path,
                   alternative_sample_id = NULL,
                   overwrite = FALSE,
                   output_path = "") {
  # Check if output already exists and the "overwrite" parameter. If the file is
  # found in disk and "overwrite" is set to F, then we read the previously
  # generated output of the function.
  if (file_test("-f", output_path) & !overwrite) {
    if (grepl("\\.fst$", output_path)) {
      sj_df <- fst::read.fst(output_path) %>% tibble::as_tibble()
    } else if (output_path != "") {
      sj_df <- readRDS(output_path) %>% tibble::as_tibble()
    }

    return(sj_df)
  }

  # If a different sample ID naming is employed between the metadata dataframe
  # and the SJ.out.tab files, you can employ the "alternative_sample_id"
  # argument to use that field as the SJ.out.tab files identifier.
  if(!is.null(alternative_sample_id)){
    sample_id_dictionary <- metadata %>% dplyr::select(all_of(c("sample_id", alternative_sample_id))) %>% tibble::deframe()
    metadata <- metadata %>%
      dplyr::mutate(sample_id = stringr::str_replace_all(sample_id, sample_id_dictionary))
  }

  # Match each sample ID in the metadata dataframe to a file in the folder
  # provided.
  sj_files <- list.files(sj_path, full.names = T)
  sj_sample_files <- foreach(i = 1:nrow(metadata)) %do%
    {
      # Use the metadata dataframe to locate the files.
      row <- metadata[i, ]

      sample_id <- row$sample_id
      sj_sample_path <- sj_files[grepl(sample_id, sj_files)]
      return(tibble::tibble(sj_sample_path, sample_id))
    } %>% dplyr::bind_rows()

  # Correct the sample IDs to match those in the metadata
  if(!is.null(alternative_sample_id)){
    sj_sample_files <- sj_sample_files %>%
      dplyr::mutate(sample_id = stringr::str_replace_all(sample_id, setNames(names(sample_id_dictionary), sample_id_dictionary)))
  }

  # Using the package "vroom", we can read all SJ.out.tab files with one call,
  # applying the same filters and merging the rows into a single dataframe.
  # Additionally, we add the "sample_id" field to each row so that we can keep
  # track of the junction's origin. The strand is automatically converted into
  # "*", "+" and "-".
  sj_df <- vroom::vroom(sj_sample_files$sj_sample_path,
    delim = "\t",
    col_names = c(
      "chr", "intron_start", "intron_end", "strand",
      "intron_motif", "intron_annotation",
      "unique_reads_junction", "multimap_reads_junction",
      "max_splice_alignment_overhang"
    ),
    col_types = "cdddddddd",
    id = "sample_id",
    progress = F
  ) %>%
    dplyr::filter(unique_reads_junction > 0) %>%
    dplyr::mutate(sample_id = unname(deframe(sj_sample_files)[sample_id])) %>%
    dplyr::mutate(strand = case_when(strand == 0 ~ "*",
      strand == 1 ~ "+",
      strand == 2 ~ "-",
      .default = as.character(strand)
    ))

  # Add junction id
  sj_df <- sj_df %>%
    dplyr::group_by(chr, intron_start, intron_end, strand) %>%
    dplyr::mutate(junction_id = dplyr::cur_group_id()) %>%
    dplyr::relocate(sample_id, junction_id) %>%
    dplyr::ungroup()

  # If an output path is provided, we write in disk the resulted dataframe.
  if (grepl("\\.fst$", output_path)) {
    sj_df %>% fst::write.fst(output_path)
  } else if (output_path != "") {
    sj_df %>% saveRDS(output_path)
    Sys.sleep(1)
  }

  return(sj_df)
}


#' Annotates junctions from first pass STAR alignment
#'
#' Source:
#' https://github.com/RHReynolds/LBD-seq-bulk-analyses/blob/main/R/create_rse_jx.R#L16
#'
#' @param sj_df dataframe, contains all the junctions extracted using
#'   \code{loadSJ}.
#' @param gtf_path character vector, path to the reference transcriptome.
#' @param annotation string, name of the source of the annotation file. Only
#'   supported annotation types are "Ensembl" and "Gencode", and any other input
#'   might no work as intended. "Gencode" is the recommended annotation for any
#'   other type than "Ensembl". Defaults to "Gencode".
#' @param prune_output boolean, whether to clean the output object from
#'   non-relevant columns. Defaults to FALSE.
#' @param overwrite boolean, whether to overwrite previously generated results
#'   from the function. If set to FALSE, the function looks for the files in
#'   disk and loads them if possible. By default, FALSE.
#' @param output_path string, path to where to store the results. Defaults to
#'   empty string so that no result is stored in disk.
#'
#' @return GRanges object.
#' @export
annotateSJ <- function(sj_df,
                       gtf_path,
                       annotation = "Ensembl",
                       prune_output = F,
                       overwrite = F,
                       output_path = "") {
  # Check if output already exists and the "overwrite" parameter. If the file is
  # found in disk and "overwrite" is set to F, then we read the previously
  # generated output of the function.
  if (file_test("-f", output_path) & !overwrite) {
    sj_annotated <- readRDS(output_path)
    return(sj_annotated)
  }

  # Load reference transcriptome and other function variables
  TxDb_ref <- loadGTF(gtf_path, annotation)
  GenomeInfoDb::seqlevelsStyle(TxDb_ref) <- "NCBI"

  valid_seqnames <- c(str_c(1:22), "X", "Y")

  # A custom version of dasper is employed that allows for multi-processing
  # annotation. To active this functionality and keep compatibility with
  # previous versions, the function looks for the variable "dasper_cores" in the
  # current environment. If you use the latest public version of dasper, then
  # this variable is simply ignored.

  # dasper_cores <- 4

  # Annotate the introns
  sj_annotated <- sj_df %>%
    dplyr::distinct(junction_id, chr, intron_start, intron_end, strand) %>%
    GenomicRanges::GRanges() %>%
    GenomeInfoDb::`seqlevelsStyle<-`("NCBI") %>%
    GenomeInfoDb::keepSeqlevels(value = valid_seqnames, pruning.mode = "tidy") %>%
    dasper::junction_annot(
      ref = TxDb_ref,
      ref_cols = c("gene_id", "tx_name", "exon_id"),
      ref_cols_to_merge = c("gene_id", "tx_name")
    )

  # Remove some columns
  if (prune_output) {
    relevant_fields <- c("junction_id", "in_ref", "gene_id_junction", "tx_name_junction", "strand_junction", "type")
    sj_annotated <- sj_annotated[, relevant_fields]
  }

  # If an output path is provided, we write in disk the resulted dataframe.
  if (output_path != "") sj_annotated %>% saveRDS(output_path)

  return(sj_annotated)
}


#' Prunes annotated junctions
#'
#' Once the junctions are annotated, this function applies the pruning process
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#' @param gtf_path string, path to where the reference transcriptome can be
#'   found on disk.
#' @param blacklist_path GRanges, contains the blacklisted regions.
#' @param overwrite boolean, whether to overwrite previously generated results
#'   from the function. If set to FALSE, the function looks for the files in
#'   disk and loads them if possible. By default, FALSE.
#' @param output_path string, path to where to store the results. Defaults to
#'   empty string so that no result is stored in disk.
#'
#' @return Pruned GRanges object.
#' @export
pruneSJ <- function(sj_annotated,
                    gtf_path,
                    blacklist_path,
                    overwrite = F,
                    output_path = ""){
  # Check if output already exists and the "overwrite" parameter. If the file is
  # found in disk and "overwrite" is set to F, then we read the previously
  # generated output of the function.
  if (file_test("-f", output_path) & !overwrite) {
    sj_annotated_pruned <- readRDS(output_path)
    return(sj_annotated_pruned)
  }

  # Load the ENCODE blacklisted region
  encode_blacklist_hg38 <- loadEncodeBlacklist(blacklist_path)

  # 1. Remove ENCODE blacklisted regions.
  # 2. Add a column to identify ambiguous genes.
  # 3. Remove short junctions.
  # 4. Add biotype percentage.
  sj_annotated_pruned <- removeEncodeBlacklistRegionsSJ(sj_annotated, encode_blacklist_hg38) %>%
    addColumnAmbiguousGenesSJ() %>%
    removeShortJunctionsSJ() %>%
    addBiotypePercentageSJ(gtf_path)

  # If an output path is provided, we write in disk the resulted dataframe.
  if (output_path != "") sj_annotated_pruned %>% saveRDS(output_path)

  return(sj_annotated_pruned)
}


#' Remove the junctions from the ENCODE blacklisted regions
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#' @param encode_blacklist_hg38 GRanges, contains the blacklisted regions.
#'
#' @return GRanges with only the junctions that do not overlap with the
#'   blacklisted regions.
#' @export
removeEncodeBlacklistRegionsSJ <- function(sj_annotated, encode_blacklist_hg38) {
  ## Look fot the overlaps between sj_annotated and the ENCODE blacklisted region
  overlaps <- GenomicRanges::findOverlaps(
    query = encode_blacklist_hg38,
    subject = sj_annotated,
    ignore.strand = F,
    type = "any"
  )

  idxs <- S4Vectors::subjectHits(overlaps)

  ## If an overlap is found, remove the junctions
  if(length(idxs) > 0) sj_annotated <- sj_annotated[-idxs, ]

  return(sj_annotated)
}


#' Remove the junctions from ambiguous genes.
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#'
#' @return GRanges with only the junctions associated to one or less genes.
#' @export
removeAmbiguousGenesSJ <- function(sj_annotated) {
  idxs <- which(sapply(sj_annotated$gene_id_junction, length) >= 2)
  if(length(idxs) > 0) sj_annotated <- sj_annotated[-idxs, ]

  return(sj_annotated)
}


#' Add field to specify if the junction is from an ambiguous gene.
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#'
#' @return GRanges with an added column for ambiguous genes.
#' @export
addColumnAmbiguousGenesSJ <- function(sj_annotated) {
  idxs <- which(sapply(sj_annotated$gene_id_junction, length) >= 2)

  GenomicRanges::mcols(sj_annotated)[["ambig_gene"]] <- F
  GenomicRanges::mcols(sj_annotated)[idxs, "ambig_gene"] <- T

  return(sj_annotated)
}


#' Remove the junctions shorter than 25bp.
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#'
#' @return GRanges with only the junctions bigger than 25 bp.
#' @export
removeShortJunctionsSJ <- function(sj_annotated) {
  idxs <- which(sj_annotated %>% BiocGenerics::width() < 25)
  if(length(idxs) > 0) sj_annotated <- sj_annotated[-idxs, ]

  return(sj_annotated)
}


#' Add biotype percentage
#'
#' Load the reference transcriptome biotype information for each transcript and
#' calculates the percentage of protein-coding transcripts for every junction.
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#' @param gtf_path string, path to where the reference transcriptome can be
#'   found on disk.
#'
#' @return GRanges object with the biotype percentage added.
#' @export
addBiotypePercentageSJ <- function(sj_annotated, gtf_path){
  # Not sure how to load the transcripts biotype and transcripts IDs directly
  # into TxDb_ref. As of now, only using rtracklayer::import has I been able to
  # load the transcripts.

  transcript_biotype_df <- rtracklayer::import(gtf_path) %>%
    tibble::as_tibble() %>%
    `colnames<-`(gsub("transcript_type", "transcript_biotype", colnames(.))) %>%
    dplyr::distinct(transcript_id, transcript_biotype, gene_id) %>%
    dplyr::filter(!is.na(transcript_id)) %>%
    dplyr::mutate(transcript_id = stringr::str_remove(transcript_id, "\\..+"),
                  gene_id = stringr::str_remove(gene_id, "\\..+")) %>%
    dplyr::distinct(transcript_id, .keep_all = T) %>%
    dplyr::select(tx_name = transcript_id,
                  tx_biotype = transcript_biotype,
                  tx_gene_id = gene_id)

  df_introns <- sj_annotated %>%
    tibble::as_tibble() %>%
    dplyr::select(seqnames, start, end, strand, junction_id, tx_name_junction, gene_id_junction) %>%
    tidyr::unnest(tx_name_junction) %>%
    dplyr::mutate(tx_name_junction = stringr::str_remove(tx_name_junction, "\\..+")) %>%
    dplyr::left_join(transcript_biotype_df, by = c("tx_name_junction" = "tx_name"))

  ## Calculate percentage of all junctions (uncomment the filter to only use unambiguous junctions)
  df_all_percentage <- df_introns %>%
    dplyr::group_by(junction_id) %>%
    # dplyr::filter(n_distinct(tx_gene_id) == 1)  %>%
    dplyr::group_by(seqnames, start, end, strand, junction_id, tx_biotype) %>%
    dplyr::distinct(tx_name_junction, .keep_all = T) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(percent = (n/sum(n))*100) %>%
    dplyr::ungroup()

  ## Filter only protein coding
  df_all_percentage_merged <- df_all_percentage %>%
    dplyr::group_by(junction_id) %>%
    tidyr::pivot_wider(id_cols = c(seqnames, start, end, strand, junction_id),
                       names_from = tx_biotype, values_from = percent) %>%
    dplyr::select(seqnames, start, end, strand, junction_id,
                  protein_coding = protein_coding) %>%
    replace(is.na(.), 0) %>%
    dplyr::ungroup()

  gr_all_percentages <- df_all_percentage_merged %>%
    GenomicRanges::GRanges()

  overlaps <- GenomicRanges::findOverlaps(query = gr_all_percentages, subject = sj_annotated,
                                          type = "equal")

  S4Vectors::mcols(sj_annotated)$protein_coding <- NA
  S4Vectors::mcols(sj_annotated)$protein_coding[S4Vectors::subjectHits(overlaps)] <- S4Vectors::mcols(gr_all_percentages)$protein_coding[S4Vectors::queryHits(overlaps)]

  return(sj_annotated)
}

#' Calculates the proportions of junctions/reads for each category
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#' @param sj_df dataframe with all the junctions from the SJ.out.tab files.
#' @param metadata dataframe, contains additional information for each sample.
#' @param category_names named list of string. The values must be the final
#'   names of the category and the names must be the values found in the
#'   sj_annotated object.
#' @param overwrite boolean, whether to overwrite previously generated results
#'   from the function. If set to FALSE, the function looks for the files in
#'   disk and loads them if possible. By default, FALSE.
#' @param output_path string, path to where to store the results. Defaults to
#'   empty string so that no result is stored in disk.
#'
#' @return dataframe with the proportions for each sample and category.
#' @export
generateJunctionProportions <- function(sj_annotated,
                                        sj_df,
                                        metadata,
                                        category_names,
                                        overwrite = F,
                                        output_path = ""){
  # Check if output already exists and the "overwrite" parameter. If the file is
  # found in disk and "overwrite" is set to F, then we read the previously
  # generated output of the function.
  if(file_test("-f", output_path) & !overwrite){
    annotated_prop_df <- readRDS(output_path)
    return(annotated_prop_df)
  }

  # 1. Extract the junctions IDs and the junction categories from
  # "sj_annotated".
  #
  # 2. Add an entry of each junction ID for each entry in the "sj_df" dataframe.
  # That is, one row per occurrence of the junction ID in the extracted
  # junctions. This also adds the number of reads for that particular sample and
  # junction.
  #
  # 3. Grouping by sample, calculate the total number of junctions, the total
  # number of counts and the proportion by each type.
  #
  # 4. Add metadata information to each sample and category for later
  # visualizations.
  annotated_prop_df <- S4Vectors::mcols(sj_annotated)[c("junction_id", "type")] %>%
    tibble::as_tibble() %>%
    dplyr::mutate(type = factor(category_names[type %>% as.character()], levels = unique(category_names))) %>%
    dplyr::left_join(sj_df %>% dplyr::select(sample_id, junction_id, reads = unique_reads_junction),
                     by = "junction_id") %>%
    dplyr::group_by(sample_id, type) %>%
    dplyr::summarise(n_junc = n(), n_count = sum(reads)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(total_junc = sum(n_junc),
                  total_count = sum(n_count),
                  prop_junc = n_junc/total_junc,
                  prop_count = n_count/total_count) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(metadata, by = "sample_id")

  # If an output path is provided, we write in disk the resulted dataframe.
  if(output_path != "") annotated_prop_df %>% saveRDS(output_path)

  return(annotated_prop_df)
}


#' Calculates the proportions of junctions/reads for each category
#'
#' Specific for the protein-coding categories.
#'
#' @param sj_annotated GRanges, contains the annotated information for each
#'   junction.
#' @param sj_df dataframe with all the junctions from the SJ.out.tab files.
#' @param metadata dataframe, contains additional information for each sample.
#' @param overwrite boolean, whether to overwrite previously generated results
#'   from the function. If set to FALSE, the function looks for the files in
#'   disk and loads them if possible. By default, FALSE.
#' @param output_path string, path to where to store the results. Defaults to
#'   empty string so that no result is stored in disk.
#'
#' @return dataframe with the proportions for each sample and category.
#' @export
generateJunctionProportionsPC <- function(sj_annotated,
                                          sj_df,
                                          metadata,
                                          category_protein_coding,
                                          overwrite = F,
                                          output_path = ""){
  # Check if output already exists and the "overwrite" parameter. If the file is
  # found in disk and "overwrite" is set to F, then we read the previously
  # generated output of the function.
  if(file_test("-f", output_path) & !overwrite){
    annotated_prop_df <- readRDS(output_path)
    return(annotated_prop_df)
  }

  # 1. Extract the junctions IDs and the junction categories from
  # "sj_annotated".
  #
  # 2. Add an entry of each junction ID for each entry in the "sj_df" dataframe.
  # That is, one row per occurrence of the junction ID in the extracted
  # junctions. This also adds the number of reads for that particular sample and
  # junction.
  #
  # 3. Grouping by sample, calculate the total number of junctions, the total
  # number of counts and the proportion by each type.
  #
  # 4. Add metadata information to each sample and category for later
  # visualizations.
  annotated_prop_df <- S4Vectors::mcols(sj_annotated)[c("junction_id", "type", "protein_coding")] %>%
    tibble::as_tibble() %>%
    dplyr::filter(type == "annotated") %>%
    dplyr::left_join(sj_df %>% dplyr::select(sample_id, junction_id, reads = unique_reads_junction),
                     by = "junction_id") %>%
    dplyr::mutate(protein_coding = case_when(protein_coding == 100 ~ "PC",
                                             protein_coding == 0 ~ "Non-PC",
                                             .default = "Other")) %>%
    dplyr::mutate(protein_coding = factor(protein_coding, levels = category_protein_coding)) %>%
    dplyr::group_by(sample_id, protein_coding) %>%
    dplyr::summarise(n_junc = n(), n_count = sum(reads)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(total_junc = sum(n_junc),
                  total_count = sum(n_count),
                  prop_junc = n_junc/total_junc,
                  prop_count = n_count/total_count) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metadata, by = "sample_id")

  # If an output path is provided, we write in disk the resulted dataframe.
  if(output_path != "") annotated_prop_df %>% saveRDS(output_path)

  return(annotated_prop_df)
}


#' Converts a dataframe to a count matrix
#'
#' @param df dataframe, must have "sample_id" and "type" fields. Additionally,
#'   it must have a field with the name provided in the parameter "field".
#' @param field string, field in the dataframe object that will be the value in
#'   the count matrix.
#' @param category_field, string, field in the dataframe object to split into
#'   different columns. That is, the different categories in our data. Defaults
#'   to "type".
#'
#' @return
#' @export matrix, rows are the categories, columns are the samples and the
#'   values are the field from the parameter "field".
df_to_cts <- function(df, field, category_field = "type"){
  cts <- df %>%
    dplyr::select(sample_id, type = !!category_field, field = !!field) %>%
    tidyr::pivot_wider(id_cols = type, names_from = sample_id, values_from = field) %>%
    tibble::column_to_rownames("type") %>%
    as.matrix()

  return(cts)
}


#' Converts a count matrix object to a dataframe
#'
#' @param cts matrix.
#' @param metadata dataframe, contains additional information for each sample.
#' @param field string, name in the final dataframe where the values from the
#'   count matrix will be stored. Defaults to "n_junc".
#' @param corrected_name string, identifier of the correction method. Defaults
#'   to "Uncorrected".
#'
#' @return dataframe with the "sample_id", the "type" or category and the matrix
#'   values as a column named after the argument "field".
#' @export
cts_to_df <- function(cts,
                      metadata,
                      field = "n_junc",
                      corrected_name = "Uncorrected"){
  cts_df <- cts %>%
    tibble::as_tibble(rownames = "type") %>%
    tidyr::pivot_longer(cols = -type, names_to = "sample_id", values_to = field) %>%
    dplyr::left_join(metadata, by = "sample_id") %>%
    dplyr::mutate(corrected = corrected_name) %>%
    dplyr::relocate(sample_id, type, !!field)

  return(cts_df)
}


#' Extract residuals from lmFit and Crumblr's object
#'
#' @param lmFit MArrayLM, linear model generated from limma.
#' @param cobj Elist, output of the function crumblr::crumblr().
#' @param design Matrix, design of the experiment.
#' @param experimental_condition string, experimental condition being studied.
#'
#' @return
#' @export
extractResidualsCrumblr <- function(lmFit, cobj, design, experimental_condition){
  library(limma)
  library(crumblr)

  # Extract the corrected proportions
  beta <- coef(lmFit)[, grepl(paste(c(experimental_condition, "\\(Intercept\\)"), collapse = "|"), colnames(coef(lmFit)))]
  x <- design[, grepl(paste(c(experimental_condition, "\\(Intercept\\)"), collapse = "|"), colnames(design))]
  resid <- residuals(lmFit, cobj)

  cobj_adjusted <- t(beta %*% t(x) + resid)

  return(cobj_adjusted)
}

#' Apply Limma contrasts
#'
#' @param lmFit MArrayLM, linear model generated from limma.
#' @param design Matrix, design of the experiment.
#' @param experimental_condition string, experimental condition being studied.
#'   It is employed to generate the contrasts. Only valid inputs are: study_arm
#'   and mutation.
#'
#' @return MArrayLM, linear model generated from limma with contrasts.
applyLimmaContrasts_PPMI <- function(lmFit,
                                     design,
                                     experimental_condition){
  # Two different contrasts can be executed based on the experimental condition:
  if(!experimental_condition %in% c("study_arm", "mutation")) stop("The experimental condition is not valid")

  if(experimental_condition == "study_arm"){
    # Boolean to identify the different experiments that are performed
    is_HCPD_experiment <- any(grepl("study_armHealthyControl", colnames(design))) & any(grepl("study_armPD", colnames(design)))
    is_GC_experiment <- any(grepl("study_armGeneticCohortPD", colnames(design))) & any(grepl("study_armGeneticCohortUnaffected", colnames(design)))

    if(is_HCPD_experiment & is_GC_experiment){
      contrasts_fit <- limma::makeContrasts(
        "gcUnaffected - HealthyControl" = study_armGeneticCohortUnaffected - study_armHealthyControl,
        "gcPD - HealthyControl" = study_armGeneticCohortPD - study_armHealthyControl,
        "gcPD - gcUnaffected" = study_armGeneticCohortPD - study_armGeneticCohortUnaffected,
        levels = design)
    }else if(is_HCPD_experiment){
      contrasts_fit <- limma::makeContrasts(
        "PD - HealthyControl" = study_armPD - study_armHealthyControl,
        levels = design)
    }else if(is_GC_experiment){
      contrasts_fit <- limma::makeContrasts(
        "gcPD - gcUnaffected" = study_armGeneticCohortPD - study_armGeneticCohortUnaffected,
        levels = design)
    }
  }

  if(experimental_condition == "mutation"){
    # Boolean to identify the different experiments that are performed
    is_HCPD_mut_experiment <- any(grepl("mutationHealthyControl", colnames(design)))
    is_GC_mut_experiment <- any(grepl("mutationLRRK2Aff", colnames(design))) & any(grepl("mutationLRRK2Unaff", colnames(design)))

    if(is_GC_mut_experiment & is_HCPD_mut_experiment){
      contrasts_fit <- limma::makeContrasts(
        "LRRK2_Unaff - HealthyControl" = mutationLRRK2Unaff - mutationHealthyControl,
        "LRRK2_PD - HealthyControl" = mutationLRRK2Aff - mutationHealthyControl,
        "LRRK2_PD - LRRK2_Unaff" = mutationLRRK2Aff - mutationLRRK2Unaff,
        "SNCA_Unaff - HealthyControl" = mutationSNCAUnaff - mutationHealthyControl,
        "SNCA_PD - HealthyControl" = mutationSNCAAff - mutationHealthyControl,
        "SNCA_PD - SNCA_Unaff" = mutationSNCAAff - mutationSNCAUnaff,
        "GBA_Unaff - HealthyControl" = mutationGBAUnaff - mutationHealthyControl,
        "GBA_PD - HealthyControl" = mutationGBAAff - mutationHealthyControl,
        "GBA_PD - GBA_Unaff" = mutationGBAAff - mutationGBAUnaff,
        levels = design)
    }else if(is_HCPD_mut_experiment){
      contrasts_fit <- limma::makeContrasts(
        "PD - HealthyControl" = mutationSporadic - mutationHealthyControl,
        levels = design)
    }else if(is_GC_mut_experiment){
      contrasts_list <- c()
      if(sum(grepl("LRRK2", colnames(design))) == 2) contrasts_list <- c(contrasts_list, "LRRK2_PD - LRRK2_Unaff" = "mutationLRRK2Aff - mutationLRRK2Unaff")
      if(sum(grepl("GBA", colnames(design))) == 2) contrasts_list <- c(contrasts_list, "GBA_PD - GBA_Unaff" = "mutationGBAAff - mutationGBAUnaff")
      if(sum(grepl("SNCA", colnames(design))) == 2) contrasts_list <- c(contrasts_list, "SNCA_PD - SNCA_Unaff" = "mutationSNCAAff - mutationSNCAUnaff")

      contrasts_fit <- limma::makeContrasts(contrasts = contrasts_list, levels = design) %>% `colnames<-`(names(contrasts_list))
    }
  }

  # Apply the contrast
  lmFit <- limma::contrasts.fit(lmFit, contrasts_fit) %>% limma::eBayes()
}


#' Extract the results from the contrasts
#'
#' @param lmFit_con MArrayLM, linear model generated from limma with contrasts.
#' @param coef string, coefficient to extract the results.
#'
#' @return dataframe with limma results.
#' @export
extractContrastResults <- function(lmFit_con, coef = NULL, rownames = "Type"){
  lmFit_con %>%
    limma::topTable(n = Inf, coef = coef) %>%
    tibble::as_tibble(rownames = rownames) %>%
    dplyr::mutate(coef = coef) %>%
    dplyr::relocate(coef)
}


#' Correct Negative Values (Deprecated)
#'
#' Looks for negative values after correction and set thems to 0. It logs the
#' number of values affected.
#'
#' @param cts_adj matrix, contain the adjusted counts.
#' @param category string, specifies the category that is being studied.
#'   Employed to display more information in the logs. Defaults to NULL.
#'
#' @return matrix with no negative values.
#' @export
correctNegativeValues <- function(cts_adj, category = NULL){
  log_text <- paste0("Number of negative values: ", sum(cts_adj < 0), " out of ", length(cts_adj), ".")
  if(!is.null(category)) log_text <- paste0("(", category, ") ", log_text)
  logger::log_info(log_text)

  cts_adj[cts_adj < 0] <- 0

  return(cts_adj)
}


#' Merge all junction dataframes (Depreca)
#'
#' Given the original junction dataframe with the counts in long format and the
#' adjusted unique junctions and junction counts (also in long format), it
#' combines the three dataframes into a single instance.
#'
#' @param annotated_prop_df dataframe, with the proportions for each sample and category.
#' @param cts_juncs_adj_df dataframe, with the adjusted unique junctions in long format.
#' @param cts_counts_adj_df dataframe, with the adjusted junction counts in long format.
#'
#' @return dataframe with the corrected and not corrected counts in long format.
#' @export
mergeJunctionDF <- function(annotated_prop_df, cts_juncs_adj_df, cts_counts_adj_df){
  annotated_prop_adj_df <- dplyr::inner_join(
    cts_juncs_adj_df %>% dplyr::select(sample_id, type, n_junc),
    cts_counts_adj_df,
    by = c("sample_id", "type")
  ) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(total_junc = sum(n_junc),
                  total_count = sum(n_count),
                  prop_junc = n_junc/total_junc,
                  prop_count = n_count/total_count,
                  .after = n_count)

  prop_df_full <- dplyr::bind_rows(
    annotated_prop_adj_df,
    annotated_prop_df %>% dplyr::mutate(corrected = "Uncorrected")
  ) %>%
    dplyr::mutate(type = factor(type, levels = levels(annotated_prop_df$type))) %>%
    dplyr::relocate(sample_id, corrected) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::ungroup()

  return(prop_df_full)
}