#' Loads the reference genome into memory (ENSEMBL)
#'
#' The different versions can be downloaded
#' \href{http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/}{here}. If the
#' reference transcriptome is not from ENSEBML, it tries to use the
#' \code{loadGTF} function to load a generic transcriptome.
#'
#' @param gtf_path path to the reference genome .gtf file.
#' @param overwrite boolean, whether to overwrite previously reference
#'   transcriptome. Defaults to FALSE.
#'
#' @return the connection to the reference genome DB.
#' @export
loadEdb <- function(gtf_path, overwrite = F) {
  if (!exists("TxDb_ref") | overwrite) {
    logger::log_info("Loading the reference genome.")
    TxDb_ref <<- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    TxDb_ref <<- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  } else {
    logger::log_info("Variable 'TxDb_ref' already loaded!")
  }

  return(TxDb_ref)
}


#' Loads the reference genome into memory (generic)
#'
#' @param gtf_path path to the reference genome .gtf file.
#' @param annotation string, name of the source of the annotation file. Only
#'   supported annotation types are "Ensembl" and "Gencode", and any other imput
#'   might no work as intended. "Gencode" is the recommended annotation for any
#'   other type than "Ensembl". Defaults to "Gencode".
#' @param overwrite boolean, whether to overwrite previously reference
#'   transcriptome. Defaults to FALSE.

#' @return TxDb-class object.
#' @export
loadGTF <- function(gtf_path, annotation = "Gencode", overwrite = F){
  if(!exists("TxDb_ref") | overwrite){
    if(annotation == "Ensembl"){
      TxDb_ref <<- loadEdb(gtf_path, overwrite)
    }else{
      logger::log_info("Loading the reference genome.")
      TxDb_ref <<- GenomicFeatures::makeTxDbFromGFF(gtf_path)
    }
  }else{
    logger::log_info("Variable 'TxDb_ref' already loaded!")
  }

  return(TxDb_ref)
}


#' Loads the ENCODE blacklisted regions into memory
#'
#' The different versions can be downloaded
#' \href{https://github.com/Boyle-Lab/Blacklist/tree/master/lists}{here}.
#'
#' @param blacklist_path path to the ENCODE blacklisted regions .bed file.
#'
#' @return the ENCODE blacklisted region in GRanges object.
#' @export
loadEncodeBlacklist <- function(blacklist_path) {
  if (!exists("encode_blacklist_hg38")) {
    logger::log_info("\t\t Loading the v2 ENCODE blacklisted regions.")
    encode_blacklist_hg38 <<- rtracklayer::import(blacklist_path) %>% GenomeInfoDb::`seqlevelsStyle<-`("NCBI")
  } else {
    logger::log_info("\t\t Variable 'encode_blacklist_hg38' is already loaded!")
  }

  return(encode_blacklist_hg38)
}

#' Symmetric difference + exclusives
#'
#' Obtains the elements in X not in Y, the elements in Y not in X, and the
#' elements exclusive to each set combined (also called symmetric difference).
#' Alternative approach to `dplyr::symdiff`, which only returns the last
#' element.
#'
#' @param x character vector.
#' @param y character vector.
#'
#' @return list of different exclusive elements: X not in Y, Y not in X and the
#'   symmetric difference between both.
#' @export
setdiff_sym <- function(x, y){
  diff1 <- setdiff(x, y)
  diff2 <- setdiff(y, x)
  diff_merge <- union(diff1, diff2)

  return(list("x.ni.y" = diff1, "y.ni.x" = diff2, "sym.diff" = diff_merge))
}
