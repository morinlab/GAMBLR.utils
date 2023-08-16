#' @title Cleanup MAF.
#'
#' @description Transform input maf columns to allow for usage of dplyr verbs.
#'
#' @details Transform input maf columns to allow for usage of dplyr verbs.
#' Allowing for a stright-forward plotting workflow as well as downstream data aggregation and manipulation.
#' This function expects a set number of columns to exist in the incoming maf in order for this to work.
#' To view the columns, see bundled data in GAMBLR.data.
#'
#' @param maf_df input MAF data frame.
#'
#' @return maf_df with transformed columns
#'
#' @import dplyr
#' @export
#'
#' @examples
#'
#' clean_maf = cleanup_maf(maf_df = GAMBLR.data::sample_data$grch37$maf)
#'
cleanup_maf = function(maf_df){

  #cleanup various columns that store text to make them useful (make numeric, drop denominator etc)
  maf_df = mutate(maf_df,EXON = gsub("/.+", "", EXON)) %>%
    mutate(EXON = as.numeric(EXON)) %>%
    mutate(INTRON = gsub("/.+", "", INTRON)) %>%
    mutate(INTRON = as.numeric(INTRON)) %>%
    mutate(CDS_position = gsub("/.+", "", CDS_position)) %>%
    mutate(CDS_position = as.numeric(as.character(CDS_position))) %>%
    mutate(cDNA_position = gsub("/.+", "", cDNA_position)) %>%
    mutate(cDNA_position = as.numeric(as.character(cDNA_position))) %>%
    mutate(Protein_position = gsub("/.+", "", Protein_position)) %>%
    mutate(Protein_position = as.numeric(as.character(Protein_position)))

  return(maf_df)
}
