#' @title Convert mutation data to a shareable format.
#'
#' @description `sanitize_maf_data` returns an oncomatrix of patient/gene data indicating only data needed to produce an oncoplot.
#'
#' @details Write an oncomatrix from a MAF File for further plotting. This is meant to be run by individuals who have access to data sets to
#' "sanitize" a subset of data for subsequent use by them or others who don't have permission to access the raw data.
#' Example: User J has full permissions for ICGC data and has read permissions on a MAF file. User B needs to make some oncoplots
#' and/or perform some statistical analysis on the frequency and assortment of mutations in that data set but don't need all the details.
#' User J can run this function on a maf file and provide the path of the output to user B.
#'
#' @param mutation_maf_path Provide either the full path to a MAF file.
#' @param mutation_maf_data Otherwise provide a data frame of the MAF data.
#' @param output_oncomatrix Optionally provide the path for your sanitized output file (otherwise it writes to the working directory).
#' @param genes_keep Specify which genes you want to remain in the output. Make sure there are no duplicated elements in the vector.
#' @param genes_drop Optionally specify which genes to drop (this doesn't mean all other genes will remain. Maftools decides that part).
#'
#' @return The full path to the oncomatrix file (a matrix with Variant_Classification or Multi_Hit indicating coding mutation status per patient).
#'
#' @importFrom utils write.table
#' @import dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' safe_oncomatrix_path <- sanitize_maf_data(
#'      mutation_maf_data = GAMBLR.data::sample_data$grch37$maf,
#'      genes_keep = c(
#'          "MYC", "ID3", "ARID1A", "FOXO1", "TP53", "FAT4", "IGLL5"
#'      )
#' )
#'
sanitize_maf_data = function(mutation_maf_path,
                             mutation_maf_data,
                             output_oncomatrix,
                             genes_keep,
                             genes_drop = c()){

  if(missing(mutation_maf_path) & missing(mutation_maf_data)){

    warning("Provide either a path to a MAF file or a data frame of mutations")
    return()
  }
  if(!missing(mutation_maf_path)){
    mutation_maf_data = fread_maf(mutation_maf_path)
  }
  #because we use oncoplot we need to ensure the user gets an oncomatrix for all the genes they care about
  #optionally we also exclude genes
  if(missing(genes_keep)){
    warning("you should provide a list of genes to retain in the output to ensure your genes of interest are included")
    genes_keep = rev(unique(lymphoma_genes$Gene))
  }
  maf_o = as.data.frame(mutation_maf_data)
  onco_matrix <- create_onco_matrix(
    maf_df = maf_o,
    genes = genes_keep
  ) %>% as.data.frame
  #writes to working directory
  write.table(
    onco_matrix,
    file = paste0(getwd(), "/onco_matrix.txt"),
    quote = F,
    sep = "\t"
  )
  if(!missing(output_oncomatrix)){
    #rename it
    file.rename("onco_matrix.txt", output_oncomatrix)
  }else{
    output_oncomatrix = paste0(getwd(), "/onco_matrix.tsv")
  }
  message(paste("your data is in:", output_oncomatrix))
  return(output_oncomatrix)
}
