#' @title Annotate copy number matrix by gene names
#' 
#' @description Annotates a copy number matrix using user-provided gene names 
#' 
#' @param these_samples_metadata The metadata for samples of interest to be included in the returned matrix.
#'   Can be created with `get_gambl_metadata` function.
#' @param gene_ids The "gene_id" column stores gene symbols (characters).
#' @param cn_matrix User can provide a matrix of CN values for the samples in the metadata.
#' See [GAMBLR.utils::segmented_data_to_cn_matrix] for more information on how to create this matrix.
#' 
#' @return A copy number matrix with sample ids as rows and gene names as columns
#' 
#' @import dplyr 
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' # Get sample metadata 
#' dlbcl_meta = get_gambl_metadata() %>% filter(pathology == "DLBCL")
#' 
#' # Consider some genes
#' gene_symbols = c("MYC", "MIR17HG", "CCND3","ID3","DDX3X", "SYNCRIP")
#' 
#' # Get cn_matrix
#' all_segments = get_cn_segments(these_samples_metadata = dlbcl_meta)
#' cn_matrix = segmented_data_to_cn_matrix(
#'    seg_data = all_segments,
#'     strategy="auto_split",
#'    these_samples_metadata = dlbcl_meta,
#'    adjust_for_ploidy=TRUE
#' )
#' 
#' annotate_cn_matrix_by_genes (these_samples_metadata = dlbcl_meta,
#'                              gene_ids = gene_symbols,
#'                              cn_matrix = cn_matrix
#' )
#' }

annotate_cn_matrix_by_genes = function(these_samples_metadata,
                                       gene_ids,
                                       cn_matrix){
  
  if(missing(these_samples_metadata)){
    if(verbose){
      print("missing these_samples_metadata")
    }
  }
  
  gene_bins = unlist(map_regions_to_bins(
    query_regions = gene_ids,
    regions = colnames(cn_matrix),
    query_type = "gene",
    first = TRUE
  ))
  
  cn_matrix = cn_matrix[these_samples_metadata$sample_id,gene_bins,drop = FALSE] %>% as.data.frame()
  
  colnames(cn_matrix) = names(gene_bins)
  
  return(cn_matrix)
}