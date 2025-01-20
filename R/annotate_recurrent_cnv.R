#' @title Annotate Recurrent CNVs.
#' 
#' @description Annotates recurrent CNVs from a data frame with CNV data.
#'
#' @details This function takes a data frame with CNVs (`seq_df`) and annotates recurrent CNVs.
#'
#' @param seg_df Data frame of sequences with start and end coordinates.
#' @param seg_file Optional argument to read sequences from file (currently not used in function).
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' my_segs = get_sample_cn_segments(this_sample_id = "HTMCP-01-06-00422-01A-01D")
#' annotated = annotate_recurrent_cnv(seg_df = my_segs)
#' }
#' 
annotate_recurrent_cnv = function(seg_df,
                                  seg_file){

  cnv_coord_df = data.frame(chrom = c("18"),
                            start = c(60000000),
                            end = c(61000000),
                            anchor = c("right"),
                            min_cn = c(3),
                            name = "18der_gain")

}
