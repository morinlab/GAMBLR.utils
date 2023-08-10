#' @title Region To Bins.
#'
#' @description Split a contiguous genomic region on a chromosome into non-overlapping bins
#'
#' @details This function takes genomic coordinates with the `chromosome`, `start`, and `end` parameters.
#' Lastly, the user can also specify the bin size with `bin_size`. Default is 20000
#'
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param start Query start coordinate of the range you are restricting to.
#' @param end Query end coordinate of the range you are restricting to.
#' @param bin_size The size of the bins, default is 2000.
#'
#' @return Data frame describing the bins various ways.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' chr8q_bins = region_to_bins(chromosome = "8",
#'                             start = 48100000,
#'                             end = 146364022,
#'                             bin_size = 20000)
#'
region_to_bins = function(chromosome,
                 start,
                 end,
                 bin_size = 2000){

  if(missing(chromosome)){
    stop("Please provide a chromosome...")
  }

  if(missing(start)){
    stop("Please provide the start coordinates...")
  }

  if(missing(end)){
    stop("Please provide the end coordinates...")
  }

  bin_df = data.frame(bin_chr = chromosome, bin_start = seq(start, end, bin_size))
  bin_df = mutate(bin_df,bin_end = bin_start+ bin_size) %>%
    dplyr::filter(bin_end<=end) %>%
    mutate(region = paste0(bin_chr,":",bin_start,"-",bin_end))

  return(bin_df)
}
