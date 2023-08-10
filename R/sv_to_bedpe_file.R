#' @title SV To BED File.
#'
#' @description Write bedpe format data frame to a file that will work with IGV and UCSC genome browser.
#'
#' @details This function takes four parameters; a data frame with SVs, formatted as bedpe with `sv_df`.
#' The `path` parameter lets the user control the output folder. The default is "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/".
#' `file_name` for specifying the output file name. Results will be written to `results/icgc_dart/misc/`.
#' Lastly, the `add_chr_prefix` lets the user control if chromosomes should be prefixed with "chr" or not.
#' The default is TRUE.
#'
#' @param sv_df data frame of bedpe formatted SV data.
#' @param path The path to the output folder. Default is "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/".
#' @param filename File name (will be written to results/icgc_dart/misc/FILENAME).
#' @param add_chr_prefix Whether to force chr to be added to chromosome names. Default is TRUE.
#'
#' @return bedpe data frame that is compatible with IGN browser.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' SVs_bedpe = sv_to_bedpe_file(sv_df = sv_dataframe,
#'                              filename = "SVs.bedpe",
#'                              add_chr_prefix = TRUE)
#' }
#'
sv_to_bedpe_file = function(sv_df,
                            path = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/",
                            filename = "my_svs.bedpe",
                            add_chr_prefix = TRUE){

  #add chr prefix if missing
  if(add_chr_prefix){
    if(!any(grepl("chr", region_sv$CHROM_A[1]))){
      sv_df = mutate(sv_df, CHROM_A = paste0("chr", CHROM_A)) %>%
        mutate(CHROM_B = paste0("chr", CHROM_B))
    }
  }
  bed_file = paste0(path, filename)
  write.table(sv_df, file = bed_file, sep = "\t", quote = F, row.names = F, col.names = F)
}
