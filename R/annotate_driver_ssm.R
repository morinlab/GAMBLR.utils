#' @title Annotate Driver SSM.
#'
#' @description Retrieve data or use the supplied data frame of genes to get/annotate MAF data indicating which rows are putative driver mutations.
#'
#' @details Provide a maf-like data frame with `maf_df` as the only required parameter. 
#' For information on how to use the additional parameters, refer to the parameter descriptions and function examples. 
#'
#' @param maf_df Data frame of MAF-format mutations.
#' @param lymphoma_type Optional keyword to find genes for annotation (e.g. "BL","DLBCL","FL") is not specifying driver genes.
#' @param driver_genes Optional vector of Hugo_Symbol of genes for coding mutations.
#' @param noncoding_regions Optional named vector of regions to use to further restrict noncoding mutations per gene.
#'
#' @return Incoming MAF, subset to putative driver mutations.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' # Get MAF-format data frame 
#' my_maf <- GAMBLR.data::sample_data$grch37$maf
#' 
#' # If you have GSC-restricted access, you can get your MAF-format data frame
#' # by using `GAMBLR.results`
#' my_meta <- GAMBLR.results::get_gambl_metadata()
#' my_maf <- GAMBLR.results::get_ssm_by_samples(
#'   these_samples_metadata = my_meta,
#'   projection = "grch37"
#' )
#' 
#' # Annotate driver SSM
#' driver_ssm_b = annotate_driver_ssm(
#'   maf_df = my_maf,
#'   lymphoma_type = "DLBCL",
#'   noncoding_regions=c("NFKBIZ"="chr3:101578206-101578365")
#' )
#'
annotate_driver_ssm = function(maf_df,
                               lymphoma_type,
                               driver_genes,
                               noncoding_regions = c("NFKBIZ" = "chr3:101578206-101578365", "HNRNPH1" = "chr5:179,045,946-179,046,683")){

  #get the gene list if missing
  if(missing(driver_genes)){
    driver_genes = GAMBLR.data::lymphoma_genes %>%
      pull(Gene)

    ngene = length(driver_genes)
    message(paste("using", ngene, "genes"))
  }
  if(missing(maf_df)){
    warning("Please specify mutations data frame (MAF-format)")
    return()
  }else{
    #this all assumes the user has left the noncoding variants in their MAF
    kept_ssm = maf_df %>%
      dplyr::filter(Variant_Classification %in% coding_vc) %>%
      dplyr::filter(Hugo_Symbol %in% driver_genes)

    for(gene in names(noncoding_regions)){
      message(paste("adding", unname(noncoding_regions[gene]), "for", gene))
      nc_ssm = maf_df %>%
      dplyr::filter(Variant_Classification == unname(noncoding_regions[gene]))

      if(!is.na(noncoding_regions[gene])){
        #also restrict to coordinates
        chunks = region_to_chunks(noncoding_regions[gene])
        nc_ssm = dplyr::filter(nc_ssm,Start_Position >= chunks$start & Start_Position <= chunks$end)
      }
      kept_ssm = bind_rows(kept_ssm, nc_ssm)
      print(dim(kept_ssm))
    }
  }
  return(kept_ssm)
}
