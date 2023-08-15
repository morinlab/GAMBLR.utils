#' @title Count SSM In A Region
#'
#' @description Count the variants in a region with a variety of filtering options.
#'
#'
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param start Query start coordinate of the range you are restricting to.
#' @param end Query end coordinate of the range you are restricting to.
#' @param these_samples_metadata A metadata table subset to the sample IDs of interest. If not provided, the function will return all metadata and regions will be returned for all samples in the metadata.
#' @param all_mutations_in_these_regions If you are calling this function many times (e.g. bins spanning a larger region), to save a ton of time you are strongly encouraged to provide the output of `get_ssm_by_region` on the entire region of interest and passing it to this function
#' @param count_by Defaults to counting all variants. Specify 'sample_id' if you want to collapse and count only one per sample
#' @param seq_type The seq_type you want back, default is genome.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' #define a region.
#' my_region = gene_to_region(gene_symbol = "MYC",
#'                            return_as = "region")
#'
#' #get meta data and subset
#' my_metadata = GAMBLR.data::gambl_metadata
#' fl_metadata = dplyr::filter(my_metadata, pathology == "FL")
#'
#' #count SSMs for the selected sample subset and defined region.
#' fl_ssm_counts_myc = count_ssm_by_region(region = my_region,
#'                                        these_samples_metadata = fl_metadata)
#'
count_ssm_by_region = function(region,
                               chromosome,
                               start,
                               end,
                               all_mutations_in_these_regions,
                               these_samples_metadata,
                               count_by,
                               seq_type = "genome"){

  if(missing(these_samples_metadata)){
    these_samples_metadata = GAMBLR.helpers::handle_metadata(this_seq_type = seq_type)
  }
  if(!missing(all_mutations_in_these_regions)){
    # function was provided the mutations already so we just need to subset it to the region of interest

    region_muts = dplyr::filter(all_mutations_in_these_regions,Start_Position >= start, Start_Position < end)
  }else if(missing(region)){
    region_muts = GAMBLR.helpers::handle_ssm_by_region(chromosome=chromosome,qstart=start,qend=end,streamlined = TRUE)
  }else{
    region_muts = GAMBLR.helpers::handle_ssm_by_region(region=region,streamlined = TRUE)
  }
  keep_muts = dplyr::filter(region_muts,Tumor_Sample_Barcode %in% these_samples_metadata$Tumor_Sample_Barcode)
  if(missing(count_by)){
    #count everything even if some mutations are from the same patient
    return(nrow(keep_muts))
  }else if(count_by == "sample_id"){
    return(length(unique(keep_muts$Tumor_Sample_Barcode)))
  }else{
    print("Not sure what to count")
  }
}
