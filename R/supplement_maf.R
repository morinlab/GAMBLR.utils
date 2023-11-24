#' @title Supplement MAF.
#'
#' @description Complement maf with missing samples.
#'
#' @details Specify the initial MAF with `incoming_maf` (to be supplemented with missing samples) and
#' give the function a filtered metadata table (with the sample IDs of interest) to the `these_samples_metadata`.
#'
#' @param incoming_maf The initial MAF data frame to be supplemented with missing samples.
#' @param these_samples_metadata The metadata data frame that contains Tumor_Sample_Barcode column with ids to be present in the complemented MAF.
#'
#' @return maf_df with complemented Tumor_Sample_Barcode and other columns ready to be used downstream.
#'
#' @export
#'
#' @examples
#' my_metadata = GAMBLR.data::gambl_metadata
#' reddy_meta = dplyr::filter(my_metadata, cohort=="dlbcl_reddy")
#'
#' small_maf = dplyr::filter(GAMBLR.data::sample_data$grch37$maf,
#'  Tumor_Sample_Barcode %in% reddy_meta$Tumor_Sample_Barcode)
#'
#' small_maf = dplyr::filter(small_maf, Hugo_Symbol == "MYC")
#'
#' complete_maf = supplement_maf(incoming_maf = small_maf,
#'                               these_samples_metadata = reddy_meta)
#'
supplement_maf <- function(incoming_maf,
                           these_samples_metadata){

  missing_sample_ids = setdiff(these_samples_metadata$Tumor_Sample_Barcode,
                               incoming_maf$Tumor_Sample_Barcode)

  missing_sample_maf = incoming_maf %>%
    dplyr::filter(Tumor_Sample_Barcode == "Imaginary Sample ID") %>%
    add_row(Tumor_Sample_Barcode = missing_sample_ids,
           Hugo_Symbol = "GARBAGE",
           Chromosome = ifelse(stringr::str_detect(incoming_maf$Chromosome[1], "chr"), "chr1", "1"),
           Start_Position = 1,
           End_Position = 1,
           Variant_Classification = "Missense_Mutation")

  full_maf = rbind(incoming_maf, missing_sample_maf)
  return(full_maf)
}
