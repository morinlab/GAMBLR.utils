#' @title Get Sample Wildcards
#'
#' @description Get wildcards for a sample_id/seq_type combination.
#'
#' @details Return sample wildcards, useful for getting wildcard information necessary for retrieving sample-level flat-files with glue.
#'
#' @param this_sample_id The sample ID of interest.
#' @param seq_type The desired seq type, e.g genome/capture.
#'
#' @return Nothing.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' get_sample_wildcards(this_sample_id = "00-15201_tumorA",
#'                      seq_type = "genome")
#'
get_sample_wildcards = function(this_sample_id,
                                seq_type){

  sample_meta = get_gambl_metadata(seq_type_filter = seq_type) %>%
    dplyr::filter(sample_id==this_sample_id)
  this_patient_id = sample_meta$patient_id
  if(sample_meta$pairing_status=="matched"){
    normal_meta = get_gambl_metadata(seq_type_filter = seq_type,tissue_status_filter = c("normal","tumour")) %>%
      dplyr::filter(patient_id==this_patient_id) %>% dplyr::filter(tissue_status=="normal")
      normal_id = normal_meta$sample_id
      return(list(tumour_sample_id=this_sample_id,
               normal_sample_id=normal_id,
               seq_type = seq_type,
               pairing_status=sample_meta$pairing_status,
               genome_build=sample_meta$genome_build,
               unix_group=sample_meta$unix_group))
  }else{
    message("This function only works with matched samples for now")
    return()
  }
}
