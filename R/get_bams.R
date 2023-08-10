#' @title Get BAMs.
#'
#' @description Get full paths for bam files for a sample or patient.
#'
#' @details Returns a list with BAM paths for tumour, normal and mrna data.
#' This function expects a sample ID (`this_sample_id`) or a patient ID (`this_patient_id`).
#'
#' @param this_sample_id Sample ID of interest.
#' @param this_patient_id patient ID of interest.
#'
#' @return A list that contains the genome_build and an igv-friendly build (igv_build), a list of bam file paths for tumour, normal and mrna data.
#'
#' @import dplyr
#' @export
#'
#' @examples
#'
#' #example 1, using a sample ID
#' bam_details = get_bams(this_sample_id = "HTMCP-01-06-00422-01A-01D")
#'
#' #example 2, using a patient ID
#' bam_details = get_bams(this_patient_id = "HTMCP-01-06-00422")
#'
get_bams = function(this_sample_id,
                    this_patient_id){

  meta = get_gambl_metadata(tissue_status_filter = c("tumour", "normal"), seq_type_filter = seq_type_filter)
  meta_tumour = get_gambl_metadata(tissue_status_filter = c("tumour"), seq_type_filter = seq_type_filter)
  meta_mrna = get_gambl_metadata(seq_type_filter = "mrna")
  #get all samples for this patient
  if(missing(this_patient_id)){
    this_patient_id = meta %>%
      dplyr::filter(sample_id == this_sample_id) %>%
      dplyr::pull(patient_id)
  }
  meta_patient = meta %>%

    dplyr::filter(patient_id == this_patient_id)

  meta_mrna_patient = meta_mrna %>%
    dplyr::filter(patient_id == this_patient_id)

  build = dplyr::pull(meta_patient, genome_build) %>%
    head(1)
  if(build == "hs37d5"){
    igv_build = "hg19"
  }else{
    igv_build = build
  }
  bam_path_pattern = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/{seq_type}_bams/{sample_id}.{genome_build}.bam"

  #tumour_genome_bams = dplyr::filter(meta_patient, seq_type == seq_type_filter & tissue_status == "tumour") %>%
  #  dplyr::pull(data_path)
  tumour_genome_bams = mutate(meta_patient,bam_path=glue::glue(bam_path_pattern)) %>% pull(bam_path)
  
  bam_details = list(igv_build = igv_build, genome_build = build, tumour_bams = tumour_genome_bams)
  normal_genome_bams = dplyr::filter(meta_patient, seq_type == seq_type_filter & tissue_status == "normal") %>%
    dplyr::pull(data_path)

  unix_group = dplyr::filter(meta_patient, seq_type == seq_type_filter & tissue_status == "tumour") %>% slice_head() %>%
    dplyr::pull(unix_group) 
  


  bam_details$pairing_status = get_gambl_metadata(seq_type_filter = "genome") %>%
    dplyr::filter(tissue_status == "tumour", patient_id == this_patient_id) %>%
    dplyr::pull(pairing_status) %>%
    unique()


  bam_details$unix_group = unix_group
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }else{
    print("No Normal")
  }
  
  rnaseq_bams = dplyr::filter(meta_mrna_patient, seq_type == "mrna") %>%
    dplyr::pull(data_path)
  if(length(rnaseq_bams)){
    bam_details$rnaseq_bams = rnaseq_bams
  }
  return(bam_details)
}
