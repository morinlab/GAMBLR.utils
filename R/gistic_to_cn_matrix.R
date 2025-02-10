



#' Make a binary CN matrix in GISTIC regions using segmented data
#'
#' @param gistic_lesions_file The all_lesions output file from GISTIC from the same pathology you are working on
#' @param these_samples_metadata Optional metadata that will be used to subset your segment data to only overlapping samples 
#' @param seg_data Data frame containing segmented copy number data (i.e. seg format)
#' @param wide_peaks Whether to use wide peaks instead of narrow peaks (FALSE)
#' @param max_CN Maximum value where CN will be truncated to limit the range
#' @param drop_inconsistent Set regions with a CN direction inconsistent with peak type to neutral/diploid
#' @param as_binary One-hot encoding (0 = no CN, 1 = CN)
#' @param scale_by_sample Adjust for overall ploidy of each sample. Default (TRUE) is a close approximation to what GISTIC reports
#' @param fill_missing_with Specify how the value will be assigned to any region not covered by a segment: Options: "diploid" or "avg_ploidy"
#' @param peak_names_from Specify how the columns for the peaks will be named in the result: either "coordinates" or "GISTIC"
#' @param generate_heatmaps Optionally generate heatmaps. Default is TRUE.
#' @param genome_build Specify the genome build if necessary
#'
#' @returns a list of data frames
#' @export
#'
#' @examples
#' 
#' all_segments = get_cn_segments()
#' dlbcl_genomes_meta = get_gambl_metadata() %>% 
#'     dplyr::filter(pathology=="DLBCL",seq_type=="genome")
#'     
#' all_out = gistic_to_cn_state_matrix("all_lesions.conf_90.txt",
#'                               dlbcl_genomes_meta,
#'                               all_segments,
#'                               as_binary = T,
#'                               scale_by_sample = T)
#'                               
#' gambl_cn_matrix_gistic_peaks = all_out$gambl_cn_matrix %>% rownames_to_column("sample_id")
#' 
#' prettyForestPlot(mutmat=gambl_cn_matrix_gistic_peaks, 
#'                         metadata=dlbcl_genomes_meta,
#'                         comparison_column = "COO_consensus",
#'                         comparison_values = c("GCB","ABC"))
#' 
gistic_to_cn_state_matrix <- function(gistic_lesions_file,
                                      these_samples_metadata,
                                      seg_data,
                                      wide_peaks = FALSE,
                                      max_CN = 6,
                                      drop_inconsistent = TRUE,
                                      as_binary = TRUE,
                                      scale_by_sample = TRUE,
                                      peak_names_from = "coordinates",
                                      fill_missing_with = "diploid",
                                      generate_heatmaps = TRUE,
                                      genome_build) {

  process_peaks <- function(lesions_regions, peak_type) {
    
    peaks <- select(lesions_regions, 1, !!sym(peak_limits_column)) %>%
      mutate(sep_region = str_remove(!!sym(peak_limits_column), "\\(.+")) %>%
      mutate(region=sep_region) %>% 
      separate(sep_region, into = c("chrom", "start", "end"), convert = TRUE) %>%
      mutate(type = str_extract(`Unique Name`, "\\S+")) %>%
      mutate(region = str_remove(region, "\\(.+")) %>%
      select(-!!sym(peak_limits_column))
    return(peaks)
  }

  lesions <- suppressMessages(read_tsv(gistic_lesions_file, col_names = TRUE)) %>%
    filter(!grepl("CN", `Unique Name`))
  lesions_regions <- select(lesions, 1:6)

  if (missing(genome_build)) {
    genome_build <- get_genome_build(seg_data)
  }

  
  peak_limits_column <- if (wide_peaks ) "Wide Peak Limits" else "Peak Limits"
  peak <- process_peaks(lesions_regions, peak_limits_column)
  
  lesions_values <- lesions %>%
    select(1, 10:(ncol(lesions) - 1)) %>%
    left_join(lesions_regions,., by = "Unique Name") %>%
    mutate(region = str_remove(!!sym(peak_limits_column), "\\(.+")) %>%
    select(-1:-6)

  lesions_values = lesions_values %>%
    pivot_longer(cols = -region, names_to = "sample_id", values_to = "magnitude")

  
  if (as_binary) {
    lesions_values <- lesions_values %>%
      left_join(peak, by = "region")  %>%
      mutate(CN = case_when(
        magnitude == 0 ~ 0,
        magnitude %in% c(1, 2) & type == "Amplification" ~ 1,
        magnitude %in% c(1, 2) & type == "Deletion" ~ 1
      ))
  } else {
    lesions_values <- lesions_values %>%
      left_join(peak, by = "region") %>%
      mutate(CN = case_when(
        magnitude == 0 ~ 2,
        magnitude == 1 & type == "Amplification" ~ 3,
        magnitude == 2 & type == "Amplification" ~ 4,
        magnitude == 1 & type == "Deletion" ~ 1,
        magnitude == 2 & type == "Deletion" ~ 0
      ))
  }

  lesions_values = lesions_values %>% select(-magnitude,-`Unique Name`:-type)


  lesions_values <- lesions_values %>%
    pivot_wider(id_cols = "sample_id",names_from = "region", values_from = "CN") %>%
   column_to_rownames("sample_id")
  if(generate_heatmaps){
    gh <- Heatmap(lesions_values, cluster_columns = TRUE)
  }
  

  regions_bed <- select(peak, chrom, start, end, `Unique Name`) %>%
    arrange(chrom, start)
  if(peak_names_from == "coordinates"){
    regions_bed = mutate(regions_bed,name = paste0(chrom,":",start,"-",end)) %>%
    select(-`Unique Name`)

  }

  regions_bed = create_bed_data(regions_bed,genome_build = genome_build)
  
  if (!missing(these_samples_metadata)) {
    gistic_peaks_binned <- segmented_data_to_cn_matrix(
      regions = regions_bed,
      strategy = "custom_regions",
      fill_missing_with = fill_missing_with,
      seg_data = seg_data,
      these_samples_metadata = these_samples_metadata
    )
  } else {
    gistic_peaks_binned <- segmented_data_to_cn_matrix(
      regions = regions_bed,
      strategy = "custom_regions",
      fill_missing_with = fill_missing_with,
      seg_data = seg_data,
      adjust_for_ploidy = scale_by_sample
    )
  }

  gistic_peaks_long <- gistic_peaks_binned %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(cols = -sample_id, names_to = "region", values_to = "CN") %>%
    separate(region, into = c("chrom", "startend"), sep = ":") %>%
    separate(startend, into = c("start", "end"), sep = "-", convert = TRUE) %>%
    mutate(CN = pmin(round(CN), max_CN))
  
  if(peak_names_from=="coordinates"){
    peak = mutate(peak,name=region)
  }else{
    peak = dplyr::rename(peak,name=`Unique Name`)
  }

  gistic_peaks_long <- left_join(gistic_peaks_long, select(peak, chrom, start, end, type, name), by = c("chrom", "start", "end"))

  if (drop_inconsistent) {
    gistic_peaks_long <- mutate(gistic_peaks_long, CN = case_when(
      type == "Deletion" & CN > 1 ~ 2,
      type == "Amplification" & CN < 3 ~ 2,
      TRUE ~ CN
    ))
    if (as_binary) {
      gistic_peaks_long <- mutate(gistic_peaks_long, CN = ifelse(CN == 2, 0, 1))
    }
  }

  gistic_peaks_wide <- gistic_peaks_long %>%
    select(-chrom,-start,-end,-type) %>% 
    pivot_wider(id_cols="sample_id",names_from = "name", values_from = "CN") %>%
    column_to_rownames("sample_id")

  if (peak_names_from == "coordinates") {
    colnames(gistic_peaks_wide) <- colnames(lesions_values)
  } else {
    colnames(lesions_values) <- colnames(gistic_peaks_wide)
  }
  if(generate_heatmaps){
    hh <- Heatmap(gistic_peaks_wide, cluster_columns = TRUE)
    return(list(
      gambl_cn_matrix = gistic_peaks_wide,
      gistic_cn_matrix = lesions_values,
      gambl_heatmap = hh,
      gistic_heatmap = gh
    ))
  }else{
    return(list(
      gambl_cn_matrix = gistic_peaks_wide,
      gistic_cn_matrix = lesions_values
    ))
  }
}