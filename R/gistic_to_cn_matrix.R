



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
#' @param missing_data_as_diploid Fill in gaps as diploid
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
#' all_out = gistic_to_cn_matrix("all_lesions.conf_90.txt",
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
gistic_to_cn_matrix = function(gistic_lesions_file,
                          these_samples_metadata,
                          seg_data,
                          wide_peaks=FALSE,
                          max_CN=6,
                          drop_inconsistent = TRUE,
                          as_binary = TRUE,
                          scale_by_sample=TRUE,
                          missing_data_as_diploid=TRUE,
                          genome_build){
  lesions = suppressMessages(read_tsv(gistic_lesions_file, col_names = TRUE)) %>% 
    filter(!grepl("CN",`Unique Name`))
  
  lesions_regions = select(lesions,1:6)

  if(missing(genome_build)){
    genome_build = get_genome_build(seg_data)
    print(genome_build)
  }
  
  if(wide_peaks){
    peak = select(lesions_regions,1,`Wide Peak Limits`)  %>% 
      mutate(region=str_remove(`Wide Peak Limits`,"\\(.+")) %>%
      mutate(duplicated_region=region) %>%
      separate(duplicated_region,into=c("chrom","start","end")) %>%
      mutate(type=str_extract(`Unique Name`,"(\\S+)")) %>%
      dplyr::rename("Peak Limits"="Wide Peak Limits")
  }else{
    peak = select(lesions_regions,1,`Peak Limits`) %>% 
      mutate(region=str_remove(`Peak Limits`,"\\(.+")) %>%
      mutate(duplicated_region=region) %>%
      separate(duplicated_region,into=c("chrom","start","end")) %>%
      mutate(type=str_extract(`Unique Name`,"(\\S+)"))
  }


  
  
  #10 to ncol(x)-1
  lasti = ncol(lesions)-1
  lesions_values = lesions[,c(1,10:lasti)]
  lesions_values = left_join(lesions_regions,lesions_values) %>%
    mutate(region=str_remove(`Peak Limits`,"\\(.+")) %>%
    select(-1:-6) %>% column_to_rownames("region") %>% t() %>% as.data.frame()
  
  lesions_values_long = rownames_to_column(lesions_values,"sample_id") %>% 
    pivot_longer(-sample_id,names_to="region",values_to="magnitude")
  
  
  if(as_binary){
    lesions_values_long = left_join(lesions_values_long,peak) %>% 
      mutate(CN=case_when(magnitude==0 ~ 0,
                          magnitude==1 & type == "Amplification" ~ 1,
                          magnitude==2 & type == "Amplification" ~ 1,
                          magnitude==1 & type == "Deletion" ~ 1,
                          magnitude==2 & type == "Deletion" ~ 1))
  }else{
    lesions_values_long = left_join(lesions_values_long,peak) %>% 
      mutate(CN=case_when(magnitude==0 ~ 2,
                          magnitude==1 & type == "Amplification" ~ 3,
                          magnitude==2 & type == "Amplification" ~ 4,
                          magnitude==1 & type == "Deletion" ~ 1,
                          magnitude==2 & type == "Deletion" ~ 0))
  }
  lesions_values = pivot_wider(lesions_values_long,names_from="region",id_cols="sample_id",values_from="CN") %>% 
    column_to_rownames("sample_id")

  gh = Heatmap(lesions_values,cluster_columns = T)

  regions_bed = select(peak,chrom,start,end,`Unique Name`) %>%
    arrange(chrom,start)
  if(genome_build=="grch37"){
    regions_bed = mutate(regions_bed,chrom = str_remove(chrom,"chr"))
  }
  if(!missing(these_samples_metadata)){
    gistic_peaks_binned = segmented_data_to_cn_matrix(
                                        regions = regions_bed,
                                        strategy="custom_regions",
                                        missing_data_as_diploid = missing_data_as_diploid,
                                        seg_data = seg_data,
                                        adjust_for_ploidy = scale_by_sample,
                                        these_samples_metadata = these_samples_metadata)
  }else{

    gistic_peaks_binned = segmented_data_to_cn_matrix(
                                        regions = regions_bed,
                                        strategy="custom_regions",
                                        missing_data_as_diploid = missing_data_as_diploid,
                                        seg_data = seg_data,
                                        adjust_for_ploidy = scale_by_sample
    )
  }
  
  
  gistic_peaks_long = gistic_peaks_binned %>% rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id,names_to="region",values_to="CN") %>%
    mutate(original_region=region) %>%
    separate(region,into=c("chrom","startend"),sep=":") %>%
    
    separate(startend,into=c("start","end"),sep="-") %>%
    mutate(CN=round(CN)) %>%
    mutate(CN=ifelse(CN>max_CN,max_CN,CN))

  #re-annotate the rows by type
  
  gistic_peaks_long = left_join(gistic_peaks_long,select(peak,chrom,start,end,type))
  #drop backwards and neutral
  if(drop_inconsistent){
    gistic_peaks_long = mutate(gistic_peaks_long,CN=case_when(type=="Deletion" & CN >1 ~ 2,
                                                              type == "Amplification" & CN < 3 ~ 2,
                                                              TRUE ~ CN)) 
    if(as_binary){
      gistic_peaks_long = mutate(gistic_peaks_long,CN=ifelse(CN==2,0,1))
    }
    
  }

  gistic_peaks_wide = pivot_wider(gistic_peaks_long,
                                  id_cols="sample_id",
                                  names_from="original_region",
                                  values_from="CN") %>% 
    column_to_rownames("sample_id")
 
  hh = Heatmap(gistic_peaks_wide,cluster_columns = T)

  return(list(gambl_cn_matrix=gistic_peaks_wide,
              gistic_cn_matrix=lesions_values,
              gambl_heatmap=hh,
              gistic_heatmap=gh))
  
}

