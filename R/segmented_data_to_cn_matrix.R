
#' Filter and summarise CN values seg_data for a single genomic region
#'
#' @description Assign the copy number value to a region based on the segment(s) that overlap it
#'
#' @details This function returns CN states for a single region specified in the format "chromosome:start-end"
#' using segmented copy number data (seg_data). It will either return a segment in the standard format (all original columns)
#' or a streamlined format with only the ID and CN columns
#'
#'
#' @param seg_data Data that will be processed for the specified region
#' @param region Region to be processed in the format "chromosome:start-end"
#' @param streamlined If TRUE, only return the ID and CN columns (default is FALSE)
#' @param weighted_average If TRUE, calculate the CN value as a weighted average of the segments
#' that overlap the region. Otherwise, use the CN value of the first segment that overlaps the region (default is TRUE)
#' @param missing_data_as_diploid If there is no data for the region assume it is diploid
#' 
#' @return data.frame
#' @export
#'
#' @examples
#' 
#' region_segs = process_cn_segments_by_region(
#'        region = "chrX:1-1000000",
#'        streamlined = FALSE, 
#'        weighted_average = T, 
#'        seg_data = seg_data)
#' 
process_cn_segments_by_region = function(seg_data,
                                         region,
                                         streamlined=FALSE,
                                         weighted_average=TRUE,
                                         missing_data_as_diploid=TRUE){
  region = gsub(",", "", region)
  split_chunks = unlist(strsplit(region, ":"))
  chromosome = split_chunks[1]
  startend = unlist(strsplit(split_chunks[2], "-"))
  qstart = as.numeric(startend[1])
  qend = as.numeric(startend[2])
  
  #Save for later 
  dummy_df = data.frame(ID=unique(seg_data$ID),CN=2,log.ratio=0)
  
  all_segs = 
    dplyr::filter(seg_data, (chrom == chromosome & start >= qstart & start <= qend)|
                    (chrom == chromosome & end > qstart & end <= qend)|
                    (chrom == chromosome & end > qend & start <= qstart)) %>%
    mutate(start=ifelse(start < qstart,qstart,start),end=ifelse(end>qend,qend,end)) %>% 
    mutate(length=end-start)
  
  all_segs = all_segs %>% mutate(CN_L = length * CN,logr_L = length*log.ratio) 
  if(weighted_average){
    all_segs = all_segs %>% 
      group_by(ID) %>%
      summarise(total_L = sum(length), log.ratio = sum(logr_L)/sum(length), 
                CN = sum(CN_L)/sum(length)) %>% ungroup() 
    
  }else{
    all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))
  }
  if(missing_data_as_diploid){
     dummy_df = dplyr::filter(dummy_df,!ID %in% all_segs$ID)
     all_segs = bind_rows(dummy_df,all_segs)
  }
  #clean up NaN values
  all_segs = mutate(all_segs,
                    CN=ifelse(is.nan(CN),2,CN),
                    log.ratio=ifelse(is.nan(log.ratio),0,log.ratio))
  if(!streamlined){
    all_segs = mutate(all_segs,
                      chrom=chromosome,
                      start=qstart,
                      end=qend,
                      LOH_flag=NA,
                      )
    all_segs = dplyr::select(all_segs, ID, chrom, start, end, LOH_flag, log.ratio, CN)

  }else{
    all_segs = dplyr::select(all_segs, ID, CN)
  }
  
  return(all_segs)
  
}


#' @title Segmented data to CN matrix
#'
#' @description Convert segmented data to a matrix of CN States for a set of regions.
#'
#' @details This function returns CN states for the specified regions using the CN data in seg_data and (optionally) assumes regions with no data are diploid.
#' For how to determine/specify the coordinates of each region, refer to the parameter descriptions and examples.
#'
#' @param seg_data A data frame of segments that will be used to infer the copy number state of each region
#' @param strategy The general strategy to define regions. Available options are: 'custom_regions','auto_split','cytobands','GISTIC' 
#' @param regions Required when strategy is set to 'custom_regions'. A data frame in bed-like format or a vector of regions in the format "chrom:start-end"
#' @param these_samples_metadata Optional metadata table to auto-subset the data to samples in that table before returning. If missing, the result will include a row for every sample in seg_data.
#' @param n_bins_split Split genome into N equally sized bins
#' @param use_cytoband_name Use cytoband names instead of region names, e.g p36.33.
#' @param missing_data_as_diploid Fill in any sample/region combinations with missing data as diploid (e.g., CN state like 2). Default is FALSE.
#' @param adjust_for_ploidy Whether to adjust for high ploidy
#' @param gistic_lesions_file Path to gistic lesions file (only needed if strategy="GISTIC")
#' @param genome_build Specify the genome build (usually not required)
#'
#' @return Copy number matrix with sample_id as rows and regions as columns.
#'
#' @import dplyr tibble stringr tidyr GenomicDistributions
#' @export
#'
#' @examples
#'
#'
#' dlbcl_genome_meta = get_gambl_metadata() %>%
#'     filter(pathology=="DLBCL",
#'     seq_type=="genome")
#' # Create the copy number matrix from this data
#' all_segments = get_cn_segments(these_samples_metadata = dlbcl_genome_meta)
#' all_states_binned = segmented_data_to_cn_matrix(
#'                                   seg_data = all_segments,
#'                                   strategy="auto_split",
#'                                   n_bins_split=2500,
#'                                   missing_data_as_diploid = T,
#'                                   these_samples_metadata = dlbcl_genome_meta)
#' 
#'
#' gistic_cn_mat = segmented_data_to_cn_matrix(
#'                              these_samples_metadata = dlbcl_genome_meta,
#'                              seg_data=all_segments,
#'                              strategy="GISTIC",
#'                              gistic_lesions_file="all_lesions.conf_90.txt")
#'
#'
segmented_data_to_cn_matrix = function(seg_data,
                            strategy="auto_split",
                            regions,
                            these_samples_metadata,
                            n_bins_split=1000,
                            use_cytoband_name = FALSE,
                            missing_data_as_diploid = FALSE,
                            adjust_for_ploidy=FALSE,
                            genome_build,
                            gistic_lesions_file){
  if(missing(these_samples_metadata)){
    print("missing these_samples_metadata")
  }
  if(missing(genome_build)){
    if (inherits(seg_data, "seg_data")) {
      genome_build <- get_genome_build(seg_data)
    } else {
      stop("specify genome_build or provide a seg_data object")
    }
  }
  if(strategy=="GISTIC"){
    gistic_processed = gistic_to_cn_matrix(seg_data=seg_data,
                                    gistic_lesions_file=gistic_lesions_file,
                                    wide_peaks=TRUE,
                                    drop_inconsistent=TRUE,
                                    scale_by_sample=adjust_for_ploidy,
                                    missing_data_as_diploid=missing_data_as_diploid)
    return(gistic_processed$gambl_cn_matrix)
  }
  if(any(missing(seg_data))){
    stop("one or more required arguments are missing. Required: seg_data, genome_build")
  }
  if(strategy=="cytobands"){
    all_cytobands = TRUE
  }else{
    all_cytobands = FALSE
  }
  
  if(adjust_for_ploidy){
    seg_data = mutate(seg_data,size=end-start,CN=ifelse(CN>5,5,CN),weighted = CN*size) %>% 
      group_by(ID) %>% mutate(average=sum(weighted)/sum(size)) %>% 
      mutate(CN=round(CN-average+2))
  }

  #retrieve the CN value for this region for every segment that overlaps it
  bed2region=function(x){
    paste0(x[1], ":", as.integer(x[2]), "-", as.integer(x[3]))
  }
  if(all_cytobands){
    message(paste("Working with cytobands for", genome_build,"This will take awhile but it does work, trust me!"))
    if(genome_build=="grch37"){
      regions_bed = GAMBLR.data::cytobands_grch37
    }else if(genome_build=="hg38"){
      regions_bed = GAMBLR.data::cytobands_hg38
    }else{
      stop(paste("Unsupported genome_build:",genome_build))
    }
    
    colnames(regions_bed) = c("chromosome_name", "start_position", "end_position", "name", "dunno")
    if(use_cytoband_name){
      regions_bed = mutate(regions_bed, region_name = paste0(str_remove(chromosome_name, pattern = "chr"), name))
      region_names = pull(regions_bed, region_name)
    }else{
      regions_bed = mutate(regions_bed, region_name = paste0(chromosome_name,":",as.integer(start_position),"-",as.integer(end_position)))
      region_names = pull(regions_bed, region_name)
    }
    regions = apply(regions_bed, 1, bed2region)
    #region_names = regions
    #use the cytobands from the circlize package (currently hg19 but can extend to hg38 once GAMBLR handles it) Has this been updated?
  }else if(strategy=="auto_split"){
    if(genome_build=="grch37"){
      length_df = GAMBLR.data::cytobands_grch37 %>% group_by(cb.chromosome) %>% slice_max(cb.end) 
      all_lengths = pull(length_df,cb.end)
      names(all_lengths) = pull(length_df,cb.chromosome)
      
      bin_df = GenomicDistributions::binChroms(binCount = n_bins_split,chromSizes=all_lengths)
      
    }else if(genome_build=="hg38"){
      length_df = GAMBLR.data::cytobands_hg38 %>% group_by(cb.chromosome) %>% slice_max(cb.end) 
      all_lengths = pull(length_df,cb.end)
      names(all_lengths) = pull(length_df,cb.chromosome)
      bin_df = GenomicDistributions::binChroms(binCount = n_bins_split,chromSizes=all_lengths)
      
    }else{
      stop(paste("genome build not supported or not specified:",genome_build))
    }
    
    regions = apply(bin_df, 1, bed2region)
    region_names = regions

  }else if(strategy=="custom_regions"){
    if("data.frame" %in% class(regions)){
      region_names = unname(unlist(regions[,4]))
      #print(class(region_names))
      regions = apply(regions, 1, bed2region)
      #print(regions)
      #print(region_names)
      names(regions) = str_replace_all(region_names,"\\s","")
    }
    #print("HERE")
    #print(regions)
    #print("====")
    
  }
  if(strategy != "cytobands"){
    if(!is.null(names(regions))){
      region_names = names(regions)
    }
  }
  if(!missing(these_samples_metadata)){
      #subset to the samples in the provided
      seg_data = dplyr::filter(seg_data,ID %in% these_samples_metadata$sample_id) 
  }

  region_segs = lapply(regions,function(x){process_cn_segments_by_region(region = x, streamlined = TRUE, weighted_average = T, seg_data = seg_data)})

  
  tibbled_data = tibble(region_segs, region_name = region_names)
  unnested_df = tibbled_data %>%
    unnest_longer(region_segs)

  seg_df = data.frame(ID = unnested_df$region_segs$ID, CN = unnested_df$region_segs$CN,region_name = unnested_df$region_name)
  
  seg_df = dplyr::rename(seg_df,sample_id=ID) 


  eg = expand_grid(sample_id = unique(seg_data$ID), region_name = as.character(unique(seg_df$region_name)))

  all_cn = left_join(eg, seg_df, by = c("sample_id" = "sample_id", "region_name" = "region_name"))

  #fill in any sample/region combinations with missing data as diploid
  if(missing_data_as_diploid){
    all_cn = mutate(all_cn, CN = replace_na(CN, 2))
  }

  cn_matrix = pivot_wider(all_cn, id_cols = "sample_id", names_from = "region_name", values_from = "CN") %>%
    column_to_rownames("sample_id")

  #order the regions the same way the user provided them for convenience
  if(any(!region_names %in% colnames(cn_matrix))){
    missing = region_names[!region_names %in% colnames(cn_matrix)]
    nmissing = length(missing)
    region_names = region_names[!region_names %in% missing]
    message(paste("missing data for",nmissing,"regions"))
  }
  cn_matrix = cn_matrix[,region_names, drop=FALSE]

  return(cn_matrix)
}
