
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
#' @param fill_missing_with Specify how the value will be assigned to any region not covered by a segment: Options: "diploid" or "avg_ploidy"
#'
#' @return data.frame
#'
#'
#' @import GAMBLR.helpers
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
                                         filler_values,
                                         verbose=FALSE){


  region = gsub(",", "", region)
  split_chunks = unlist(strsplit(region, ":"))
  chromosome = split_chunks[1]
  startend = unlist(strsplit(split_chunks[2], "-"))
  qstart = as.numeric(startend[1])
  qend = as.numeric(startend[2])

  if(verbose){
    print(region)
  }
  

  all_segs =
    dplyr::filter(seg_data, (chrom == chromosome & start >= qstart & start <= qend)|
                    (chrom == chromosome & end > qstart & end <= qend)|
                    (chrom == chromosome & end > qend & start <= qstart)) %>%
    mutate(start=ifelse(start < qstart,qstart,start),end=ifelse(end>qend,qend,end)) %>%
    mutate(length=end-start) %>% 
    dplyr::filter(length>0)
  if(nrow(all_segs)==0){
    if(verbose){
      print(paste("NO values found for",region))
    }
    #there is no data from any patients within this region!
    dummy_full = mutate(filler_values,
      chrom=chromosome,start=qstart,end=qend)
    return(dummy_full)
  }
  all_segs = all_segs %>% mutate(CN_L = length * CN,logr_L = length*log.ratio)
  if(weighted_average){
    all_segs = all_segs %>%
      group_by(ID) %>%
      summarise(total_L = sum(length), log.ratio = sum(logr_L)/sum(length),
                CN = sum(CN_L)/sum(length)) %>% ungroup()
    if(verbose){
      print(head(all_segs))
    }

  }else{
    all_segs = dplyr::mutate(all_segs, CN = 2*2^log.ratio)
  }
  if(!missing(filler_values)){
    filler_values$chrom = chromosome
    filler_values$start = qstart
    filler_values$end = qend
    filler_values = dplyr::filter(filler_values,!ID %in% all_segs$ID)
    all_segs = bind_rows(filler_values,all_segs)
  }
  



  #check for NaN values
  n_sample_nan = nrow(all_segs %>% filter(is.nan(CN)))
  if (n_sample_nan>0) {
    stop("NAN values found")
    
  }

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
#' @param fill_missing_with Fill in any sample/region combinations with missing data as diploid or the average ploidy ("diploid" or "avg_ploidy")
#' @param adjust_for_ploidy Whether to adjust for high ploidy
#' @param max_cn_allowed CN values higher than this will be set to this value. Default 6.
#' @param gistic_lesions_file Path to gistic lesions file (only needed if strategy="GISTIC")
#' @param rounded Set to FALSE if you want the raw averaged copy number values per region
#' @param verbose Set to TRUE for more messages
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
                            fill_missing_with = "diploid",
                            max_CN_allowed = 6,
                            adjust_for_ploidy=FALSE,
                            genome_build,
                            gistic_lesions_file,
                            rounded = TRUE,
                            verbose = FALSE){
  if(!is.numeric(max_CN_allowed)){
    stop("max_CN_allowed must be a numeric value")
  }
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
  if (fill_missing_with=="diploid") {
    
    
    dummy_df = data.frame(ID=unique(seg_data$ID),CN=2,log.ratio=0)
    if(verbose){
        print("Fill with diploid")
        print(head(dummy_df))
    }
  } else if  (fill_missing_with=="avg_ploidy") {
    if(!"dummy_segment" %in% colnames(seg_data)){
      stop("avg_ploidy mode only available when there's a dummy_segment column")
    }
    dummy_df = seg_data %>% 
      mutate(
          length = end - start + 1,
          CN_seg = CN * length,
          logr_seg = log.ratio * length
        ) %>%
        group_by(ID) %>%
        summarise(
          CN = sum(CN_seg) / sum(length),
          log.ratio = sum(logr_seg) / sum(length)
        ) # actual average per base
    print(head(dummy_df))
  }else if(fill_missing_with == "nothing"){
    #do nothing
    print("NOT filling gaps")
  }else{
    stop("fill_missing_with must be 'nothing', 'diploid' or 'avg_ploidy'")
  }
  if(strategy=="GISTIC"){
    # extract the regions from the GISTIC file and ignore the data from the samples GISTIC was run on
    # The matrix returned is 0/1 based on the direction (gain peak or loss peak) so cannot be used directly 
    # as a CN state matrix in this application
    gistic_processed = gistic_to_cn_state_matrix(seg_data=seg_data,
                                    gistic_lesions_file=gistic_lesions_file,
                                    wide_peaks=TRUE,
                                    drop_inconsistent=TRUE,
                                    scale_by_sample=adjust_for_ploidy,
                                    fill_missing_with=fill_missing_with,
                                    these_samples_metadata = these_samples_metadata,
                                    peak_names_from = "coordinates",
                                    generate_heatmaps = FALSE)
    #fill values for these regions using the data provided to the function

    peak_regions = colnames(gistic_processed$gistic_cn_matrix)

    region_processed = process_regions(regions_list=peak_regions,
                               projection = genome_build,sort=T)

    regions_bed=region_processed$regions_bed
    filled = segmented_data_to_cn_matrix(seg_data = seg_data,
                                         strategy = "custom_regions",
                                         regions = regions_bed,
                                         fill_missing_with = fill_missing_with,
                                         adjust_for_ploidy = adjust_for_ploidy,
                                         these = these_samples_metadata,
                                         genome_build = genome_build,
                                         verbose=verbose)
    return(filled)
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
    seg_data = mutate(seg_data,
      size=end-start,
      CN=ifelse(CN>max_CN_allowed,
              max_CN_allowed,CN),
              weighted = CN*size) %>%
      group_by(ID) %>% 
      mutate(average=sum(weighted)/sum(size)) %>%
      mutate(CN=2+CN-average)
  }

  #retrieve the CN value for this region for every segment that overlaps it
  bed2region=function(x){
    paste0(x[1], ":", as.integer(x[2]), "-", as.integer(x[3]))
  }
  if(all_cytobands){
    message(paste("Working with cytobands for",
                 genome_build,
                 "This will take awhile but it does work, trust me!"))
    if(genome_build=="grch37"){
      regions_bed = GAMBLR.data::cytobands_grch37
    }else if(genome_build=="hg38"){
      regions_bed = GAMBLR.data::cytobands_hg38
    }else{
      stop(paste("Unsupported genome_build:",genome_build))
    }

    colnames(regions_bed) = c("chromosome_name", "start_position", "end_position", "name", "dunno")
    if(use_cytoband_name){
      regions_bed = mutate(regions_bed,
                           region_name = paste0(str_remove(chromosome_name,
                                                           pattern = "chr"),
                                                           name))
      region_names = pull(regions_bed, region_name)
    }else{
      regions_bed = mutate(regions_bed,
                          region_name = paste0(chromosome_name,":",
                                               as.integer(start_position),
                                               "-",as.integer(end_position)))
      region_names = pull(regions_bed, region_name)
    }
    regions = apply(regions_bed, 1, bed2region)
    
    names(regions) = region_names
    if(verbose){
      print(head(regions))

    }
    #use the cytobands from the circlize package (currently hg19 but can extend to hg38 once GAMBLR handles it) Has this been updated?
  }else if(strategy=="auto_split"){
    if(genome_build=="grch37"){
      length_df = GAMBLR.data::cytobands_grch37 %>%
        mutate(cb.chromosome = factor(cb.chromosome,
                                     levels=unique(cb.chromosome))) %>%
        group_by(cb.chromosome) %>% slice_max(cb.end)
      all_lengths = pull(length_df,cb.end)
      names(all_lengths) = pull(length_df,cb.chromosome)

      bin_df = GenomicDistributions::binChroms(binCount = n_bins_split,
                                               chromSizes=all_lengths)

    }else if(genome_build=="hg38"){
      length_df = GAMBLR.data::cytobands_hg38 %>%
        mutate(cb.chromosome = factor(cb.chromosome,
                                      levels=unique(cb.chromosome))) %>%
        group_by(cb.chromosome) %>% slice_max(cb.end)
      all_lengths = pull(length_df,cb.end)
      names(all_lengths) = pull(length_df,cb.chromosome)
      bin_df = GenomicDistributions::binChroms(binCount = n_bins_split,
                                               chromSizes=all_lengths)

    }else{
      stop(paste("genome build not supported or not specified:",genome_build))
    }
    if(verbose){
      print(paste("created",nrow(bin_df),"bins"))
    }
    regions = apply(bin_df, 1, bed2region)
    region_names = regions

  }else if(strategy=="custom_regions"){
    if("data.frame" %in% class(regions)){
      region_names = unname(unlist(regions[,4]))
      regions = apply(regions, 1, bed2region)
      names(regions) = str_replace_all(region_names,"\\s","")
    }else{

      region_names = names(regions)
      if(is.null(region_names)){
        region_names = regions
      }

    }

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
  if(fill_missing_with=="nothing"){
    region_segs = lapply(regions,
                       function(x){
                        process_cn_segments_by_region(region = x,
                                                      streamlined = TRUE,
                                                      weighted_average = T,
                                                      seg_data = seg_data,
                                                      verbose=verbose)})
  }else{
    if(verbose){
      print("running process_cn_segments_by_region")
    }
    region_segs = lapply(regions,
                       function(x){
                        process_cn_segments_by_region(region = x,
                                                      streamlined = TRUE,
                                                      weighted_average = T,
                                                      seg_data = seg_data,
                                                      filler_values = dummy_df,
                                                      verbose=verbose)})
  }
  


  tibbled_data = tibble(region_segs, region_name = region_names)
  unnested_df = tibbled_data %>%
    unnest_longer(region_segs)

  seg_df = data.frame(ID = unnested_df$region_segs$ID,
                      CN = unnested_df$region_segs$CN,
                      region_name = unnested_df$region_name)

  seg_df = dplyr::rename(seg_df,sample_id=ID)


  eg = expand_grid(sample_id = unique(seg_data$ID),
                   region_name = as.character(unique(seg_df$region_name)))

  all_cn = left_join(eg,
                     seg_df,
                     by = c("sample_id" = "sample_id",
                            "region_name" = "region_name"))
  n_na = dplyr::filter(all_cn,is.na(CN)) %>% nrow()
  print(paste(n_na,"rows with NA values"))
  #fill in any sample/region combinations with missing data as diploid
  #if(fill_missing_with=="diploid"){
  #  all_cn = mutate(all_cn, CN = replace_na(CN, 2))
  #}else if(fill_missing_with == "nothing"){
  #  #leave NA values
  #  print("matrix will include NA values!")
  #}else if(fill_missing_with == "avg_ploidy"){
  #  #Helper function to fill in NA values with the average value for the same sample
  #  fill_na_with_cn <- function(large_df, small_df) {
  #    # Check that both data frames have the same row names
  #    if (!identical(rownames(large_df), rownames(small_df))) {
  #      not_in_dummy = rownames(large_df)[!rownames(large_df) %in% rownames(small_df)]
  #      nnot = length(not_in_dummy)
  #      print(paste(nnot,head(not_in_dummy)))
  #      stop("The row names of the two data frames do not match.")
  #    }
  #    # Loop over each row
  #    for (r in rownames(large_df)) {
  #      # Find which elements in the current row are NA
  #      na_idx <- is.na(large_df[r, ])
  #  
  #      # If any NA exists, replace them with the CN value from small_df
  #      if (any(na_idx)) {
  #        large_df[r, na_idx] <- small_df[r, "CN"]
  #      }
  #    }  
  #  return(large_df)
  #  }
  #  print(head(all_cn))
  #  x = dplyr::filter(all_cn,is.na(CN)) %>% head()
  #  print(x)
  #  all_cn = fill_na_with_cn(all_cn,column_to_rownames(dummy_df,"ID"))
  #}

  cn_matrix = pivot_wider(all_cn,
                          id_cols = "sample_id",
                          names_from = "region_name",
                          values_from = "CN") %>%
    column_to_rownames("sample_id")
  if(verbose){
    print(paste("regions:",length(region_names)))
  }
  #order the regions the same way the user provided them for convenience
  if(any(!region_names %in% colnames(cn_matrix))){
    missing = region_names[!region_names %in% colnames(cn_matrix)]
    nmissing = length(missing)
    region_names = region_names[!region_names %in% missing]
    message(paste("missing data in all samples for", nmissing, "regions"))
  }
  cn_matrix = cn_matrix[,region_names, drop=FALSE]

  cn_matrix[cn_matrix>max_CN_allowed]=max_CN_allowed
  if(rounded){
    cn_matrix = round(cn_matrix)
  }
  return(cn_matrix)
}
