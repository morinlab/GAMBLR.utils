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
#' @param genome_build Specify the genome build (usually not required)
#'
#' @return Copy number matrix with sample_id as rows and regions as columns.
#'
#' @import dplyr circlize tibble stringr tidyr GenomicDistributions
#' @export
#'
#' @examples
#'
#'
#' dlbcl_genome_meta = get_gambl_metadata() %>%
#'     filter(pathology=="DLBCL",
#'     seq_type=="genome")
#' # Create the copy number matrix from this data
#' all_segments = get_cn_segments()
#' all_states_binned = segmented_data_to_cn_matrix(
#'                                   seg_data = all_segments,
#'                                   strategy="auto_split",
#'                                   n_bins_split=2500,
#'                                   missing_data_as_diploid = T,
#'                                   these_samples_metadata = dlbcl_genome_meta)
#' 
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
                            genome_build){
  if(missing(genome_build)){
    if (inherits(seg_data, "seg_data")) {
      genome_build <- get_genome_build.seg_data(seg_data)
    } else {
      stop("specify genome_build or provide a seg_data object")
    }
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
      regions_bed = circlize::read.cytoband(species = "hg19")$df
    }else if(genome_build=="hg38"){
      regions_bed = circlize::read.cytoband(species = "hg38")$df
    }else{
      stop(paste("Unsupported genome_build:",genome_build))
    }
    
    colnames(regions_bed) = c("chromosome_name", "start_position", "end_position", "name", "dunno")
    if(use_cytoband_name){
      regions_bed = mutate(regions_bed, region_name = paste0(str_remove(chromosome_name, pattern = "chr"), name))
      region_names = pull(regions_bed, region_name)
    }else{
      #region_names = pull(regions_bed, region_name)
    }
    regions = apply(regions_bed, 1, bed2region)
    #use the cytobands from the circlize package (currently hg19 but can extend to hg38 once GAMBLR handles it) Has this been updated?
  }else if(strategy=="auto_split"){
    all_len = circlize::read.chromInfo()$chr.len
    bin_df = GenomicDistributions::binChroms(binCount = n_bins_split,chromSizes=all_len[c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")])
    regions = apply(bin_df, 1, bed2region)
  }else if(strategy=="custom_regions"){
    if(!class(regions)=="character"){
      regions = apply(regions_bed, 1, bed2region)
    }
    
  }
  if(!use_cytoband_name & !missing(regions)){
    region_names = names(regions)
  }
  
  if(!missing(these_samples_metadata)){
      #subset to the samples in the provided
      seg_data = dplyr::filter(seg_data,ID %in% these_samples_metadata$sample_id) 
  }
  region_segs = lapply(regions,function(x){get_cn_segments(region = x, streamlined = TRUE, weighted_average = T, seg_data = seg_data)})
  
  tibbled_data = tibble(region_segs, region_name = region_names)
  unnested_df = tibbled_data %>%
    unnest_longer(region_segs)

  seg_df = data.frame(ID = unnested_df$region_segs$ID, CN = unnested_df$region_segs$CN,region_name = unnested_df$region_name)
  
  seg_df = dplyr::rename(seg_df,sample_id=ID) 
  

  if(missing(these_samples_metadata) & !missing(seg_data)){
    eg = expand_grid(sample_id = unique(seg_data$ID), region_name = as.character(unique(seg_df$region_name)))
  }else{
    meta_arranged = these_samples_metadata %>%
      dplyr::select(sample_id, pathology, lymphgen) %>%
      arrange(pathology, lymphgen)
    eg = expand_grid(sample_id = pull(meta_arranged, sample_id), region_name = as.character(unique(seg_df$region_name)))
    
  }
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
