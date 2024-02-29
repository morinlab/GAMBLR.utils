#' @title Build UCSC browser track hub
#' 
#' @description 
#' 
#' @details 
#'
#' @param regions_bed A BED-format table with the regions you want to retrive SSMs 
#'   from. The columns 1, 2 and 3 must be chromosome names, start positions and 
#'   end positions, respectively. 
#' @param these_sample_ids A vector of sample_id that you want results for.
#' @param these_samples_metadata A metadata table (with sample IDs in a column) to 
#'   subset the samples of interest.
#' @param this_seq_type The seq type you want results for. Possible values are "genome"
#'   (default), "capture", or "mrna".
#' @param projection The projection genome build. One of "grch37" (default) or "hg38".
#' @param hub_dir Path to the directory where you want to build your track hub. 
#'   If this directory does not exist, it will be created. 
#' @param as_bigbed A Boolean value. If TRUE (default), custom tracks are saved as
#'   bigBed (.bb) files. BED files otherwise.
#' @param splitColumnName An single string to indicate which metadata column is used 
#'   to split the data into custom track files. Default is "pathology".
#' @param hub_name A string with the hub name (without spaces) to fill in the `hub` 
#'   field of the hub.txt file. Defoult is `basename(hub_dir)`.
#' @param shortLabel A string with the short hub label (maximum of 17 characters; 
#'   spaces are allowed) to fill in the `shortLabel` field of the hub.txt file. 
#'   Defoult is `basename(hub_dir)`.
#' @param longLabel A string with the long hub label (maximum of 80 characters; 
#'   spaces are allowed) to fill in the `longLabel` field of the hub.txt file. 
#'   Defoult is `basename(hub_dir)`.
#' @param email A string with the contact email to fill in the `email` field of 
#'   the hub.txt file. Default is `rdmorin@sfu.ca`.
#' @param visibility A string that controls the track visibility mode. Possible 
#'   values are "pack", "dense", "full", or "squish" (default).
#' @param bigDataUrl_dir The web location directory of the track's data file. The 
#'   default is `file.path("https://github.com/morinlab/LLMPP/blob/main/hubs", hub_name)`.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' 
build_browser_hub <- function(regions_bed,
                              these_sample_ids = NULL,
                              these_samples_metadata = NULL,
                              this_seq_type = "genome",
                              projection = "grch37",
                              hub_dir = "./my_hub",
                              as_bigbed = TRUE,
                              splitColumnName = "pathology",
                              hub_name = basename(hub_dir),
                              shortLabel = basename(hub_dir),
                              longLabel = basename(hub_dir),
                              email = "rdmorin@sfu.ca",
                              visibility = "squish",
                              bigDataUrl_dir = file.path("https://github.com/morinlab/LLMPP/blob/main/hubs", hub_name)){
  
  # create output directory and sub-directories
  dir.create(hub_dir, showWarnings = FALSE)
  track_dir <- file.path(hub_dir, projection)
  dir.create(track_dir, showWarnings = FALSE)
  
  # get metadata with the dedicated helper function
  these_samples_metadata = id_ease(these_samples_metadata = these_samples_metadata,
                                   these_sample_ids = these_sample_ids,
                                   verbose = FALSE,
                                   this_seq_type = this_seq_type)
  
  # save regions_bed to a bb file
  temp_bed <- tempfile(pattern = "regionsBed_", fileext = ".bed")
  regions_bed <- select(regions_bed, 1,2,3) %>% 
    arrange( .[[1]], .[[2]] )
  write.table(regions_bed, temp_bed, quote = FALSE, sep = "\t", row.names = FALSE, 
              col.names = FALSE)
  if(projection == "grch37"){
    chr_arms <- GAMBLR.data::chromosome_arms_grch37 %>% 
      mutate(chromosome = paste0("chr", chromosome))
  }else if(projection == "hg38"){
    chr_arms <- GAMBLR.data::chromosome_arms_hg38
  }else{
    stop("projection parameter must be \"grch37\" or \"hg38\".")
  }
  chr_sizes <- chr_arms %>%
    dplyr::filter(arm == "q") %>%
    dplyr::select(chromosome, end) %>%
    rename(size = "end")
  temp_chr_sizes <- tempfile(pattern = "chrom.sizes_")
  write.table(chr_sizes, temp_chr_sizes, quote = FALSE, sep = "\t", row.names = FALSE, 
              col.names = FALSE)
  bedtobigbed <- GAMBLR.helpers::check_config_value(config::get("dependencies")$bedToBigBed)
  bigbed_conversion = file.path(track_dir, "regions.bb") %>% 
    gettextf("%s %s %s %s", bedtobigbed, temp_bed, temp_chr_sizes, .)
  system(bigbed_conversion)
  unlink(temp_bed)
  unlink(temp_chr_sizes)
  
  # get maf data from the specified regions and metadata/samples
  maf_data <- get_ssm_by_regions(regions_bed = regions_bed, this_seq_type = this_seq_type,
                                 these_samples_metadata = these_samples_metadata, 
                                 projection = projection, streamlined = FALSE, 
                                 basic_columns = TRUE) %>% 
    suppressMessages
  
  # split maf table according to splitColumnName
  if(!is.null(splitColumnName)){
    maf_data <- select(these_samples_metadata, sample_id, all_of(splitColumnName)) %>% 
      distinct(sample_id, .keep_all = TRUE) %>% 
      left_join(maf_data, ., join_by(Tumor_Sample_Barcode == sample_id))
    maf_data <- select(maf_data, -all_of(splitColumnName)) %>% 
      split(maf_data[[splitColumnName]])
    track_names = names(maf_data)
    these_samples_metadata = split(these_samples_metadata, these_samples_metadata[[splitColumnName]]) %>% 
      "["(track_names)
  }else{
    maf <- list(maf)
    track_names = "all"
    these_samples_metadata <- list(these_samples_metadata)
  }
  
  # convert and save track files
  track_names_file <- paste0(track_names, ifelse(as_bigbed, ".bb", ".bed"))
  mapply(function(maf_i, track_names_file_i, meta_i){
    maf_to_custom_track(
      maf_i,
      meta_i,
      seq_type = this_seq_type,
      output_file = file.path(track_dir, track_names_file_i),
      as_bigbed = as_bigbed,
      projection = projection
    )
  }, maf_data, track_names_file, these_samples_metadata) %>% 
    invisible
  
  # create hub.txt file
  hub_file <- file.path(hub_dir, "hub.txt")
  
  sink(hub_file)
  
  cat( paste0("hub ", hub_name, "\n") )
  cat( paste0("shortLabel ", shortLabel, "\n") )
  cat( paste0("longLabel ", longLabel, "\n") )
  cat( "useOneFile on\n" )
  cat( paste0("email ", email, "\n") )
  cat( "\n" )
  cat( paste0("genome ", replace(projection, projection == "grch37", "hg19"), "\n") )
  
  for(i in seq_along(track_names)){
    cat( "\n" )
    cat( paste0("track ", hub_name, "_", track_names[i], "\n") )
    cat( paste0("shortLabel ", hub_name, " ", track_names[i], "\n") )
    cat( paste0("longLabel ", hub_name, " ", track_names[i], "\n") )
    cat( paste0("visibility ", visibility, "\n") )
    cat( paste0("priority ", i, "\n") )
    cat( paste0("type ", ifelse(as_bigbed, "bigBed", "bed"), "\n") )
    cat( "itemRgb on\n" )
    sub("/$", "", bigDataUrl_dir) %>% 
      file.path(track_names_file[i]) %>% 
      { cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
  }
  
  sink()
}

