#' @title Build UCSC browser track hub
#' 
#' @description Create a directory that contains the track hub files. They are: 
#' a hub.txt file (marked with `useOneFile on`), and a subdirectory (named as the 
#' projection in use) with the custom tracks to be visualized in UCSC browser. The 
#' custom track `regions` contains the regions from where SSMs were retrieved. All 
#' other custom tracks contain SSMs from samples separated by the `splitColumnName` 
#' parameter, where each file is named according to the value in `splitColumnName` 
#' that it refers to. 
#' 
#' @details `build_browser_hub` will create a custom track file for each combination of 
#' `these_seq_types` and `splitColumnName` (if mutations could be found). Custom 
#' track files are named as `<a_seq_type_value>_<a_splitColumnName_value>.bb/bed`. 
#' 
#' The `bigDataUrl` field of a track in the hub.txt file is defined in 
#' the following way:
#' 
#' ```
#' file.path(bigDataUrl_base, hub_dir, projection, 
#'           paste0(track_file_names[i], "?raw=true\n"))
#' ```
#'
#' where `track_file_names[i]` is the custom track file name.
#'
#' @param regions_bed A BED-format table with the regions you want to retrieve SSMs 
#'   from. The columns 1, 2 and 3 must be chromosome names, start positions and 
#'   end positions, respectively. The default is the `GAMBLR.data::grch37_ashm_regions`
#'   object which contains public-available SSMs from aSHM regions, but any other 
#'   regions can be provided.
#' @param these_sample_ids A vector of sample IDs that you want results for.
#' @param these_samples_metadata A metadata table (with sample IDs in a column) to 
#'   subset the samples of interest.
#' @param these_seq_types A vector of one or more seq types you want results for. 
#'   Possible values are "genome", "capture", or "mrna". If more than 
#'   one seq type is provided, for each value in `splitColumnName`, the function 
#'   creates a separate track for each provided seq type. The default is 
#'   `c("genome", "capture")`. See the **Details** section for more information.
#' @param projection The projection genome build. One of "grch37" (default) or "hg38".
#' @param local_web_host_dir Path to the directory where is your local copy of the
#'   web host to be used. For example, if the hub should be hosted on GitHub, 
#'   `local_web_host_dir` should be the path to your local copy of the repository 
#'   directory. Default is NULL.
#' @param hub_dir Path to the directory (inside a web host) where you want to build 
#'   your track hub. If this directory does not exist, it will be created. 
#' @param as_bigbed A Boolean value. If TRUE (default), custom tracks are saved as
#'   bigBed (.bb) files. BED files otherwise.
#' @param splitColumnName An single string to indicate which metadata column is used 
#'   to split the MAF data into custom track files. Default is "pathology".
#' @param hub_name A string with the hub name (without spaces) to fill in the `hub` 
#'   field of the hub.txt file. Default is `basename(hub_dir)`.
#' @param shortLabel A string with the short hub label (maximum of 17 characters 
#'   recommended; spaces are allowed) to fill in the `shortLabel` field of the 
#'   hub.txt file. Default is `basename(hub_dir)`.
#' @param longLabel A string with the long hub label (maximum of 80 characters 
#'   recommended; spaces are allowed) to fill in the `longLabel` field of the 
#'   hub.txt file. Default is `basename(hub_dir)`.
#' @param email Optional parameter. A string with the contact email to fill in the 
#'   `email` field of the hub.txt file. If not provided, hub.txt will not have this 
#'   field.
#' @param visibility A string that controls the track visibility mode. Possible 
#'   values are "pack", "dense", "full", or "squish" (default).
#' @param bigDataUrl_base A string with the base path of the web location of the 
#'   tracks' data files (`bigDataUrl` fields from the hub.txt file). For example, 
#'   if the hub should be hosted on GitHub, `bigDataUrl_base` should be something 
#'   like "https://github.com/morinlab/LLMPP/blob/main" (the default). See how the 
#'   track's paths from the bigDataUrl fields are set in the **Details** section.
#'
#' @return Nothing.
#' 
#' @import dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' # create a track hub in LLMPP GitHub repo
#' 
#' library(GAMBLR.data)
#' 
#' local_web_host_dir = "~/repos/LLMPP"
#' hub_dir = "hubs/ashm_test"
#' 
#' my_meta = get_gambl_metadata() %>% 
#'   filter(pathology %in% c("BL", "DLBCL", "FL"))
#' 
#' build_browser_hub(
#'   these_samples_metadata = my_meta,
#'   these_seq_types = c("genome", "capture"),
#'   projection = "grch37",
#'   local_web_host_dir = local_web_host_dir,
#'   hub_dir = hub_dir,
#'   splitColumnName = "pathology",
#'   longLabel = "Public aSHM mutations separated by pathologies"
#' )
#' }
#' 
build_browser_hub = function(regions_bed = GAMBLR.data::grch37_ashm_regions,
                             these_sample_ids = NULL,
                             these_samples_metadata = NULL,
                             these_seq_types = c("genome", "capture"),
                             projection = "grch37",
                             local_web_host_dir = NULL,
                             hub_dir = "my_hub",
                             as_bigbed = TRUE,
                             splitColumnName = "pathology",
                             hub_name = basename(hub_dir),
                             shortLabel = basename(hub_dir),
                             longLabel = basename(hub_dir),
                             email,
                             visibility = "squish",
                             bigDataUrl_base = "https://github.com/morinlab/LLMPP/blob/main"){
  
  # check some provided parameter 
  stopifnot("`these_seq_types` must be one or more of \"genome\", \"capture\" or \"mrna\"." = 
              all(these_seq_types %in% c("genome", "capture", "mrna")))
  stopifnot("`projection` must be one of \"grch37\" or \"hg38\"." = 
              projection %in% c("grch37", "hg38"))
  stopifnot("`visibility` must be one of \"pack\", \"dense\", \"full\", or \"squish\"." = 
              visibility %in% c("pack", "dense", "full", "squish"))
  
  # get metadata with the dedicated helper function (for each seq type)
  these_seq_types = setNames(these_seq_types, these_seq_types)
  these_samples_metadata = lapply(these_seq_types, function(these_seq_types_i){
    id_ease(these_samples_metadata = these_samples_metadata,
            these_sample_ids = these_sample_ids,
            verbose = FALSE,
            this_seq_type = these_seq_types_i)
  }) %>% 
    suppressMessages
  stopifnot("`splitColumnName` must be a column name contained in the metadata." = 
              splitColumnName %in% names(these_samples_metadata[[1]]))
  
  # if dir paths contain a forward slash at the end, remove it. 
  bigDataUrl_base = sub("/$", "", bigDataUrl_base)
  local_web_host_dir = sub("/$", "", local_web_host_dir)
  hub_dir = sub("/$", "", hub_dir)
  
  # full path to hub dir
  if(is.null(local_web_host_dir)){
    hub_dir_full_path = hub_dir
  }else{
    hub_dir_full_path = file.path(local_web_host_dir, hub_dir)
  }
  
  # create output directory and sub-directories
  dir.create(hub_dir_full_path, showWarnings = FALSE)
  track_dir = file.path(hub_dir_full_path, projection)
  dir.create(track_dir, showWarnings = FALSE)
  
  # save regions_bed to a bed/bb file
  regions_bed = dplyr::select(regions_bed, 1,2,3) %>% 
    arrange( .[[1]], .[[2]] )
  if(as_bigbed){
    temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
    write.table(regions_bed, temp_bed, quote = FALSE, sep = "\t", row.names = FALSE, 
                col.names = FALSE)
    if(projection == "grch37"){
      chr_arms = GAMBLR.data::chromosome_arms_grch37 %>% 
        mutate(chromosome = paste0("chr", chromosome))
    }else{ # so projection is hg38
      chr_arms = GAMBLR.data::chromosome_arms_hg38
    }
    chr_sizes = chr_arms %>%
      dplyr::filter(arm == "q") %>%
      dplyr::select(chromosome, end) %>%
      rename(size = "end")
    temp_chr_sizes = tempfile(pattern = "chrom.sizes_")
    write.table(chr_sizes, temp_chr_sizes, quote = FALSE, sep = "\t", row.names = FALSE, 
                col.names = FALSE)
    bedtobigbed = GAMBLR.helpers::check_config_value(config::get("dependencies")$bedToBigBed)
    regions_bed_file = file.path(track_dir, "regions.bb")
    bigbed_conversion = gettextf("%s %s %s %s", bedtobigbed, temp_bed, temp_chr_sizes, 
                                 regions_bed_file)
    system(bigbed_conversion)
    unlink(temp_bed)
    unlink(temp_chr_sizes)
  }else{
    regions_bed_file = file.path(track_dir, "regions.bed")
    write.table(regions_bed, regions_bed_file, quote = FALSE, sep = "\t", row.names = FALSE, 
                col.names = FALSE)
  }
  
  # get maf data from the specified regions and metadata/samples (for each seq type)
  maf_data = mapply(function(these_seq_types_i, these_samples_metadata_i){
    get_ssm_by_regions(regions_bed = regions_bed, this_seq_type = these_seq_types_i,
                       these_samples_metadata = these_samples_metadata_i, 
                       projection = projection, streamlined = FALSE, 
                       basic_columns = TRUE) %>% 
      suppressMessages
  }, these_seq_types, these_samples_metadata, SIMPLIFY = FALSE)
  
  # split maf table according to splitColumnName
  if(!is.null(splitColumnName)){
    maf_data = mapply(function(maf_data_i, these_samples_metadata_i){
      maf_data_i = dplyr::select(these_samples_metadata_i, sample_id, all_of(splitColumnName)) %>% 
        distinct(sample_id, .keep_all = TRUE) %>% 
        left_join(maf_data_i, ., join_by(Tumor_Sample_Barcode == sample_id))
      dplyr::select(maf_data_i, -all_of(splitColumnName)) %>% 
        split(maf_data_i[[splitColumnName]])
    }, maf_data, these_samples_metadata)
    track_names = lapply(maf_data, names)
    these_samples_metadata = mapply(function(these_samples_metadata_i, track_names_i){
      split(these_samples_metadata_i, these_samples_metadata_i[[splitColumnName]]) %>% 
        "["(track_names_i)
    }, these_samples_metadata, track_names)
  }else{
    maf_data = lapply(maf_data, list)
    track_names = lapply(these_seq_types, function(x) list("all"))
    these_samples_metadata = lapply(these_samples_metadata, list)
  }
  track_names = mapply(paste, these_seq_types, track_names, sep = "_")
  
  # convert and save track files
  track_file_names = lapply(track_names, paste0, ifelse(as_bigbed, ".bb", ".bed"))
  mapply(
    function(maf_i, meta_i, these_seq_types_i, track_names_file_i){
      mapply(maf_to_custom_track, 
             maf_data = maf_i,
             these_samples_metadata = meta_i,
             this_seq_type = these_seq_types_i,
             output_file = file.path(track_dir, track_names_file_i),
             as_bigbed = as_bigbed,
             projection = projection)
    },
    maf_data, these_samples_metadata, these_seq_types, track_file_names
  ) %>% 
    invisible
  
  
  ### create hub.txt file
  
  # open file
  hub_file = file.path(hub_dir_full_path, "hub.txt")
  sink(hub_file)
  
  # write header
  cat( paste0("hub ", hub_name, "\n") )
  cat( paste0("shortLabel ", shortLabel, "\n") )
  cat( paste0("longLabel ", longLabel, "\n") )
  cat( "useOneFile on\n" )
  if(!missing(email)){
    cat( paste0("email ", email, "\n") )
  }
  cat( "\n" )
  cat( paste0("genome ", replace(projection, projection == "grch37", "hg19"), "\n") )
  
  # write info of the track which stores regions from where ssms are retrieved
  cat( "\n" )
  cat( paste0("track ", hub_name, "_regions\n") )
  cat( paste0("shortLabel ", hub_name, " regions\n") )
  cat( paste0("longLabel Regions where mutations were taken to hub ", hub_name, "\n") )
  cat( paste0("visibility ", visibility, "\n") )
  cat( paste0("priority 1\n") )
  cat( paste0("type ", ifelse(as_bigbed, "bigBed", "bed"), "\n") )
  file.path(bigDataUrl_base, hub_dir, projection, basename(regions_bed_file)) %>% 
    { cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
  
  # write info of the tracks that store ssms split by splitColumnName
  track_names = unlist(track_names)
  track_file_names = unlist(track_file_names)
  for(i in seq_along(track_names)){
    cat( "\n" )
    cat( paste0("track ", hub_name, "_", track_names[i], "\n") )
    cat( paste0("shortLabel ", hub_name, " ", track_names[i], "\n") )
    cat( paste0("longLabel ", hub_name, " ", track_names[i], "\n") )
    cat( paste0("visibility ", visibility, "\n") )
    cat( paste0("priority ", i+1, "\n") )
    cat( paste0("type ", ifelse(as_bigbed, "bigBed", "bed"), "\n") )
    cat( "itemRgb on\n" )
    file.path(bigDataUrl_base, hub_dir, projection, track_file_names[i]) %>% 
      { cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
  }
  
  # close file
  sink()
}
