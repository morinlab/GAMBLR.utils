#' @title Build UCSC browser track hub
#' 
#' @description Create a directory that contains the track hub files. They are: 
#' a hub.txt file (named as `<projection>_hub.txt` and marked with `useOneFile on`), 
#' and a subdirectory (named as the projection in use) with the custom tracks to 
#' be visualized in UCSC browser. The custom track `regions` contains the regions 
#' from where SSMs were retrieved. All other custom tracks contain SSMs from samples 
#' separated by the `splitColumnName` parameter, where each file is named according 
#' to the value in `splitColumnName` that it refers to. 
#' 
#' @details `build_browser_hub` will create a custom track file for each combination of 
#' `these_seq_types` and `splitColumnName` (if mutations can be found). Custom 
#' track files are named as `<a_seq_type_value>_<a_splitColumnName_value>.bb`. 
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
#' @param maf_data A MAF-format data frame with the SSMs to be visualized in the 
#'   UCSC browser. It can be filtered by regions and samples provided by parameters 
#'   `regions_bed` and `these_sample_ids`/`these_samples_metadata`, respectively.
#' @param regions_bed A BED-format table with the regions you want to retrieve SSMs 
#'   from. The columns 1, 2 and 3 must be chromosome names, start positions and 
#'   end positions, respectively. If not provided, the function defaults to aSHM 
#'   regions, *i.e.* one of `GAMBLR.data::grch37_ashm_regions` or `GAMBLR.data::hg38_ashm_regions`,
#'   which is chosen automatically depending on the provided `projection`.
#' @param these_sample_ids A vector of sample IDs that you want results for.
#' @param these_samples_metadata A metadata table (with sample IDs in a column) to 
#'   subset the samples of interest.
#' @param these_seq_types A vector of one or more seq types you want results for. 
#'   Possible values are "genome" or "capture". If more than 
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
#' @param contact_email Required parameter. A string with the contact email to fill 
#'   in the `email` field of the hub.txt file.
#' @param visibility A string that controls the track visibility mode. Possible 
#'   values are "pack", "dense", "full", or "squish" (default).
#' @param bigDataUrl_base A string with the base path of the web location of the 
#'   tracks' data files (`bigDataUrl` fields from the hub.txt file). For example, 
#'   if the hub should be hosted on GitHub, `bigDataUrl_base` should be something 
#'   like "https://github.com/morinlab/LLMPP/blob/main" (the default). See how the 
#'   track's paths from the bigDataUrl fields are set in the **Details** section.
#' @param bedToBigBed_path Path to your local `bedToBigBed` UCSC tool. If missing, 
#'   `GAMBLR.helpers::check_config_value` is called internally and the `bedToBigBed` 
#'   path is obtained from the `config.yml` file saved in the current working 
#'   directory.
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
#'   longLabel = "Public aSHM mutations separated by pathologies",
#'   contact_email = "rdmorin@sfu.ca"
#' )
#' }
#' 
build_browser_hub = function(maf_data,
                             regions_bed,
                             these_sample_ids = NULL,
                             these_samples_metadata = NULL,
                             these_seq_types = c("genome", "capture"),
                             projection = "grch37",
                             local_web_host_dir = NULL,
                             hub_dir = "my_hub",
                             splitColumnName = "pathology",
                             hub_name = basename(hub_dir),
                             shortLabel = basename(hub_dir),
                             longLabel = basename(hub_dir),
                             contact_email,
                             visibility = "squish",
                             bigDataUrl_base = "https://github.com/morinlab/LLMPP/blob/main",
                             bedToBigBed_path){
  
  # check some provided parameter 
  stopifnot("`these_seq_types` must be one or more of \"genome\" or \"capture\"." = 
              all(these_seq_types %in% c("genome", "capture")))
  stopifnot("`projection` must be one of \"grch37\" or \"hg38\"." = 
              projection %in% c("grch37", "hg38"))
  stopifnot("`visibility` must be one of \"pack\", \"dense\", \"full\", or \"squish\"." = 
              visibility %in% c("pack", "dense", "full", "squish"))
  stopifnot("`contact_email` must be provided. UCSC browser will require this field to upload the track hub." = 
              !missing(contact_email))
  
  if(missing(bedToBigBed_path)){
    bedToBigBed_path = tryCatch(
      GAMBLR.helpers::check_config_value(config::get("dependencies")$bedToBigBed),
      error=function(e){
        k = paste0("Since a `bedToBigBed_path` wasn't provided, `build_browser_hub` tries to use a config.yml file to get the bedToBigBed path. However...\n", e)
        stop(k, call. = FALSE)
      }
    )
  }else{
    stopifnot("`bedToBigBed_path` points to a non-existent file." = 
                file.exists(bedToBigBed_path))
  }
  
  # get metadata with the dedicated helper function (for each seq type)
  these_seq_types = setNames(these_seq_types, these_seq_types)
  these_samples_metadata = lapply(these_seq_types, function(these_seq_types_i){
    id_ease(these_samples_metadata = these_samples_metadata,
            these_sample_ids = these_sample_ids,
            verbose = FALSE,
            this_seq_type = these_seq_types_i)
  }) %>% 
    suppressMessages
  
  # check provided splitColumnName parameter 
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
  
  # define regions_bed
  if(!missing(regions_bed)){
    if(is.null(regions_bed)){
      warning("`regions_bed = NULL` is treated as it was not provided. The default aSHM regions of the provided `projection` are used.")
      rm(regions_bed)
    }
  }
  if(missing(regions_bed)){
    if(projection == "grch37"){
      regions_bed = GAMBLR.data::grch37_ashm_regions
    }else{
      regions_bed = GAMBLR.data::hg38_ashm_regions
    }
  }
  
  # save regions_bed to a bb file
  regions_bed = dplyr::select(regions_bed, 1,2,3) %>% 
    arrange( .[[1]], .[[2]] )
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
  regions_bb_file = file.path(track_dir, "regions.bb")
  bigbed_conversion = gettextf("%s %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, 
                               regions_bb_file)
  system(bigbed_conversion)
  unlink(temp_bed)
  unlink(temp_chr_sizes)
  
  # get maf data from the specified regions and metadata/samples (for each seq type)
  if(missing(maf_data)){
    maf_data = mapply(function(these_seq_types_i, these_samples_metadata_i){
      get_ssm_by_regions(regions_bed = regions_bed,
                         this_seq_type = these_seq_types_i,
                         these_samples_metadata = these_samples_metadata_i, 
                         projection = projection,
                         streamlined = FALSE, 
                         basic_columns = TRUE) %>% 
        suppressMessages
    }, these_seq_types, these_samples_metadata, SIMPLIFY = FALSE)
  }else{
    maf_data = mapply(function(these_seq_types_i, these_samples_metadata_i){
      get_ssm_by_regions(maf_data = maf_data,
                         regions_bed = regions_bed,
                         this_seq_type = these_seq_types_i,
                         these_samples_metadata = these_samples_metadata_i, 
                         projection = projection,
                         streamlined = FALSE, 
                         basic_columns = TRUE) %>% 
        suppressMessages
    }, these_seq_types, these_samples_metadata, SIMPLIFY = FALSE)
  }
  
  # split maf table according to splitColumnName
  if(!is.null(splitColumnName)){
    maf_data = mapply(function(maf_data_i, these_samples_metadata_i){
      maf_data_i = dplyr::select(these_samples_metadata_i, sample_id, all_of(splitColumnName)) %>% 
        distinct(sample_id, .keep_all = TRUE) %>% 
        left_join(maf_data_i, ., join_by(Tumor_Sample_Barcode == sample_id))
      dplyr::select(maf_data_i, -all_of(splitColumnName)) %>% 
        split(maf_data_i[[splitColumnName]])
    }, maf_data, these_samples_metadata, SIMPLIFY = FALSE)
    track_names = lapply(maf_data, names)
    these_samples_metadata = mapply(function(these_samples_metadata_i, track_names_i){
      split(these_samples_metadata_i, these_samples_metadata_i[[splitColumnName]]) %>% 
        "["(track_names_i)
    }, these_samples_metadata, track_names, SIMPLIFY = FALSE)
  }else{
    maf_data = lapply(maf_data, list)
    track_names = lapply(these_seq_types, function(x) list("all"))
    these_samples_metadata = lapply(these_samples_metadata, list)
  }
  
  # check seq types that ssms could be found
  is_there_muts = lengths(track_names) > 0
  if( !any(is_there_muts) ){
    stop("No SSMs were found. Reset your parameters.")
  }
  maf_data = maf_data[is_there_muts]
  these_samples_metadata = these_samples_metadata[is_there_muts]
  these_seq_types = these_seq_types[is_there_muts]
  track_names = track_names[is_there_muts] %>% 
    mapply(paste, these_seq_types, ., sep = "_", SIMPLIFY = FALSE)
  
  # convert and save track files
  track_file_names = lapply(track_names, paste0, ".bb")
  mapply(
    function(maf_i, meta_i, these_seq_types_i, track_names_file_i){
      mapply(maf_to_custom_track, 
             maf_data = maf_i,
             these_samples_metadata = meta_i,
             this_seq_type = these_seq_types_i,
             output_file = file.path(track_dir, track_names_file_i),
             as_bigbed = TRUE,
             projection = projection,
             bedToBigBed_path = bedToBigBed_path)
    },
    maf_data, these_samples_metadata, these_seq_types, track_file_names
  ) %>% 
    invisible
  
  
  ### create hub.txt file
  
  # open file
  hub_file = paste0(projection, "_hub.txt") %>% 
    file.path(hub_dir_full_path, .)
  sink(hub_file)
  
  # write header
  cat( paste0("hub ", hub_name, "\n") )
  cat( paste0("shortLabel ", shortLabel, "\n") )
  cat( paste0("longLabel ", longLabel, "\n") )
  cat( "useOneFile on\n" )
  cat( paste0("email ", contact_email, "\n") )
  cat( "\n" )
  cat( paste0("genome ", replace(projection, projection == "grch37", "hg19"), "\n") )
  
  # write info of the track which stores regions from where ssms are retrieved
  cat( "\n" )
  cat( "track regions\n" )
  cat( "shortLabel regions\n" )
  cat( "longLabel Regions of interest from where mutations were retrieved\n" )
  cat( paste0("visibility ", visibility, "\n") )
  cat( "priority 1\n" )
  cat( paste0("type bigBed\n") )
  file.path(bigDataUrl_base, hub_dir, projection, basename(regions_bb_file)) %>% 
    { cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
  
  # write info of the tracks that store ssms split by splitColumnName
  track_names = unlist(track_names)
  track_file_names = unlist(track_file_names)
  for(i in seq_along(track_names)){
    cat( "\n" )
    cat( paste0("track ", track_names[i], "\n") )
    cat( paste0("shortLabel ", track_names[i], "\n") )
    cat( paste0("longLabel ", track_names[i], "\n") )
    cat( paste0("visibility ", visibility, "\n") )
    cat( paste0("priority ", i+1, "\n") )
    cat( paste0("type bigBed 9\n") )
    cat( "itemRgb on\n" )
    file.path(bigDataUrl_base, hub_dir, projection, track_file_names[i]) %>% 
      { cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
  }
  
  # close file
  sink()
}
