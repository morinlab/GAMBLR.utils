#' @title Generate md5 Hash For Samples
#'
#' @description Generate an md5 hash for a set of samples to help ensure reproducibility
#'
#' @details This function can accept a wide range of formatted sample IDs to create an md5 hash.
#' For example, if the user is working with an already subset metadata table (with sample IDs of interest),
#' The user can give this table to the function with `these_sampels_metadata`.
#' As an alternative, sample IDs can also be provided as a vector of characters with `these_samples` parameter.
#' Another option is to use defined sample sets (GAMBL) with `sample_set_name`.
#' As a final option, the user can also provide a data frame with samples IDs instead of loading them from the GAMBL repo,
#' This is achieved with calling the `sample_sets_df` parameter.
#'
#'
#' @param these_samples_metadata Optionally provide a metadata table or any data frame with a column named sample_id that has been subset to the samples you're working with.
#' @param these_samples Optionally provide a vector of sample_id you are working with.
#' @param sample_set_name Optionally provide the name of a sample set in GAMBL and the function will load the samples from that set and provide the hash.
#' @param sample_sets_df Optionally provide a data frame of the sample sets instead of relying on/loading the local file from the GAMBL repo.
#'
#' @return The md5 hash of the ordered set of sample_id.
#'
#' @import digest dplyr readr
#' @export
#'
get_samples_md5_hash = function(these_samples_metadata,
                                these_samples,
                                sample_set_name,
                                sample_sets_df){

  if(!missing(these_samples_metadata)){
    collapsed = dplyr::select(these_samples_metadata,sample_id) %>%
      arrange() %>%
      pull() %>% paste(.,collapse=",")

      digested = digest::digest(collapsed,serialize = FALSE)
  }else if(!missing(these_samples)){
    digested = digest::digest(paste(these_samples[order(these_samples)],collapse=","),serialize=FALSE)
  }else if(!missing(sample_set_name)){
    #load the sample set table and pull the samples based on its contents and the name provided
    sample_sets_file = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("sample_sets")$default))
    if(missing(sample_sets_df)){
      sample_sets = suppressMessages(read_tsv(sample_sets_file))
    }else{
      sample_sets = sample_sets_df
    }
    setname = as.symbol(sample_set_name)
    collapsed = dplyr::filter(sample_sets,!!setname==1) %>%
      dplyr::select(sample_id) %>%
      arrange() %>% pull() %>% paste(.,collapse=",")
    digested = digest::digest(collapsed,serialize=FALSE)
  }
  return(digested)
}
