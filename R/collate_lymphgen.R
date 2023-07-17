#' @title Collate Lymphgen.
#'
#' @description Expand a sample_table (metadata) horizontally with different flavours of lymphgen data.
#'
#' @details This function takes a sample table (metadata) and adds different flavours of lymphgen data.
#' It is possible to call this function with an already subset metadata table (with sample IDs of interest) with `these_samples_metadata`.
#' If this is done, the function will join the lymphgen data with this table. Currently, the only supported `lymphgen_version` is "default".
#' For more information refer to the function examples.
#'
#' @param these_samples_metadata Optional parameter with metadata filtered for sample_ids of interest. If provided, this function will join lymphgen with this metadata, regardless of tidy TRUE/FALSE.
#' @param lymphgen_version Version of selected lymphgen, default is "default".
#' @param tidy Boolean parameter, set to TRUE for tidy format (i.e long format with no columns dropped). Default is FALSE, which returns the data in a wide format, keeping both the original Subtype. Prediction and tidied LymphGen values and puts the values from each "flavour" in its own column.
#'
#' @return A df with lymphgen information.
#'
#' @import dplyr tidyr readr stringr
#' @export
#'
#' @examples
#' this_meta = get_gambl_metadata()
#' dlbcl_meta = dplyr::filter(this_meta, pathology == "DLBCL")
#'
#' wide_lymphgen = collate_lymphgen(these_samples_metadata = dlbcl_meta,
#'                                  lymphgen_version = "default",
#'                                  tidy = FALSE)
#'
collate_lymphgen = function(these_samples_metadata,
                            lymphgen_version = "default",
                            tidy = FALSE){

  #TODO Update the key in the config to match the version once updated, as discussed on PR.
  if(lymphgen_version == "default"){
    lymphgen_template = check_config_value(config::get("results_versioned")$lymphgen_template$default)
  }else{
    stop("Currently, only lymphgen_version = default is accepted")
  }

  #repo base
  repo_base = check_config_value(config::get("repo_base"))
  flavours = check_config_value(config::get("results_merged_wildcards")$lymphgen_template)
  flavour = str_split(flavours, pattern = ",")
  flavour = unlist(flavour)
  lymphgen_path = paste0(repo_base, lymphgen_template)

  load_lymphgen = function(flavour, lymphgen_path){
    lg_path = glue::glue(lymphgen_path)
    if(!file.exists(lg_path)){ #ignore missing flavours.
      return()
    }
    lg_df = suppressMessages(read_tsv(lg_path)) %>%
      mutate(flavour = flavour) #append the flavour in its own column called "flavour".
    return(lg_df)
  }

  lymphgen_results = lapply(flavour, load_lymphgen, lymphgen_path = lymphgen_path)
  lymphgen_results = bind_rows(lymphgen_results) #get lymphgen results tables stacked on top of each other, with the results from each flavour identified by the `flavour` column.
  lymphgen_results = tidy_lymphgen(lymphgen_results, lymphgen_column_in = "Subtype.Prediction", lymphgen_column_out = "LymphGen")
  colnames(lymphgen_results)[1] = "sample_id"

  if(!tidy){
    lymphgen_untidy = lymphgen_results %>%
      select(sample_id, Subtype.Prediction, LymphGen, flavour) %>%
      pivot_wider(names_from = flavour,
                  values_from = c(Subtype.Prediction, LymphGen),
                  names_glue = "{.value}_{flavour}")

      if(!missing(these_samples_metadata)){
        meta_data = these_samples_metadata
        lymphgen_untidy = left_join(meta_data, lymphgen_untidy)
      }

      return(lymphgen_untidy)

    }else{
    if(!missing(these_samples_metadata)){
      lymphgen_results = left_join(these_samples_metadata, lymphgen_results)
    }

    return(lymphgen_results)
  }
}
