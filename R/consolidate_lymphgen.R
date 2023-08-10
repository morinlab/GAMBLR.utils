#' @title Consolidate Lymphgen.
#'
#' @description Replace the lymphgen column in the incoming metadata with classification for additional samples.
#'
#' @details Supplement the "lymphgen" column of the metadata with classification for additional samples.
#' Expects at least to have columns "patient_id" to bind on, and "lymphgen" to supplement the data on.
#'
#' @param sample_table Input data frame with metadata.
#' @param derived_data_path Optional argument specifying the path to a folder with files following the pattern *lymphgen.txt.
#' @param verbose Default is TRUE.
#'
#' @return A data frame with a supplemented lymphGen column.
#'
#' @import dplyr purrr readr
#' @export
#'
#' @examples
#' metadata = get_gambl_metadata()
#' consolidate_lymphgen(sample_table = metadata)
#'
consolidate_lymphgen = function(sample_table,
                                derived_data_path = "",
                                verbose = TRUE){

  if (derived_data_path == "") {
    path_to_files = check_config_value(config::get("derived_and_curated"))
    project_base = check_config_value(config::get("project_base"))
    derived_data_path = paste0(project_base, path_to_files)
    if (verbose) {
      message(
        paste0(
          "No external data path was provided, using default path ",
          derived_data_path
        )
      )
    }
  } else{
    if (verbose) {
      message(paste0("Using the specified path ", derived_data_path))
    }
  }

  lymphgen_files = dir(derived_data_path, pattern = "*lymphgen.txt")
  if (length(lymphgen_files) > 0) {
    if (verbose) {
      message(paste0(
        "Found these file(s) with lymphgen information: ",
        lymphgen_files
      ))
    }
  } else{
    if (verbose) {
      message(
        paste0(
          "No file(s) with lymphgen information were found at the path",
          lymphgen_files
        )
      )
      message(
        "If you expected the data to be present at the specified location, please ensure they follow naming convention *.lymphgen.txt"
      )
    }
  }

  for (f in lymphgen_files) {
    full = paste0(project_base, path_to_files, f)
    this_data = suppressMessages(read_tsv(full, comment = "#"))
    sample_table =
      sample_table %>% left_join(this_data, by = 'patient_id', suffix = c(".X", ".Y")) %>%
      split.default(gsub('.[XY]', '', names(.))) %>%
      map_dfc(~ if (ncol(.x) == 1)
        .x
        else
          mutate(.x, !!sym(gsub(
            '.X', '', names(.x)[1]
          )) := coalesce(!!!syms(names(
            .x
          ))))) %>%
      select(colnames(sample_table))
    sample_table = tidy_lymphgen(sample_table,
                                 lymphgen_column_in = "lymphgen",
                                 lymphgen_column_out = "lymphgen",
                                 relevel=TRUE)
  }

  return(sample_table)
}
