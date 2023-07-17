#' @title Annotate Hotspots.
#'
#' @description Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations.
#'
#' @details This function takes an already loaded MAF data frame with the `mutation_maf` parameter.
#' The user can then control the minimum number of recurrences for mutations to be included with `recurrance_min`,
#' The default is 5. `analysis_base` controls the base name go hotspot output directory.
#' Lastly, `p_thresh` sets the p value threshold, default is 0.05.
#'
#' @param mutation_maf A data frame in MAF format.
#' @param recurrence_min minimum number of recurrences for mutation to be included, default is 5.
#' @param analysis_base Base name for hot spot output directory.
#' @param p_thresh P value threshold, default is 0.05.
#'
#' @return The same data frame with one additional column "hot_spot".
#'
#' @import dplyr tidyr readr
#' @export
#'
#' @examples
#' my_metadata = get_gambl_metadata()
#' all_coding_ssm = get_coding_ssm(these_samples_metadata = my_metadata,
#'                                 projection = "grch37",
#'                                 seq_type = "genome")
#'
#' hot_ssms = annotate_hotspots(all_coding_ssm)
#'
annotate_hotspots = function(mutation_maf,
                             recurrence_min = 5,
                             analysis_base = c("FL--DLBCL", "BL--DLBCL"),
                             p_thresh = 0.05){

  hotspot_info = list()
  for(abase in analysis_base){
    base_path = check_config_value(config::get("repo_base"))

    clust_full_path = paste0(base_path, check_config_value(config::get("results_versioned")$oncodriveclustl$clusters))
    clust_full_path = glue::glue(clust_full_path)
    all_full_path = paste0(base_path, check_config_value(config::get("results_versioned")$oncodriveclustl$elements))
    all_full_path = glue::glue(all_full_path)
    clust_hotspot = suppressMessages(readr::read_tsv(clust_full_path))
    all_hotspot = suppressMessages(readr::read_tsv(all_full_path))

  clustered_hotspots = clust_hotspot %>%
    dplyr::select(-RANK) %>%
    dplyr::filter(N_SAMPLES > recurrence_min & P < p_thresh)

    arranged = clustered_hotspots %>%
      separate_rows(COORDINATES, convert = TRUE) %>%
      group_by(SYMBOL, MAX_COORD) %>%
      arrange(COORDINATES)

    mins = arranged %>%
      slice_head() %>%
      dplyr::rename("START" = "COORDINATES")

    maxs = arranged %>%
      slice_tail() %>%
      dplyr::rename("END" = "COORDINATES")

    hotspot_ranges = left_join(mins, dplyr::select(maxs, c(MAX_COORD, END)), by = c("SYMBOL", "MAX_COORD"))
    hotspot_info[[abase]] = hotspot_ranges
  }
  merged_hotspot = do.call("rbind", hotspot_info) %>%
    ungroup()

  long_hotspot = merged_hotspot %>%
    dplyr::select(MAX_COORD, CHROMOSOME, START, END) %>%
    pivot_longer(c(START, END), names_to = "which", values_to = "COORDINATE") %>%
      dplyr::select(-which)

  #again take highest and lowest value for each MAX_COORD
  starts = long_hotspot %>%
    group_by(MAX_COORD) %>%
    arrange(COORDINATE) %>%
    slice_head()

  ends = long_hotspot %>%
    group_by(MAX_COORD) %>%
    arrange(COORDINATE) %>%
    slice_tail()

  long_hotspot = bind_rows(starts, ends)
  filled_coords = long_hotspot %>%
    group_by(MAX_COORD) %>%
    arrange(MAX_COORD, COORDINATE) %>%
    complete(COORDINATE = seq(COORDINATE[1], COORDINATE[2])) %>%
    tidyr::fill(CHROMOSOME, .direction = "up") %>%
    dplyr::rename("Start_Position" = "COORDINATE") %>%
    dplyr::rename("Chromosome" = "CHROMOSOME") %>%
    ungroup()

  filled_coords = mutate(filled_coords, hot_spot = TRUE)
  #just the ssms that match these coordinates!
  hot_ssms = left_join(mutation_maf, filled_coords, by = c("Chromosome", "Start_Position"))
  return(hot_ssms)
}
