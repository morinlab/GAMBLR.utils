#' @title Calculate Mutation Frequency By Sliding Window.
#'
#' @description Count the number of mutations in a sliding window across a region for all samples.
#'
#' @details This function is called to return the mutation frequency for a given region, for all GAMBL samples. Regions are specified with the `this_region`parameter.
#' Alternatively, the region of interest can also be specified by calling the function with `chromosome`, `start_pos`, and `end_pos` parameters.
#' It is also possible to return a plot of the created bins. This is done with setting `plot_type = TRUE`.
#' There are a collection of parameters available for further customizing the return, for more information, refer to the parameter descriptions and examples.
#' This function is unlikely to be used directly in most cases. See [GAMBLR::get_mutation_frequency_bin_matrix] instead.
#'
#' @param this_region Genomic region in bed format.
#' @param chromosome Chromosome name in region.
#' @param start_pos Start coordinate of region.
#' @param end_pos End coordinate of region.
#' @param metadata Data frame containing sample ids and column with annotated data for the 2 groups of interest. All other columns are ignored. Currently, function exits if asked to compare more than 2 groups.
#' @param seq_type The seq_type you want back, default is genome.
#' @param slide_by Slide size for sliding window, default is 100.
#' @param window_size Size of sliding window, default is 1000.
#' @param plot_type Set to TRUE for a plot of your bins. By default no plots are made.
#' @param sortByColumns Which of the metadata to sort on for the heatmap
#' @param return_format Return format of mutations. Accepted inputs are "long" and "long-simple". Default is "long-simple".
#' @param min_count_per_bin Minimum counts per bin, default is 3.
#' @param return_count Boolean statement to return count. Default is FALSE.
#' @param drop_unmutated This may not currently work properly. Default is FALSE.
#' @param classification_column Only used for plotting, default is "lymphgen"
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flat-files (only works for streamlined data, not full MAF details). Default is FALSE.
#' @param mode Only works with indexed flat-files. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#' @return Count matrix.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr cowplot ggplot2
#' @export
#'
#' @examples
#' chr11_mut_freq = calc_mutation_frequency_sliding_windows(this_region = "chr11:69455000-69459900",
#'                                                          slide_by = 10,
#'                                                          window_size = 10000)
#'
calc_mutation_frequency_sliding_windows = function(this_region,
                                                   chromosome,
                                                   start_pos,
                                                   end_pos,
                                                   metadata,
                                                   seq_type,
                                                   slide_by = 100,
                                                   window_size = 1000,
                                                   plot_type = "none",
                                                   sortByColumns = "pathology",
                                                   return_format = "long-simple",
                                                   min_count_per_bin = 3,
                                                   return_count = FALSE,
                                                   drop_unmutated = FALSE,
                                                   classification_column = "lymphgen",
                                                   from_indexed_flatfile = FALSE,
                                                   mode = "slms-3"){

  max_region = 1000000
  if(missing(metadata)){
    metadata = get_gambl_metadata()
  }
  if(missing(this_region)){
    this_region = paste0(chromosome, ":", start_pos, "-", end_pos)
  }else{
    chunks = region_to_chunks(this_region)
    chromosome = chunks$chromosome
    start_pos = as.numeric(chunks$start)
    end_pos = as.numeric(chunks$end)
  }
  region_size = end_pos - start_pos
  if(region_size < max_region){
    message(paste("processing bins of size", window_size, "across", region_size, "bp region"))
  }else{
    message(paste("CAUTION!\n", region_size, "exceeds maximum size recommended by this function."))
  }
  windows = data.frame(start = seq(start_pos, end_pos, by = slide_by)) %>%
    mutate(end = start + window_size - 1)
  #use foverlaps to assign mutations to bins
  windows.dt = as.data.table(windows)
  region_ssm = GAMBLR::get_ssm_by_region(region = this_region, streamlined = FALSE, seq_type=seq_type, from_indexed_flatfile = from_indexed_flatfile, mode = mode) %>%
    dplyr::rename(c("start" = "Start_Position", "sample_id" = "Tumor_Sample_Barcode")) %>%
    mutate(mutated = 1)

  region.dt = region_ssm %>%
    dplyr::mutate(start = as.numeric(as.character(start)), end = start + 1, end = as.numeric(as.character(end))) %>%
    dplyr::relocate(start, .before=end) %>%
    as.data.table()

  setkey(windows.dt, start, end)
  setkey(region.dt, start, end)
  windows_overlap = foverlaps(windows.dt, region.dt) %>%
    dplyr::filter(!is.na(start)) %>%
    dplyr::rename(c("window_start" = "i.start", "mutation_position" = "start")) %>%
    dplyr::select(-i.end, -end, -mutation_position) %>%
    as.data.frame()

  windows_tallied_full = windows_overlap %>%
    group_by(sample_id, window_start) %>%
    tally() %>%
    dplyr::filter(n >= min_count_per_bin) %>%
    arrange(sample_id) %>%
    as.data.frame()
  windows_tallied = windows_tallied_full

  all_samples = pull(metadata, sample_id) %>%
    unique()

  num_samples = length(all_samples)
  lg_cols = get_gambl_colours("lymphgen")
  path_cols = get_gambl_colours("pathology")
  annos = data.frame(window_start = rep(start_pos,num_samples), sample_id = factor(all_samples))
  annos = left_join(annos, metadata, by = "sample_id")
  windows_tallied = left_join(metadata, windows_tallied, by = "sample_id")
  windows_tallied = arrange(windows_tallied, lymphgen)
  windows_tallied$classification = factor(windows_tallied[,classification_column], levels = unique(windows_tallied[,classification_column]))
  if(drop_unmutated){
    windows_tallied = windows_tallied %>%
      dplyr::filter(!is.na(n))
  }
  if(classification_column == "lymphgen"){
    windows_tallied = arrange(windows_tallied, pathology, lymphgen)
    annos = arrange(annos, pathology, lymphgen)
  }else{
    windows_tallied = arrange(windows_tallied, classification)
    annos = arrange(annos, classification)
  }
  annos$sample_id = factor(annos$sample_id, levels = unique(annos$sample_id))
  windows_tallied$sample_id = factor(windows_tallied$sample_id, levels = unique(windows_tallied$sample_id))
  if(plot_type %in% c("points", "point")){
    #add a bin at position 1 for pathology
    windows_tallied = dplyr::filter(windows_tallied, !is.na(window_start))
    p = ggplot2::ggplot() +
                 geom_point(data = annos, aes(x = window_start, y = sample_id, colour = pathology)) +
                 geom_point(data = windows_tallied, aes(x = window_start, y = sample_id, colour = classification)) +
                 theme(axis.text = element_text(size = 4)) +
                 theme(axis.text.y = element_blank()) +
                 scale_colour_manual(values = c(lg_cols, path_cols))
  }else if(plot_type %in% c("tile", "tiled")){
    print(annos)
    p = ggplot() +
        geom_point(data = annos, aes(x = window_start, y = sample_id, colour = lymphgen)) +
        geom_tile(data = windows_tallied, aes(x = window_start, y = sample_id, fill = n)) +
        scale_fill_gradient(low = "orange", high = "red", na.value = NA) +
        theme_cowplot() +
        scale_colour_manual(values = c(lg_cols, path_cols)) +
        theme(axis.text.y = element_blank())
  }
  if(plot_type != "none"){
    print(p)
  }
  if(return_count){
    windows_tallied = mutate(windows_tallied, bin = paste0(window_start, "_", chromosome)) %>%
      mutate(mutated = n)
  }else{
    #binary mutated or not
    windows_tallied = mutate(windows_tallied, bin = paste0(window_start, "_", chromosome)) %>%
    mutate(mutated = 1)
  }
  a = dplyr::select(windows_tallied, sample_id, bin, mutated)
  completed = complete(a, sample_id, bin, fill = list(mutated = 0))
  widened = pivot_wider(completed, names_from = sample_id, values_from = mutated)
  if(return_format == "long"){
    return(windows_tallied)
  }else if(return_format == "long-simple"){
    #just return the columns needed for completing and making a wide matrix for many regions
    windows_simple = dplyr::select(windows_tallied, sample_id, bin, mutated)
    return(windows_simple)
  }else{
    return(widened)
  }
}
