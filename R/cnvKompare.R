#' @title Compare segmented data for multiple samples.
#'
#' @description `cnvKompare` returns a list in variable data formats allowing to evaluate concordance of CNV data between multiple samples.
#'
#' @details This function will compare CNV data between samples with multiple time points. It can also handle same-sample comparison
#' between different CNV callers if the sample ID is specified in unique fashion. For groups with more than 2 samples,
#' optionally pairwise comparisons can be performed. The comparison is made based on the internally calculated score,
#' which reflects the percentage of each cytoband covered by CNV (rounded to the nearest 5\%) and its absolute CN. Optionally,
#' the heatmap of cnvKompare scores can be returned. In addition, the function will return all concordant and discordant cytobands.
#' Finally, the time series plot of CNV log ratios will be returned for all lymphoma genes, with further functionality to subset
#' it to a panel of genes of interest.
#'
#' @param patient_id Specify patient_id to retrieve sample ids from GAMBL metadata.
#' @param these_sample_ids Optionally, specify sample ids for comparison.
#' @param this_seg Optional input data frame of seg file. Must adhere to seg format.
#' @param seg_path Optionally, specify the path to a local seg file. Must adhere to seg format.
#' @param genes_of_interest Provide specific genes to be displayed on the time-series plot.
#' @param projection Argument specifying the projection of seg file, which will determine coordinates of the cytobands. Default is grch37, but hg38 is also accepted.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param ignore_cytoband_labels Cytobands to be ignored. By default, "acen", "gvar", "stalk" are excluded.
#' @param max_overlap For a time-series plot, how many maximum overlapping points are allowed?
#' @param min_concordance Integer value from 0 to 100 to indicate the minimum required similarity between cytobands to be considered concordant. The default is 90 (90 %).
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is FALSE.
#' @param return_heatmap Boolean argument specifying whether to return a heatmap of cnvKompare scores. Default is TRUE.
#' @param compare_pairwise Boolean argument specifying whether to perform pairwise comparisons if there are more than 2 time points in the group. Default is TRUE.
#' @param show_x_labels Optional boolean parameter for hiding/showing x axis labels, default is TRUE.
#'
#' @return A list of overall and pairwise percent concordance, concordant and discordant cytobands, comparison heatmap of cnvKompare scores, and time series ggplot object.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr circlize ggplot2 ggrepel readr tibble GAMBLR.helpers
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#'
#' #get segs for comparison
#' my_segs = get_sample_cn_segments(these_sample_ids = c("02-13135T", "04-28140T"))
#'
#' my_list = cnvKompare(this_seg = my_segs,
#'                      these_sample_ids = c("02-13135T", "04-28140T"),
#'                      genes_of_interest = c("EZH2", "TP53", "MYC", "CREBBP","GNA13"),
#'                      show_x_labels = FALSE)
#'
cnvKompare = function(patient_id,
                      these_sample_ids,
                      this_seg,
                      seg_path,
                      genes_of_interest,
                      projection = "grch37",
                      this_seq_type = "genome",
                      ignore_cytoband_labels = c("acen", "gvar", "stalk"),
                      max_overlap = 20,
                      min_concordance = 90,
                      exclude_sex = FALSE,
                      return_heatmap = TRUE,
                      compare_pairwise = TRUE,
                      show_x_labels = TRUE){

  # initialize output list
  output = list()

  # convert min concordance pct to fraction
  min_concordance = min_concordance/100

  # check that sample identifiers are provided
  if (missing(patient_id) & missing(these_sample_ids)) {
    stop("Please provide patient id or sample ids for comparison.")
  }

  # retrieve sample ids if only patient id is specified
  if (missing(these_sample_ids)) {
    these_sample_ids = GAMBLR.helpers::handle_metadata(this_seq_type = this_seq_type)
    these_sample_ids = dplyr::filter(these_sample_ids, patient_id == {{ patient_id }})
    these_sample_ids = pull(these_sample_ids, sample_id)
    message(paste0(
      "Found ",
      length(these_sample_ids),
      " samples for patient ",
      patient_id,
      " ..."
    ))
  }

  # get cytobands
  if (projection %in% c("hg19", "grch37")) {
    cytobands = circlize::read.cytoband(species = "hg19")$df %>%
      mutate(V1 = gsub("chr", "", V1))
  } else if (projection %in% c("hg38", "grch38")) {
    cytobands = circlize::read.cytoband(species = "hg38")$df
  } else {
    stop("Please specify one of hg19, grch37, hg38, or grch38 projections.")
  }
  cytobands = cytobands %>%
    `names<-`(c("cb.chromosome", "cb.start", "cb.end", "cb.name", "label")) %>%
    dplyr::filter(!label %in% ignore_cytoband_labels)
  if (exclude_sex) {
    cytobands = dplyr::filter(cytobands,!grepl("X|Y", cb.chromosome))
  }
  cytobands = as.data.table(cytobands)
  setkey(cytobands, cb.chromosome, cb.start, cb.end)

  # get the multi-sample seg file
  if (!missing(seg_path)) {
    message(paste0("Reading seg file from: ", seg_path))
    these_samples_seg = suppressMessages(read_tsv(seg_path))
  } else if (!missing(this_seg)) {
    message("Using supplied seg file")
    these_samples_seg = this_seg
  } else {
    message("Retreiving the CNV data using GAMBLR ...")
    these_samples_seg = get_sample_cn_segments(these_sample_ids = these_sample_ids,
                                               projection = projection,
                                               this_seq_type = this_seq_type)
  }

  #add the CN column, if not already there
  if(!("CN" %in% colnames(these_samples_seg))){
    these_samples_seg = dplyr::mutate(these_samples_seg, CN = (2 * 2 ^ log.ratio))
  }

  these_samples_seg = these_samples_seg  %>%
    dplyr::filter(ID %in% these_sample_ids) %>% # if user-provided seg, ensure only samples of comparison are present
    relocate(ID, .after = last_col())
  if (exclude_sex) {
    these_samples_seg = dplyr::filter(these_samples_seg,!grepl("X|Y", chrom))
  }

  these_samples_seg = as.data.table(these_samples_seg)
  setkey(these_samples_seg, chrom, start, end)

  # overlap seg with cytoband regions
  # if segment extends beyond cytoband, cut it at the cytoband coordinates
  cytoband_overlap =
    foverlaps(cytobands,
              these_samples_seg,
              nomatch = 0) %>%
    as.data.frame %>%
    group_by(cb.name) %>%
    mutate(
      start = ifelse(start > cb.start, start, cb.start),
      end = ifelse(end < cb.end, end, cb.end)
    ) %>%
    ungroup

  # calculate % of cytoband covered by CNV and concordance score
  message("Calculating CNV concordance ...")
  for_output =
    cytoband_overlap %>%
    mutate(
      band_length = cb.end - cb.start,
      cnv_length = end - start,
      # round the % of cytoband covered to the nearest 5%
      pct_covered = (cnv_length / band_length * 100),
      pct_covered = 5*round(pct_covered/5)
    ) %>%
    # name each cytoband as chr_cytoband
    mutate(name = paste0(cb.chromosome, "_", cb.name)) %>%
    # calculate % covered by cytoband and it's CN
    group_by(ID, name) %>%
    mutate(pct_covered = sum(pct_covered),
           CN = mean(CN)) %>%
    ungroup() %>%
    arrange(name, ID) %>%
    distinct(ID, CN, pct_covered, name, .keep_all = TRUE) %>%
    # designate score of %covered*CN to infer intra-sample concordance
    mutate(score = CN * pct_covered)

  # overall concordance
  concordance =
    for_output %>%
    select(ID, name, score) %>%
    spread(., ID, score) %>%
    column_to_rownames("name")

  overall_concordance = ifelse((rowSums(concordance) >= ((concordance[, 1] * ncol(concordance))*min_concordance) &
                                  rowSums(concordance) <= ((concordance[, 1] * ncol(concordance))*(1+(1-min_concordance)))),
                               "YES",
                               "NO")
  overall_concordance = overall_concordance[!is.na(overall_concordance)]
  overall_concordance_pct = round(((
    sum(overall_concordance == "YES") / length(overall_concordance)
  ) * 100), 2)
  output$overall_concordance_pct = overall_concordance_pct

  # return cytobands consistent across samples
  concordant_cytobands =
    for_output %>%
    # output-specific
    select(ID, cb.chromosome, cb.start, cb.end, name, score, log.ratio) %>%
    dplyr::filter(name %in% names(overall_concordance[overall_concordance == "YES"]))

  output$concordant_cytobands = concordant_cytobands

  # return cytobands discordant across samples
  discordant_cytobands =
    for_output %>%
    # output-specific
    select(ID,
           cb.chromosome,
           cb.start,
           cb.end,
           name,
           pct_covered,
           log.ratio,
           score) %>%
    dplyr::filter(name %in% names(overall_concordance[overall_concordance == "NO"]))

  output$discordant_cytobands = discordant_cytobands

  # heatmap of cnvKompare scores
  if (return_heatmap) {
    message("Building heatmap ...")
    hMap <- concordance %>%
        as.data.frame %>%
        rownames_to_column("x") %>%
        pivot_longer(!x, names_to = "y", values_to = "score") %>%
        ggplot(aes(x = x, y = y, fill = score)) +
        geom_tile() +
            scale_fill_gradientn(
                colors = hcl.colors(10, "Blue-Yellow")
            ) +
        guides(
            fill = guide_colourbar(
                title = "cnvKompare score"
            )
        )  +
        labs(
            y = "",
            x = "Cytoband"
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.text.y = element_text(
                face = "bold",
                color = "black"
            ),
            axis.text.x = element_blank(),
            axis.title = element_text(
                face = "bold"
            ),
            panel.background = element_blank()
        )

    output$Heatmap = hMap
  }

  # VAF-like plot
  # genes
  if (projection %in% c("hg19", "grch37")) {
    for_plot_lg = GAMBLR.data::grch37_lymphoma_genes_bed %>%
      as.data.table()
  } else if (projection %in% c("hg38", "grch38")) {
    for_plot_lg = GAMBLR.data::hg38_lymphoma_genes_bed %>%
      as.data.table()
  }
  # did user specify particular genes of interest to display on the plot?
  if (!missing(genes_of_interest)) {
    message("Subsetting lymphoma genes to specified genes of interest ...")
    for_plot_lg = dplyr::filter(for_plot_lg, hgnc_symbol %in% genes_of_interest)
  }
  setkey(for_plot_lg, chromosome_name, start_position, end_position)

  # CNV data
  for_plot = as.data.table(for_output)
  setkey(for_plot, cb.chromosome, start, end)

  # generate plot
  time_plot =
    foverlaps(for_plot_lg,
              for_plot,
              nomatch = NULL) %>%
    as.data.frame %>%
    select(ID, hgnc_symbol, log.ratio) %>%
    ggplot(., aes(x = ID, y = log.ratio, group = hgnc_symbol)) +
    geom_line(aes(color = hgnc_symbol)) +
    ggrepel::geom_label_repel(
      aes(label = hgnc_symbol, color = hgnc_symbol),
      nudge_x = 0.1,
      na.rm = TRUE,
      label.size = NA,
      fill = NA,
      segment.color = "transparent",
      max.overlaps = max_overlap
    ) +
    GAMBLR.helpers::theme_Morons() +
    theme(legend.position = "none",
          axis.title.x = element_blank())

  output$time_plot = time_plot

  # for groups with >2 samples, make pairwise comparisons
  if (compare_pairwise & length(these_sample_ids) > 2) {
    message("Performing pairwise comparisons ...")

    # generate all possible combinations
    possible_combinations = apply(combn(these_sample_ids, 2), 2, paste, collapse =
                                    '--')

    for (combination in possible_combinations) {
      # samples in this pairwise comparison
      these_samples = unlist(strsplit(combination, split = "--"))

      # pct concordance in this pair
      this_concordance = concordance[, these_samples]
      this_concordance = ifelse((this_concordance[, 1] >= (this_concordance[, 2])*min_concordance) &
                                  this_concordance[, 1] <= (this_concordance[, 2])*(1+(1-min_concordance)),
                                "YES",
                                "NO")
      names(this_concordance) = rownames(concordance)
      this_concordance = this_concordance[!is.na(this_concordance)]
      this_concordance_pct = round(((
        sum(this_concordance == "YES") / length(this_concordance)
      ) * 100), 2)
      output$pairwise_comparisons[[combination]]$pairwise_concordance_pct = this_concordance_pct

      # return cytobands consistent in this pair
      concordant_cytobands =
        for_output %>%
        dplyr::filter(ID %in% these_samples) %>%
        # output-specific
        select(ID, cb.chromosome, cb.start, cb.end, name, score) %>%
        dplyr::filter(name %in% names(this_concordance[this_concordance == "YES"]))

      output$pairwise_comparisons[[combination]]$concordant_cytobands = concordant_cytobands

      # return cytobands discordant in this pair
      discordant_cytobands =
        for_output %>%
        dplyr::filter(ID %in% these_samples) %>%
        # output-specific
        select(ID,
               cb.chromosome,
               cb.start,
               cb.end,
               name,
               pct_covered,
               log.ratio,
               score) %>%
        dplyr::filter(name %in% names(this_concordance[this_concordance == "NO"]))

      output$pairwise_comparisons[[combination]]$discordant_cytobands = discordant_cytobands

    }
  }
  message("DONE!")
  return(output)
}
