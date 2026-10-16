#' @title Calculate proportion of genome altered by CNV.
#'
#' @description `calculate_pga` returns a data.frame with estimated proportion
#'      of genome altered for each sample.
#'
#' @details This function calculates the percent of genome altered (PGA) by CNV.
#' It takes into account the total length of sample's CNV and relates it to the
#' total genome length to return the proportion affected by CNV. In addition,
#' there is an option to calculate PGA per each chromosome individually by
#' setting the argument "per_chromosome" to TRUE.
#' The input is expected to be a seg file. The path to a local SEG file can also
#' be provided. If The custom seg file is provided, the minimum required columns
#' are sample (or ID), chrom, start, end, and log.ratio (or CN with
#' `value_scale = "CN"`).
#' The function can work with either individual or multi-sample seg files.
#' The telomeres are always excluded from calculation, and centromeres/sex
#' chromosomes can be optionally included or excluded. Segments are clipped to
#' the regions measured, so each base counts at most once.
#'
#' A segment counts as altered when its value differs from the sample's
#' baseline by at least `cutoff`. With `centering = "none"` (default) the
#' baseline is diploid (log.ratio 0 or CN 2). With "median" or "mode", each
#' sample is first centred on its own baseline, exactly as in
#' [segmented_data_to_cnv_status], so that a sample whose whole genome is
#' shifted (e.g. by ploidy or normalisation) isn't counted as altered
#' throughout. Use the same centring as the CNV calls when using PGA as a
#' covariate for them (e.g. the burden in [pairwise_event_association]).
#'
#' @param this_seg Input data frame of seg file.
#' @param seg_path Optionally, specify the path to a local seg file.
#' @param projection Argument specifying the projection of seg file, which will
#'      determine chr prefix, chromosome coordinates, and genome size. Default
#'      is grch37, but hg38 is also accepted.
#' @param cutoff The minimum difference from the sample's baseline for a
#'      segment to be considered as CNV, on the scale of `value_scale`.
#'      Default NULL: 0.56 for "log_ratio" (about 1 copy) and 0.5 for "CN".
#'      This value is expected to be a positive number and applies to both
#'      deletions and amplifications.
#' @param exclude_sex Boolean argument specifying whether to exclude sex
#'      chromosomes from calculation (and from the baseline). Default is TRUE.
#' @param exclude_centromeres Boolean argument specifying whether to exclude
#'      centromeres from calculation: TRUE (default) measures the chromosome
#'      arms only, FALSE measures each chromosome from the start of its p arm
#'      to the end of its q arm.
#' @param per_chromosome Boolean argument definint whether the PGA should be
#'      calculated per total genome size (default) or per each individual
#'      chromosome.
#' @param value_scale "log_ratio" (default; uses the log.ratio column) or "CN"
#'      (uses the CN column).
#' @param centering Baseline each sample is centred on before calling
#'      altered segments: "none" (default; diploid), "median" (the sample's
#'      length-weighted median, as in [segmented_data_to_cnv_status]) or
#'      "mode" (the dominant CN level from [segmented_data_to_ploidy_estimate];
#'      requires `value_scale = "CN"`).
#'
#' @return A data frame with sample_id, PGA (the fraction of the genome
#'      altered), PGA_gain and PGA_loss (the fractions gained and lost,
#'      relative to the sample's baseline). With `per_chromosome = TRUE`, the
#'      same per chromosome: chr<N>_pga, then chr<N>_pga_gain and
#'      chr<N>_pga_loss.
#'
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' seg <- get_sample_cn_segments(
#'     these_sample_ids = c(
#'         "14-36022T",
#'         "DOHH-2"
#'     )
#' )
#'
#' total_pga <- calculate_pga(this_seg = seg)
#' per_chromosome_pga <- calculate_pga(
#'      this_seg = seg,
#'      per_chromosome = TRUE
#' )
#'
#' # centred like segmented_data_to_cnv_status, e.g. as a burden covariate
#' centred_pga <- calculate_pga(this_seg = seg, centering = "median",
#'                              cutoff = 0.1)
#' }
#'
calculate_pga = function(this_seg,
                         seg_path,
                         projection = "grch37",
                         cutoff = NULL,
                         exclude_sex = TRUE,
                         exclude_centromeres = TRUE,
                         per_chromosome = FALSE,
                         value_scale = "log_ratio",
                         centering = "none") {
  # check for required argument
  if (missing(this_seg) & missing (seg_path)) {
    stop("Please provide the data frame of seg file or path to the local seg.")
  }
  if (!value_scale %in% c("log_ratio", "CN")) {
    stop("value_scale must be \"log_ratio\" or \"CN\"")
  }
  if (is.null(cutoff)) {
    cutoff = ifelse(value_scale == "log_ratio", 0.56, 0.5)
  }
  value_col = ifelse(value_scale == "log_ratio", "log.ratio", "CN")

  # ensure the specified projection is correct and define chromosome coordinates
  if (projection == "grch37") {
    chr_coordinates = GAMBLR.data::chromosome_arms_grch37
  } else if (projection == "hg38") {
    chr_coordinates = GAMBLR.data::chromosome_arms_hg38
  } else {
    stop(
      "You specified projection that is currently not supported. Please provide seg files in either hg38 or grch37."
    )
  }

  # exclude sex chromosomes
  if (exclude_sex) {
    chr_coordinates = chr_coordinates %>%
      dplyr::filter(!grepl("X|Y", chromosome))
  }

  # with centromeres included, measure each chromosome from the start of its
  # p arm to the end of its q arm; otherwise measure the arms only
  if (!exclude_centromeres) {
    chr_coordinates = chr_coordinates %>%
      group_by(chromosome) %>%
      mutate(start = min(start),
             end = max(end)) %>%
      ungroup %>%
      distinct(chromosome, start, end)
  }

  # total size of genome in this projection
  genome_size = chr_coordinates %>%
    mutate(size = end - start) %>%
    summarise(genome_size = sum(size)) %>%
    pull(genome_size)

  # prepare for the overlaps
  chr_coordinates = chr_coordinates  %>%
    dplyr::rename("arm_start" = "start",
           "arm_end" = "end",
           "chrom" = "chromosome") %>%
    dplyr::select(chrom, arm_start, arm_end)
  chr_coordinates$chrom = gsub("chr", "", chr_coordinates$chrom)
  # measured length of each chromosome, for per-chromosome PGA
  chr_lengths = chr_coordinates %>%
    group_by(chrom) %>%
    summarise(chrom_len = sum(arm_end - arm_start), .groups = "drop")

  # work out the seg file
  if (!missing(seg_path)) {
    message(paste0("Reading thhe seg file from ", seg_path))
    this_seg = suppressMessages(read_tsv(seg_path))
  }
  this_seg = as.data.frame(this_seg)
  if (!"sample" %in% colnames(this_seg) && "ID" %in% colnames(this_seg)) {
    this_seg = dplyr::rename(this_seg, sample = ID)
  }
  required_cols = c("sample", "chrom", "start", "end", value_col)
  if (!all(required_cols %in% colnames(this_seg))) {
    stop(paste("this_seg must have columns:",
               paste(required_cols, collapse = ", "),
               "(ID is accepted in place of sample)"))
  }

  # preserve the sample ids to account later for those with 0 PGA
  sample_set = this_seg %>% distinct(sample)

  # chromosome names are compared without any prefix
  this_seg$chrom = gsub("chr", "", this_seg$chrom)

  # exclude sex chromosomes
  if (exclude_sex) {
    this_seg = this_seg %>%
      dplyr::filter(!grepl("X|Y", chrom))
  }
  this_seg = this_seg[!is.na(this_seg[[value_col]]), , drop = FALSE]

  # centre each sample on its baseline, as segmented_data_to_cnv_status does
  centre = .cnv_sample_centre(
    data.frame(ID = this_seg$sample, chrom = this_seg$chrom,
               start = this_seg$start, end = this_seg$end),
    this_seg[[value_col]], this_seg$end - this_seg$start,
    value_scale, centering, autosomes_only = exclude_sex
  )
  deviation = this_seg[[value_col]] - centre
  this_seg$direction = ifelse(deviation > 0, "gain", "loss")
  this_seg = this_seg[abs(deviation) >= cutoff, , drop = FALSE]

  # prepare for the overlaps
  this_seg = inner_join(
    this_seg %>% dplyr::select(sample, chrom, start, end, direction),
    chr_coordinates,
    by = "chrom",
    relationship = "many-to-many"
  )

  # what are the segments that overlap? count only the overlapping part
  affected_regions = this_seg %>%
    dplyr::filter(start <= arm_end & arm_start <= end) %>%
    dplyr::mutate(size = pmin(end, arm_end) - pmax(start, arm_start)) %>%
    arrange(sample, chrom, start)

  # altered length per sample (and chromosome) for gains, losses and both
  by_direction = function(x, keys){
    x %>%
      group_by(across(all_of(c(keys, "direction")))) %>%
      summarise(total = sum(size), .groups = "drop") %>%
      pivot_wider(names_from = direction, values_from = total,
                  values_fill = 0) %>%
      dplyr::mutate(gain = if ("gain" %in% names(.)) gain else 0,
                    loss = if ("loss" %in% names(.)) loss else 0,
                    both = gain + loss)
  }

  if(per_chromosome){
   affected_regions = by_direction(affected_regions, c("sample", "chrom")) %>%
    left_join(chr_lengths, by = "chrom")
   chrom_order = unique(chr_coordinates$chrom)
   wide = list()
   for (d in c("both", "gain", "loss")) {
     suffix = c(both = "_pga", gain = "_pga_gain", loss = "_pga_loss")[[d]]
     w = affected_regions %>%
       dplyr::mutate(value = round(.data[[d]] / chrom_len, 3),
                     chrom = factor(chrom, levels = chrom_order)) %>%
       dplyr::select(sample, chrom, value) %>%
       pivot_wider(names_from = chrom, values_from = value,
                   names_expand = TRUE, values_fill = 0)
     colnames(w)[-1] = paste0("chr", colnames(w)[-1], suffix)
     wide[[d]] = w
   }
   affected_regions = wide$both %>%
     left_join(wide$gain, by = "sample") %>%
     left_join(wide$loss, by = "sample")
   # now add any samples that can have 0 PGA
   affected_regions = base::merge(sample_set,
                                 affected_regions,
                                 all.x = TRUE)
   affected_regions[is.na(affected_regions)] <- 0
   affected_regions = affected_regions %>%
    rename("sample_id" = "sample") %>%
    as.data.frame()
  } else {
   affected_regions = by_direction(affected_regions, "sample") %>%
     dplyr::transmute(sample,
                      PGA = round(both / genome_size, 3),
                      PGA_gain = round(gain / genome_size, 3),
                      PGA_loss = round(loss / genome_size, 3))
   # now add any samples that can have 0 PGA
   affected_regions = base::merge(sample_set,
                                 affected_regions,
                                 all.x = TRUE) %>%
    replace_na(list(PGA = 0, PGA_gain = 0, PGA_loss = 0)) %>%
    rename("sample_id" = "sample") %>%
    as.data.frame()
  }

  return(affected_regions)
}
