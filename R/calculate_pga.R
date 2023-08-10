#' @title Calculate proportion of genome altered by CNV.
#'
#' @description `calculate_pga` returns a data.frame with estimated proportion of genome altered for each sample.
#'
#' @details This function calculates the percent of genome altered (PGA) by CNV. It takes into account the total length of
#' sample's CNV and relates it to the total genome length to return the proportion affected by CNV. The input is expected to be a seg file.
#' The path to a local SEG file can be provided instead. If The custom seg file is provided, the minimum required columns are
#' sample, chrom, start, end, and log.ratio. The function can work with either individual or multi-sample seg files. The telomeres are always
#' excluded from calculation, and centromeres/sex chromosomes can be optionally included or excluded.
#'
#' @param this_seg Input data frame of seg file.
#' @param seg_path Optionally, specify the path to a local seg file.
#' @param projection Argument specifying the projection of seg file, which will determine chr prefix, chromosome coordinates, and genome size. Default is grch37, but hg38 is also accepted.
#' @param cutoff The minimum log.ratio for the segment to be considered as CNV. Default is 0.56, which is 1 copy. This value is expected to be a positive float of log.ratio for both deletions and amplifications.
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is TRUE.
#' @param exclude_centromeres Boolean argument specifying whether to exclude centromeres from calculation. Default is TRUE.
#'
#' @return data frame
#'
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' sample_seg = get_sample_cn_segments(this_sample_id = "14-36022T")
#' sample_seg = dplyr::rename(sample_seg, "sample" = "ID")
#'
#' calculate_pga(this_seg = sample_seg)
#'
#' calculate_pga(this_seg = sample_seg,
#'               exclude_sex = FALSE)
#'
#' one_sample = get_sample_cn_segments(this_sample_id = "14-36022T")
#' one_sample = dplyr::rename(one_sample, "sample" = "ID")
#'
#' another_sample = get_sample_cn_segments(this_sample_id = "BLGSP-71-21-00243-01A-11E")
#' another_sample = dplyr::rename(another_sample, "sample" = "ID")
#'
#' multi_sample_seg = rbind(one_sample, another_sample)
#'
#' calculate_pga(this_seg = multi_sample_seg)
#'
calculate_pga = function(this_seg,
                         seg_path,
                         projection = "grch37",
                         cutoff = 0.56,
                         exclude_sex = TRUE,
                         exclude_centromeres = TRUE) {
  # check for required argument
  if (missing(this_seg) & missing (seg_path)) {
    stop("Please provide the data frame of seg file or path to the local seg.")
  }

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

  # does the user's seg file contain centromeres?
  if (exclude_centromeres) {
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
    rename("arm_start" = "start",
           "arm_end" = "end",
           "chrom" = "chromosome")

  # work out the seg file
  if (!missing(seg_path)) {
    message(paste0("Reading thhe seg file from ", seg_path))
    this_seg = suppressMessages(read_tsv(seg_path))
  }

  # preserve the sample ids to account later for those with 0 PGA
  sample_set = this_seg %>% distinct(sample)

  this_seg = this_seg %>%
    dplyr::filter(abs(log.ratio) >= cutoff) %>%
    dplyr::relocate(sample, .after = last_col())

  # ensure consistent chromosome prefixing
  if (projection == "grch37") {
    this_seg$chrom = gsub("chr", "", this_seg$chrom)
  } else {
    this_seg$chrom = gsub("chr", "", this_seg$chrom) # if there is a mish-mash of prefixes, strip them all
    this_seg$chrom = paste0("chr", this_seg$chrom)
  }

  # exclude sex chromosomes
  if (exclude_sex) {
    this_seg = this_seg %>%
      dplyr::filter(!grepl("X|Y", chrom))
  }

  # prepare for the overlaps
  this_seg = inner_join(
    this_seg,
    chr_coordinates,
    by = "chrom",
    relationship = "many-to-many"
  )

  # what are the segments that overlap?
  this_seg = this_seg %>%
    dplyr::filter(start <= arm_end & arm_start <= end) %>%
    arrange(sample, chrom, start)

  # calculate total length of CNV
  affected_regions = this_seg %>%
    dplyr::mutate(size = end - start) %>%
    group_by(sample) %>%
    summarise(total = sum(size))

  affected_regions$PGA = affected_regions$total / genome_size

  # now add any samples that can have 0 PGA
  affected_regions = base::merge(sample_set,
                                 affected_regions,
                                 all.x = TRUE) %>%
    replace_na(list(PGA = 0))

  affected_regions = affected_regions %>%
    select(-total) %>%
    `names<-`(c("sample_id", "PGA")) %>%
    as.data.frame() %>%
    dplyr::mutate(PGA = round(PGA, 3))

  return(affected_regions)

}
