#' @title Adjust ploidy for samples with CNV data.
#'
#' @description `adjust_ploidy` returns a seg file with log.ratios adjusted to the overall sample ploidy.
#'
#' @details This function adjusts the ploidy of the sample using the percent of genome altered (PGA). The PGA is calculated internally, but can also be optionally provided as data frame
#' if calculated from other sources. Only the samples above the threshold-provided PGA will have ploidy adjusted. The function can work with either individual or
#' multi-sample seg file. The telomeres are always excluded from calculation, and sex chromosomes can be optionally included or excluded. The supported projections are grch37 and hg38.
#' The chromosome prefix is handled internally per projection and does not need to be consistent.
#'
#' @param this_seg Input data frame of seg file.
#' @param seg_path Optionally, specify the path to a local seg file.
#' @param projection Argument specifying the projection of seg file, which will determine chr prefix and genome size. Default is grch37, but hg38 is also accepted.
#' @param pga If PGA is calculated through other sources, the data frame with columns sample_id and PGA can be provided in this argument.
#' @param pga_cutoff Minimum PGA for the sample to adjust ploidy. Default is 0.05 (5 %).
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is TRUE.
#' @param return_seg Boolean argument specifying whether to return a data frame in seg-consistent format, or a raw data frame with all step-by-step transformations. Default is TRUE.
#'
#' @return A data frame in seg-consistent format with ploidy-adjusted log ratios.
#'
#' @import dplyr tidyr readr
#' @export
#'
#' @examples
#' all_sample_seg = GAMBLR.data::sample_data$grch37$seg
#' all_sample_seg = dplyr::rename(all_sample_seg, "sample" = "ID")
#'
#' sample_seg = dplyr::filter(all_sample_seg, sample == "02-13135T")
#' adjust_ploidy(this_seg = sample_seg)
#'
#' multi_sample_seg = dplyr::filter(all_sample_seg, sample  %in% c("02-13135T", "SU-DHL-4"))
#' adjust_ploidy(this_seg = multi_sample_seg)
#'
adjust_ploidy = function(this_seg,
                         seg_path,
                         projection = "grch37",
                         pga,
                         pga_cutoff = 0.05,
                         exclude_sex = TRUE,
                         return_seg = TRUE) {

  # ensure the specified projection is correct
  # this is only needed to adjust ploidy if the pre-adjusted ploidy is not provided
  if (!projection %in% c("grch37", "hg38") & missing(pga)) {
    stop(
      "You specified projection that is currently not supported. Please provide seg files in either hg38 or grch37."
    )
  }

  # if the seg is a local file, read it in
  if (!missing(seg_path)) {
    message(paste0("Reading thhe seg file from ", seg_path))
    this_seg = suppressMessages(read_tsv(seg_path))
  }

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

  # if PGA is called with custom parameters or obtained elsewhere, use it - but if not, calculate on-the-fly
  if (missing(pga)) {
    message("Calculating PGA ...")
    pga = calculate_pga(this_seg = this_seg,
                        projection = projection,
                        exclude_sex = exclude_sex)
  } else {
    # ensure the column is named "sample_id" if it came outside GAMBLR
    if (!"sample_id" %in% colnames(pga)) {
      stop("Please ensure the column with sample ids in your PGA data frame is named sample_id.")
    }
  }

  # add the PGA information to the seg file
  this_seg = left_join(this_seg,
                       pga,
                       by = c("sample" = "sample_id"))

  # By how much we should adjust the ploidy?
  this_seg = this_seg %>%
    group_by(sample) %>% # account for multi-sample seg files
    dplyr::mutate(
      adjust = mean((2 * 2 ^ log.ratio)),
      # convert log.ratio to absolute CN and find average
      adjust = PGA * adjust,
      # what is the ploidy of genome affected by CNV?
      neutral = (1 - PGA) * 2,
      # how much of the genome is not affected by CNV? Assume it's in diploid state
      adjust = adjust + neutral,
      # overall genome's ploidy status
      adjust = abs(adjust - 2)
    ) %>% # by how much the ploidy should be adjusted? From average sample ploidy take out the diploid state
    dplyr::mutate(need_to_adjust = ifelse((PGA > pga_cutoff &
                                             !log.ratio == 0), "TRUE", "FALSE")) # We only need to adjust ploidy if the PGA is above cut-off and segment is not diploid

  # Adjust ploidy
  this_seg = this_seg %>%
    ungroup() %>%
    mutate(cn = (2*2^log.ratio), .before = log.ratio) %>% # convert log.ratios to absolute CN
    mutate(new_cn = ifelse(need_to_adjust == "TRUE", # always round to the integer, but adjust only if needed
                           round(cn-adjust),
                           round(cn))) %>%
    mutate(new_log.ratio = ifelse(need_to_adjust == "TRUE", # log transform the new CN states
                                  (log(abs(new_cn),2)-1),
                                  log.ratio)) %>%
    mutate(new_log.ratio = ifelse(new_log.ratio == -Inf,-10, new_log.ratio)) # deal with the -Inf for low negative numbers

  # this allows user to see all the transformations if they wish or return the standard seg file
  if (return_seg) {
    message("Returning the seg file with ploidy-adjusted CN ...")
    this_seg = this_seg %>%
      mutate(log.ratio = new_log.ratio) %>% # assign the new log.ratio
      select(sample, chrom, start, end, LOH_flag, log.ratio) # select only columns of the seg file
  }


  return(this_seg)
}
