#' @title Estimate per-sample ploidy from segmented copy number data
#'
#' @description Estimate a baseline ploidy for each sample in a set of copy
#' number segments. This is the baseline used by the `adjust_for_ploidy`
#' option of [segmented_data_to_cn_matrix] and
#' `GAMBLR.viz::segments_to_exact_cn_matrix()`, exposed so it can also be
#' used directly, e.g. as a numeric metadata column for plotting.
#'
#' @details Two estimators are available:
#'
#' `"mode"` (default): the dominant CN level -- whichever (rounded) CN level
#' covers the most total genomic length -- is found as a length-weighted
#' mode, not a length-weighted mean. The mean isn't robust to a skewed CN
#' distribution: a sample that's genuinely diploid at baseline but carries
#' a substantial burden of gains (very common) gets its mean pulled above 2
#' by those gains alone, even though most of its genome (by length) is
#' still truly CN=2. Shifting by that inflated mean then pushes the
#' actually-neutral majority of the genome into "loss" territory.
#' The dominant level is *identified* by rounding, but the value returned is
#' the continuous, length-weighted mean CN of the segments inside that
#' dominant level, not the rounded integer label. Using the rounded label
#' would silently skip correction for exactly the samples that most need
#' it: a low-signal sample whose caller-assigned neutral CN is consistently
#' ~1.9 still has its dominant level round to 2, which would give a
#' baseline of 2 and a shift of 0. Ties between levels are broken in favour
#' of the lower CN level.
#'
#' `"mean"`: the length-weighted mean CN across the genome. This is the
#' more conventional definition of ploidy and, unlike `"mode"`, reflects
#' the overall burden of gains and losses, so it is the better choice when
#' ploidy itself is the quantity of interest (e.g. testing whether specific
#' CNVs are associated with higher ploidy).
#'
#' @param seg_data A data frame of segments with columns ID, chrom, start,
#' end and CN.
#' @param method Either "mode" (default) or "mean". See Details.
#' @param autosomes_only Exclude chrX/chrY segments from the calculation.
#' Default TRUE: sex chromosome copy number reflects patient sex, not tumour
#' biology, and would bias the estimate for roughly half of a mixed-sex
#' cohort.
#' @param max_CN_allowed CN values higher than this are set to this value
#' before estimating. Default 6.
#'
#' @return A data frame with columns sample_id and ploidy, one row per
#' sample with at least one segment remaining after the autosome filter.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' seg_data = data.frame(
#'   ID = c("s1", "s1", "s1", "s2", "s2"),
#'   chrom = c("1", "2", "X", "1", "2"),
#'   start = c(1, 1, 1, 1, 1),
#'   end = c(2e8, 1e8, 1e8, 2e8, 1e8),
#'   CN = c(2, 3, 1, 4, 4)
#' )
#' segmented_data_to_ploidy_estimate(seg_data)
#' segmented_data_to_ploidy_estimate(seg_data, method = "mean")
segmented_data_to_ploidy_estimate = function(seg_data,
                                             method = "mode",
                                             autosomes_only = TRUE,
                                             max_CN_allowed = 6){
  required_cols = c("ID", "chrom", "start", "end", "CN")
  if(!all(required_cols %in% colnames(seg_data))){
    stop(paste(
      "seg_data must have columns:", paste(required_cols, collapse = ", ")
    ))
  }
  if(!method %in% c("mode", "mean")){
    stop("method must be 'mode' or 'mean'")
  }
  seg_data = as.data.frame(seg_data) %>%
    mutate(
      CN = ifelse(CN > max_CN_allowed, max_CN_allowed, CN),
      .size = end - start
    )
  if(autosomes_only){
    seg_data = filter(
      seg_data, !grepl("^(chr)?[XY]$", chrom, ignore.case = TRUE)
    )
  }
  if(method == "mean"){
    ploidy_df = seg_data %>%
      group_by(ID) %>%
      summarise(ploidy = stats::weighted.mean(CN, .size), .groups = "drop")
  }else{
    ploidy_df = seg_data %>%
      mutate(CN_level = round(CN)) %>%
      group_by(ID, CN_level) %>%
      summarise(
        total_length = sum(.size),
        level_CN = stats::weighted.mean(CN, .size),
        .groups = "drop"
      ) %>%
      arrange(CN_level) %>%
      group_by(ID) %>%
      slice_max(total_length, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(ID, ploidy = level_CN)
  }
  ploidy_df %>%
    select(sample_id = ID, ploidy) %>%
    as.data.frame()
}
