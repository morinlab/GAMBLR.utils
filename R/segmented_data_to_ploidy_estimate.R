#' @title Estimate per-sample ploidy from segmented copy number data
#'
#' @description Estimate a baseline ploidy for each sample in a set of copy
#' number segments. This is the baseline used by the `adjust_for_ploidy`
#' option of [segmented_data_to_cn_matrix] and
#' `GAMBLR.viz::segments_to_exact_cn_matrix()`, exposed so it can also be
#' used directly, e.g. as a numeric metadata column for plotting.
#'
#' @details The methods are those of [add_segment_baseline] (see
#' [baseline_methods]); "mode_rounded" gives the dominant state as a whole
#' number. Samples whose segments all come from controlfreec get 2 (see
#' [add_segment_baseline]). When the segments already have a `baseline_cn` column of the
#' method (e.g. region segments from `GAMBLR.results::get_cn_segments`,
#' whose baselines are genome-wide), those values are returned instead of
#' an estimate from the segments; a "mode" column also gives
#' "mode_rounded".
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
#' end and CN. Placeholder rows (`dummy_segment == 1`) are left out.
#' @param method "mode" (default), "mode_rounded" or "mean". See Details.
#' @param autosomes_only Use autosomes (1-22) only: sex chromosomes and any
#' other contigs are left out. Default TRUE: sex chromosome copy number reflects patient sex, not tumour
#' biology, and would bias the estimate for roughly half of a mixed-sex
#' cohort.
#' @param max_CN_allowed CN values higher than this are set to this value
#' before estimating. Default 20, as for every baseline in GAMBLR (see
#' [add_segment_baseline]); it only matters for "mean" or a dominant state
#' above the cap.
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
                                             max_CN_allowed = 20){
  required_cols = c("ID", "chrom", "start", "end", "CN")
  if(!all(required_cols %in% colnames(seg_data))){
    stop(paste(
      "seg_data must have columns:", paste(required_cols, collapse = ", ")
    ))
  }
  method = match.arg(method, baseline_methods())
  # baselines already on the segments (e.g. computed genome-wide for region
  # segments by GAMBLR.results::get_cn_segments) are used as they are
  existing = .existing_baseline(seg_data, method)
  if(!is.null(existing)){
    first = !duplicated(as.character(seg_data$ID))
    return(data.frame(sample_id = as.character(seg_data$ID)[first], ploidy = existing[first],
                      stringsAsFactors = FALSE))
  }
  if("baseline_method" %in% colnames(seg_data)){
    message("the segments' baseline_cn is \"",
            paste(unique(stats::na.omit(seg_data$baseline_method)), collapse = ", "),
            "\", so the \"", method, "\" baseline is estimated from the segments instead")
  }
  # placeholder rows never count towards ploidy
  seg_data = as.data.frame(seg_data)
  seg_data = seg_data[.real_segments(seg_data), , drop = FALSE] %>%
    mutate(
      CN = ifelse(CN > max_CN_allowed, max_CN_allowed, CN),
      .size = end - start
    )
  if(autosomes_only){
    seg_data = seg_data[.is_autosome(seg_data$chrom), , drop = FALSE]
  }
  seg_data = seg_data[!is.na(seg_data$CN), , drop = FALSE]
  if(nrow(seg_data) == 0) return(data.frame(sample_id = character(0), ploidy = numeric(0)))
  # the same estimator as add_segment_baseline, including fixed baselines
  # (controlfreec samples are 2)
  ploidy = .baseline_from(as.character(seg_data$ID), seg_data$CN, seg_data$.size, method,
                          max_cn = max_CN_allowed)
  fixed = .fixed_baselines(seg_data)
  ploidy[names(fixed)] = fixed
  data.frame(sample_id = names(ploidy), ploidy = as.numeric(ploidy), stringsAsFactors = FALSE)
}
