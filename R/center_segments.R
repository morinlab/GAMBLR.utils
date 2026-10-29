#' @title Centre copy number segments on each sample's baseline
#'
#' @description Shift each sample's copy number so that its baseline
#' (neutral) copy number becomes 2, e.g. so that a tetraploid
#' sample's unchanged genome isn't shown as a gain, or to export centred
#' data for viewing in IGV with [seg_data_to_igv].
#'
#' @details The centre of each sample is its baseline copy number from
#' `method` (see [add_segment_baseline]):
#' \itemize{
#'   \item "mode" (default): the mean copy number within the sample's
#'     dominant copy number state. This also removes fractional offsets (a
#'     sample whose neutral level sits at 2.3 is centred at 2.3), as
#'     `adjust_for_ploidy = TRUE` does in [segmented_data_to_cn_matrix].
#'   \item "mode_rounded": the dominant state as a whole number, so only
#'     whole-copy shifts are made.
#'   \item "mean": the genome-wide mean copy number.
#' }
#' A `baseline_cn` column of that method already on the segments is used
#' (e.g. from `GAMBLR.results::get_cn_segments(baseline_method = ...)`);
#' otherwise it is computed with [add_segment_baseline], which needs
#' whole-genome segments.
#' Copy number is shifted by whole copies (`CN - centre + 2`, as for ploidy
#' adjustment elsewhere in GAMBLR); the uncentred values are kept in
#' `CN_uncentered` and the centre in `centre_cn`. Placeholder rows
#' (`dummy_segment == 1`) don't count towards the centre.
#'
#' The log ratio is a separate measurement, not a transform of CN: callers
#' report it relative to the sample's own average (a tetraploid sample's
#' neutral segments are already near 0) and it is compressed by tumour
#' purity, so a shift derived from the CN baseline would be wrong. It is
#' left as it is unless `center_log_ratio = TRUE`, which centres it on its
#' own scale: by the length-weighted mean log ratio of the sample's
#' segments in its baseline copy number state (where it is neutral by CN),
#' kept in `log.ratio_centre`, with the uncentred values in
#' `log.ratio_uncentered`. Samples with no segments in that state (e.g.
#' controlfreec, which doesn't report neutral segments) are left as they
#' are.
#'
#' @param seg_data Segments with ID, chrom, start, end and CN (and
#' optionally log.ratio), e.g. from `GAMBLR.results::get_cn_segments`.
#' @param method "mode" (default), "mode_rounded" or "mean". See Details.
#' The centre is estimated from autosomes; chrX/chrY are still shifted.
#' @param center_log_ratio Also centre log.ratio, on its own scale (see
#' Details). Default FALSE: log.ratio is left as it is.
#'
#' @return `seg_data` with a centred CN column, plus `CN_uncentered`,
#' `centre_cn` and the `baseline_cn`/`baseline_method` it was centred on
#' (and, with `center_log_ratio = TRUE`, a centred log.ratio with
#' `log.ratio_uncentered` and `log.ratio_centre`).
#'
#' @export
#'
#' @examples
#' # one segment per autosome, spanning each whole chromosome (a baseline
#' # needs segments covering at least half of the autosomal genome)
#' arms = GAMBLR.data::chromosome_arms_grch37
#' chrom_end = tapply(arms$end, arms$chromosome, max)[as.character(1:22)]
#' seg = data.frame(ID = "s1", chrom = as.character(1:22), start = 1,
#'                  end = as.numeric(chrom_end), CN = c(3, rep(4, 21)))
#' centred = center_segments(seg)
#' head(centred[, c("ID", "chrom", "CN", "CN_uncentered", "centre_cn")])
center_segments = function(seg_data, method = "mode", center_log_ratio = FALSE){
  method = match.arg(method, baseline_methods())
  if(!all(c("ID", "chrom", "start", "end", "CN") %in% colnames(seg_data))){
    stop("seg_data needs columns ID, chrom, start, end and CN (centring uses the absolute ",
         "copy number; CN isn't derived from log.ratio)")
  }
  if("CN_uncentered" %in% colnames(seg_data)){
    stop("seg_data is already centred (it has a CN_uncentered column)")
  }
  seg_data = add_segment_baseline(seg_data, method = method)
  centre = as.numeric(seg_data$baseline_cn)
  # (add_segment_baseline stops on samples without a baseline)
  centre[centre <= 0] = 2
  if(center_log_ratio){
    if(!"log.ratio" %in% colnames(seg_data)) stop("center_log_ratio = TRUE needs a log.ratio column")
    # from the uncentred CN, so the neutral segments are found in its own state
    lr_centre = .log_ratio_centre(seg_data, centre)
  }
  seg_data$CN_uncentered = seg_data$CN
  seg_data$centre_cn = centre
  seg_data$CN = pmax(as.numeric(seg_data$CN) - centre + 2, 0)
  if(center_log_ratio){
    seg_data$log.ratio_uncentered = seg_data$log.ratio
    seg_data$log.ratio_centre = lr_centre
    seg_data$log.ratio = as.numeric(seg_data$log.ratio) - lr_centre
  }
  seg_data
}

# Each segment's log ratio centre: the length-weighted mean log ratio of its
# sample's autosomal segments in the baseline copy number state (the state
# the CN baseline rounds to), i.e. where the sample is neutral by CN. 0 for
# samples with no such segments (e.g. controlfreec, which doesn't report
# neutral segments)
.log_ratio_centre = function(seg_data, centre){
  s = as.data.frame(seg_data)
  lr = as.numeric(s$log.ratio)
  len = as.numeric(s$end) - as.numeric(s$start)
  neutral = .is_autosome(s$chrom) & .real_segments(s) & !is.na(lr) & !is.na(s$CN) &
    round(as.numeric(s$CN)) == round(centre)
  id = as.character(s$ID)
  by_sample = tapply(lr[neutral] * len[neutral], id[neutral], sum) /
    tapply(len[neutral], id[neutral], sum)
  out = as.vector(by_sample[id])
  none = unique(id[is.na(out)])
  if(length(none) > 0){
    message(length(none), " sample(s) have no segments in their baseline state, so their ",
            "log.ratio is left as it is: ", .list_ids(none))
  }
  out[is.na(out)] = 0
  out
}
