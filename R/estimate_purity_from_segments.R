#' @title Estimate tumour purity from absolute copy number segments
#'
#' @description Recover the purity behind a pipeline's absolute copy number
#' from the relationship between its observed log ratios and its copy number,
#' e.g. for PureCN, whose segments keep the uncorrected log ratio next to the
#' purity- and ploidy-corrected CN.
#'
#' @details For a segment with copy number C in a tumour of purity p, the
#' observed ratio is `2^log.ratio = k * (p * C + 2 * (1 - p))`, where k
#' absorbs the ploidy and the centring of the log ratios. A length-weighted
#' regression of `2^log.ratio` on C over a sample's autosomal segments (C up
#' to `max_cn`) gives slope `k * p` and intercept `2 * k * (1 - p)`, so
#' `p = 2 * slope / (intercept + 2 * slope)`. On PureCN segments this
#' recovers PureCN's own purity (median error about 0.02); it inherits that
#' fit's errors.
#'
#' The estimate is flagged (`flag`, otherwise "ok") when it is unreliable:
#' \itemize{
#'   \item "one_state": a single copy number state, so no slope.
#'   \item "few_altered": less than `min_altered` of the genome away from CN
#'     2; many noisy CN 2 segments then flatten the slope.
#'   \item "shifted": less than 10% of the genome at CN 2 (e.g. a near-
#'     triploid or near-tetraploid fit).
#'   \item "derived_log_ratio": an almost perfect fit with purity near 1, as
#'     when the log ratio was computed from CN (e.g. log2(CN / 2)), so it
#'     carries no purity information.
#'   \item "out_of_range": purity outside (0, 1].
#' }
#'
#' @param seg_data Segments with ID, chrom, start, end, CN and log.ratio, of
#' one or more samples (and pipelines, see `by`).
#' @param by Optional column to estimate separately for, e.g. "pipeline".
#' @param max_cn Copy number states above this are left out. Default 6.
#' @param min_altered Fraction of the genome that must be away from CN 2 for
#' the estimate not to be flagged "few_altered". Default 0.05.
#'
#' @return A data frame with one row per sample (and `by` value): sample_id,
#' purity, slope, intercept, r_squared, n_states, altered_fraction,
#' cn2_fraction and flag.
#'
#' @export
#'
#' @examples
#' segs = data.frame(ID = "s1", chrom = as.character(1:4), start = 1, end = 5e7,
#'                   CN = c(2, 1, 3, 2), log.ratio = log2(c(1, 0.7, 1.3, 1)))
#' estimate_purity_from_segments(segs)
estimate_purity_from_segments = function(seg_data, by = NULL, max_cn = 6, min_altered = 0.05){
  segs = as.data.frame(seg_data)
  missing_cols = setdiff(c("ID", "chrom", "start", "end", "CN", "log.ratio", by), colnames(segs))
  if(length(missing_cols)) stop("seg_data lacks: ", paste(missing_cols, collapse = ", "))
  segs = segs[.real_segments(segs) & .is_autosome(gsub("^chr", "", as.character(segs$chrom))) &
                !is.na(segs$CN) & !is.na(segs$log.ratio), , drop = FALSE]
  segs$len = as.numeric(segs$end) - as.numeric(segs$start)
  segs = segs[segs$len > 0, , drop = FALSE]
  key = if(is.null(by)) segs$ID else paste(segs$ID, segs[[by]], sep = "\r")
  out = do.call(rbind, lapply(split(seq_len(nrow(segs)), key), function(rows){
    d = segs[rows, , drop = FALSE]
    fit = .purity_fit(2^d$log.ratio, d$CN, d$len, max_cn, min_altered)
    row = data.frame(sample_id = d$ID[1], fit, stringsAsFactors = FALSE)
    if(!is.null(by)) row = cbind(row[, 1, drop = FALSE], stats::setNames(data.frame(d[[by]][1]), by), row[, -1])
    row
  }))
  rownames(out) = NULL
  out
}

#' @title Convert relative copy number to absolute with a reference segmentation
#'
#' @description Calibrate a relative segmentation (e.g. gatk_cnv or CNVkit,
#' whose CN is only a rescaled log ratio) to absolute copy number, sample by
#' sample, with an absolute segmentation of the same sample (e.g. PureCN's),
#' so copy number features can be called on CN with the same rules as for
#' purity-corrected pipelines.
#'
#' @details For each sample, the relative segments are intersected with the
#' reference segments, and `2^log.ratio` of the relative segmentation is
#' regressed on the reference CN over the overlaps (autosomes, length-
#' weighted, CN up to `max_cn`), as in [estimate_purity_from_segments]. Each
#' relative segment's CN is then `(2^log.ratio - intercept) / slope`
#' (at least 0). This absorbs the relative pipeline's own centring and
#' compression as well as purity. The fit's implied purity is reported too.
#'
#' Samples whose fit is flagged (see [estimate_purity_from_segments]), or
#' with no reference segments, are not calibrated: their segments keep their
#' CN and `cn_type` "relative". The exception is "shifted" (little of the
#' genome at CN 2, e.g. a near-tetraploid reference): the fit is linear in the
#' reference's CN whatever its ploidy, so these are calibrated, with the flag
#' kept as a warning; a wrong whole-genome-doubled reference fit is copied. Only use a reference whose fit you trust,
#' e.g. PureCN segmentations that passed [compare_segmentations].
#'
#' @param seg_data Relative segments (ID, chrom, start, end, log.ratio) of
#' one or more samples, from one pipeline.
#' @param reference Absolute segments (ID, chrom, start, end, CN) of the same
#' samples, from one pipeline.
#' @param max_cn Reference copy number states above this are left out of the
#' fit. Default 6.
#' @param min_altered As for [estimate_purity_from_segments]. Default 0.05.
#' @param digits Round the calibrated CN to this many digits. Default 2.
#'
#' @return `seg_data` with CN replaced by the calibrated copy number for the
#' calibrated samples, the original in `CN_relative`, and `cn_type`
#' ("absolute" when calibrated, else "relative"). The fit per sample
#' (sample_id, calibrated, purity, slope, intercept, r_squared, n_states,
#' altered_fraction, cn2_fraction, flag) is in `attr(result,
#' "calibration")`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' segs = GAMBLR.results::get_cn_segments(meta, pipeline = "all")
#' gatk = calibrate_relative_cn(dplyr::filter(segs, pipeline == "gatk_cnv"),
#'                              reference = dplyr::filter(segs, pipeline == "purecn"))
#' attr(gatk, "calibration")
#' }
calibrate_relative_cn = function(seg_data, reference, max_cn = 6, min_altered = 0.05, digits = 2){
  genome_build = attr(seg_data, "genome_build")
  segs = as.data.frame(seg_data)
  ref = as.data.frame(reference)
  missing_cols = setdiff(c("ID", "chrom", "start", "end", "log.ratio"), colnames(segs))
  if(length(missing_cols)) stop("seg_data lacks: ", paste(missing_cols, collapse = ", "))
  missing_cols = setdiff(c("ID", "chrom", "start", "end", "CN"), colnames(ref))
  if(length(missing_cols)) stop("reference lacks: ", paste(missing_cols, collapse = ", "))
  norm_chrom = function(x) gsub("^chr", "", as.character(x))
  ref = ref[.real_segments(ref) & !is.na(ref$CN) & .is_autosome(norm_chrom(ref$chrom)), , drop = FALSE]
  if(!"CN" %in% colnames(segs)) segs$CN = NA_real_
  segs$CN_relative = segs$CN
  segs$cn_type = rep("relative", nrow(segs))

  ids = unique(segs$ID)
  fits = do.call(rbind, lapply(ids, function(id){
    s = segs[segs$ID == id & .real_segments(segs) & !is.na(segs$log.ratio) &
               .is_autosome(norm_chrom(segs$chrom)), c("chrom", "start", "end", "log.ratio"), drop = FALSE]
    r = ref[ref$ID == id, c("chrom", "start", "end", "CN"), drop = FALSE]
    empty = data.frame(sample_id = id, calibrated = FALSE, purity = NA_real_, slope = NA_real_,
                       intercept = NA_real_, r_squared = NA_real_, n_states = 0L, altered_fraction = NA_real_,
                       cn2_fraction = NA_real_, flag = "no_reference", stringsAsFactors = FALSE)
    if(nrow(s) == 0 || nrow(r) == 0) return(empty)
    s$chrom = norm_chrom(s$chrom); r$chrom = norm_chrom(r$chrom)
    pieces = merge(s, r, by = "chrom", suffixes = c("", "_ref"))
    pieces$len = pmin(pieces$end, pieces$end_ref) - pmax(pieces$start, pieces$start_ref)
    pieces = pieces[pieces$len > 0, , drop = FALSE]
    if(nrow(pieces) == 0) return(empty)
    fit = .purity_fit(2^pieces$log.ratio, pieces$CN, pieces$len, max_cn, min_altered)
    # "shifted" (little of the genome at CN 2, e.g. a near-tetraploid
    # reference) doesn't affect the fit, which is linear in the reference CN
    data.frame(sample_id = id, calibrated = fit$flag %in% c("ok", "shifted"), fit, stringsAsFactors = FALSE)
  }))
  rownames(fits) = NULL

  ok = fits[fits$calibrated, , drop = FALSE]
  rows = segs$ID %in% ok$sample_id
  at = match(segs$ID[rows], ok$sample_id)
  segs$CN[rows] = round(pmax(0, (2^segs$log.ratio[rows] - ok$intercept[at]) / ok$slope[at]), digits)
  segs$cn_type[rows] = "absolute"
  n_skipped = sum(!fits$calibrated)
  if(n_skipped > 0){
    message(n_skipped, " of ", nrow(fits), " sample(s) not calibrated (see attr(, \"calibration\")$flag): ",
            paste(names(table(fits$flag[!fits$calibrated])), table(fits$flag[!fits$calibrated]), collapse = ", "))
  }
  if(!is.null(genome_build) && genome_build %in% c("grch37", "hg38")) segs = create_seg_data(segs, genome_build)
  attr(segs, "calibration") = fits
  segs
}

# Length-weighted fit of ratio = k * (p * CN + 2 * (1 - p)) over pieces of a
# sample's genome, with the flags of estimate_purity_from_segments
.purity_fit = function(ratio, cn, len, max_cn, min_altered){
  keep = !is.na(ratio) & !is.na(cn) & cn <= max_cn & len > 0
  ratio = ratio[keep]; cn = cn[keep]; len = len[keep]
  total = sum(len)
  altered = if(total > 0) sum(len[round(cn) != 2]) / total else NA_real_
  cn2 = if(total > 0) sum(len[round(cn) == 2]) / total else NA_real_
  n_states = length(unique(round(cn)))
  out = data.frame(purity = NA_real_, slope = NA_real_, intercept = NA_real_, r_squared = NA_real_,
                   n_states = n_states, altered_fraction = round(altered, 4), cn2_fraction = round(cn2, 4),
                   flag = "one_state", stringsAsFactors = FALSE)
  if(n_states < 2) return(out)
  fit = stats::lm(ratio ~ cn, weights = len)
  slope = unname(stats::coef(fit)[2]); intercept = unname(stats::coef(fit)[1])
  r2 = summary(fit)$r.squared
  purity = 2 * slope / (intercept + 2 * slope)
  out$purity = round(purity, 4); out$slope = slope; out$intercept = intercept; out$r_squared = round(r2, 4)
  out$flag = if(!is.finite(purity) || purity <= 0 || purity > 1.02) "out_of_range" else
    if(r2 > 0.999 && purity > 0.98) "derived_log_ratio" else
      if(altered < min_altered) "few_altered" else
        if(cn2 < 0.1) "shifted" else "ok"
  out
}
