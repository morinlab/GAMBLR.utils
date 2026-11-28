#' @title QC metrics for segmented copy number data
#'
#' @description Compute per-sample quality metrics for segmented copy number
#' data, for deciding which samples' copy number to trust (e.g. to set
#' `cn_available` to FALSE in the metadata given to
#' [assemble_classifier_features]). Only metrics are returned; there are no
#' thresholds or flags yet.
#'
#' @details Metrics are computed on the autosomes (unless `autosomes_only =
#' FALSE`), on copy number: the CN column, or `2 * 2^log.ratio` with
#' `value_scale = "log_ratio"` or when there is no CN column. Placeholder
#' (dummy) segments are ignored.
#' \itemize{
#'   \item `n_segments`, `covered_mb`, `segments_per_100mb`: segmentation
#'     density, the last normalised by the genome the segments cover.
#'   \item `short_segment_fraction`: fraction of segments shorter than
#'     `short_segment` bp (default 1 Mb).
#'   \item `median_step`: median absolute copy number change between adjacent
#'     segments of a chromosome. Callers differ: over-segmented profiles from
#'     a continuous caller (e.g. cnvkit) have many small steps, while
#'     over-segmented integer copy number (e.g. PureCN) jumps between states,
#'     so compare it within a caller.
#'   \item `baseline_mode`, `mean_cn`: the sample's dominant copy number (as
#'     [segmented_data_to_ploidy_estimate] with method "mode") and its
#'     length-weighted mean copy number.
#'   \item `supplied_baseline`, `baseline_delta`: the `baseline_cn` column of
#'     `seg_data` (e.g. from GAMBLR.results) and its difference from
#'     `baseline_mode`; NA without that column.
#'   \item `altered_fraction`: fraction of the covered genome at least
#'     `altered_change` copies (default 0.5) from `baseline_mode`. High in noisy
#'     profiles, but also in genuinely aneuploid tumours.
#'   \item `deep_loss_fraction`: fraction of the covered genome below
#'     `deep_loss` copies (default 0.5). Homozygous deletions are small, so a
#'     large fraction means a failed profile.
#'   \item `deep_loss_logratio`: length-weighted median log ratio of those
#'     segments (needs both CN and log.ratio). Near 0, the CN column
#'     contradicts the log ratios (e.g. a failed absolute fit); strongly
#'     negative, the coverage itself is extreme.
#'   \item With `compare_seg_data` (e.g. another caller's segments for the
#'     same samples): `compare_agreement`, the fraction of `bin_size` bins
#'     covered by both where the two agree on loss, neutral or gain, each
#'     relative to its own `baseline_mode` (so a doubled ploidy solution
#'     doesn't count as disagreement); `compare_correlation`, the Pearson
#'     correlation across those bins of the two profiles' copy number, each
#'     minus its own `baseline_mode` (a continuous measure, unaffected by
#'     differences in ploidy or in how purity compresses changes; NA when
#'     either profile is flat); and `compare_bins`, the number of bins
#'     compared.
#' }
#'
#' @param seg_data Segmented data with columns ID, chrom, start, end, and CN
#' and/or log.ratio.
#' @param compare_seg_data Optional second set of segments for the same
#' samples, in the same format, to compare with.
#' @param value_scale "CN" (default when there is a CN column) or "log_ratio".
#' Used for both `seg_data` and `compare_seg_data`.
#' @param autosomes_only Use autosomes (1-22) only. Default TRUE.
#' @param short_segment Length (bp) below which a segment counts as short.
#' Default 1e6.
#' @param deep_loss Copy number below which a segment counts as a deep
#' (homozygous) loss for `deep_loss_fraction`. Default 0.5.
#' @param altered_change Copy number change from the baseline counted as
#' altered (also the loss/gain threshold for `compare_agreement`). Default
#' 0.5.
#' @param bin_size Bin size (bp) for `compare_agreement`. Default 1e6.
#'
#' @return A data frame with one row per sample (sample_id, then the metrics).
#'
#' @export
#'
#' @examples
#' segs = data.frame(ID = "s1", chrom = c("1", "1", "2"),
#'                   start = c(1, 5e7, 1), end = c(5e7, 2.4e8, 2.4e8),
#'                   CN = c(2, 1, 2), log.ratio = c(0, -1, 0))
#' segmented_data_qc(segs)
#' \dontrun{
#' purecn = get_cn_segments(meta, pipeline = "purecn", projection = "grch37")
#' cnvkit = get_cn_segments(meta, pipeline = "cnvkit", projection = "grch37")
#' qc = segmented_data_qc(purecn, compare_seg_data = cnvkit)
#' }
segmented_data_qc = function(seg_data,
                             compare_seg_data = NULL,
                             value_scale = NULL,
                             autosomes_only = TRUE,
                             short_segment = 1e6,
                             deep_loss = 0.5,
                             altered_change = 0.5,
                             bin_size = 1e6){
  segs = .qc_prepare_segments(seg_data, value_scale, autosomes_only)
  value_scale = attr(segs, "value_scale")
  # samples in the order they first appear in seg_data
  samples = intersect(unique(as.character(as.data.frame(seg_data)$ID)), segs$ID)
  qc = do.call(rbind, lapply(split(segs, segs$ID)[samples], function(s){
    len = s$end - s$start
    mode = .qc_mode(s$value, len)
    steps = unlist(lapply(split(s$value, s$chrom), function(v) abs(diff(v))))
    data.frame(
      sample_id = s$ID[1],
      n_segments = nrow(s),
      covered_mb = sum(len) / 1e6,
      segments_per_100mb = 100 * nrow(s) / (sum(len) / 1e6),
      short_segment_fraction = mean(len < short_segment),
      median_step = if(length(steps)) stats::median(steps) else NA_real_,
      baseline_mode = mode,
      mean_cn = sum(s$value * len) / sum(len),
      supplied_baseline = if("baseline_cn" %in% colnames(s)) s$baseline_cn[1] else NA_real_,
      altered_fraction = sum(len[abs(s$value - mode) >= altered_change]) / sum(len),
      deep_loss_fraction = sum(len[s$value < deep_loss]) / sum(len),
      deep_loss_logratio = if("log.ratio" %in% colnames(s) && any(s$value < deep_loss))
        .qc_weighted_median(s$log.ratio[s$value < deep_loss], len[s$value < deep_loss])
        else NA_real_,
      stringsAsFactors = FALSE)
  }))
  qc$baseline_delta = qc$supplied_baseline - qc$baseline_mode
  qc = qc[, c("sample_id", "n_segments", "covered_mb", "segments_per_100mb",
              "short_segment_fraction", "median_step", "baseline_mode",
              "mean_cn", "supplied_baseline", "baseline_delta",
              "altered_fraction", "deep_loss_fraction", "deep_loss_logratio")]

  if(!is.null(compare_seg_data)){
    other = .qc_prepare_segments(compare_seg_data, value_scale, autosomes_only)
    agree = do.call(rbind, lapply(samples, function(id){
      .qc_bin_agreement(segs[segs$ID == id, ], other[other$ID == id, ],
                        bin_size, altered_change)
    }))
    qc$compare_agreement = agree$agreement
    qc$compare_correlation = agree$correlation
    qc$compare_bins = agree$bins
  }
  rownames(qc) = NULL
  qc
}

# segments with a `value` column in copies, sorted, without dummy rows
.qc_prepare_segments = function(seg_data, value_scale, autosomes_only){
  segs = as.data.frame(seg_data)
  if(is.null(value_scale)){
    value_scale = if("CN" %in% colnames(segs)) "CN" else "log_ratio"
  }
  if(!value_scale %in% c("CN", "log_ratio")){
    stop("value_scale must be \"CN\" or \"log_ratio\"")
  }
  value_col = if(value_scale == "CN") "CN" else "log.ratio"
  required = c("ID", "chrom", "start", "end", value_col)
  if(!all(required %in% colnames(segs))){
    stop("seg_data must have columns: ", paste(required, collapse = ", "))
  }
  segs = segs[.real_segments(segs), , drop = FALSE]
  segs$chrom = sub("^chr", "", as.character(segs$chrom))
  if(autosomes_only) segs = segs[.is_autosome(segs$chrom), , drop = FALSE]
  # numeric, since genome-wide sums of integer lengths overflow
  segs$start = as.numeric(segs$start)
  segs$end = as.numeric(segs$end)
  segs$value = if(value_scale == "CN") as.numeric(segs$CN) else 2 * 2^as.numeric(segs$log.ratio)
  segs = segs[!is.na(segs$value) & segs$end > segs$start, , drop = FALSE]
  segs = segs[order(segs$ID, segs$chrom, segs$start), , drop = FALSE]
  attr(segs, "value_scale") = value_scale
  segs
}

.qc_weighted_median = function(x, w){
  ok = !is.na(x)
  x = x[ok]
  w = w[ok]
  if(length(x) == 0) return(NA_real_)
  o = order(x)
  x[o][which(cumsum(w[o]) >= sum(w) / 2)[1]]
}

# length-weighted mode of copy number, on values rounded to 0.1 copies
.qc_mode = function(value, len){
  w = tapply(len, round(value, 1), sum)
  as.numeric(names(w)[which.max(w)])
}

# agreement on loss / neutral / gain and correlation of copy number (each
# relative to its profile's own mode) in fixed bins, using the segment
# covering each bin's midpoint
.qc_bin_agreement = function(a, b, bin_size, change){
  if(nrow(a) == 0 || nrow(b) == 0){
    return(data.frame(agreement = NA_real_, correlation = NA_real_, bins = 0L))
  }
  mode_a = .qc_mode(a$value, a$end - a$start)
  mode_b = .qc_mode(b$value, b$end - b$start)
  state = function(s, mode, chrom, mid){
    s = s[s$chrom == chrom, , drop = FALSE]
    i = findInterval(mid, s$start)
    v = rep(NA_real_, length(mid))
    ok = i > 0
    ok[ok] = mid[ok] <= s$end[i[ok]]
    v[ok] = s$value[i[ok]]
    v - mode
  }
  calls = do.call(rbind, lapply(intersect(unique(a$chrom), unique(b$chrom)), function(chrom){
    top = max(a$end[a$chrom == chrom], b$end[b$chrom == chrom])
    mid = seq(bin_size / 2, top, by = bin_size)
    data.frame(a = state(a, mode_a, chrom, mid), b = state(b, mode_b, chrom, mid))
  }))
  both = !is.na(calls$a) & !is.na(calls$b)
  call = function(d) ifelse(d <= -change, -1, ifelse(d >= change, 1, 0))
  a = calls$a[both]
  b = calls$b[both]
  flat = length(a) < 3 || stats::sd(a) == 0 || stats::sd(b) == 0
  data.frame(agreement = if(any(both)) mean(call(a) == call(b)) else NA_real_,
             correlation = if(flat) NA_real_ else stats::cor(a, b),
             bins = sum(both))
}
