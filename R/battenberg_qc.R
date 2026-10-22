#' @title Allele-specific QC for Battenberg summaries
#'
#' @description Annotate the samples in a [summarize_battenberg] result with
#' allele-specific QC flags, a status and the reasons for it, using
#' configurable thresholds. The data are not changed or removed; use
#' [battenberg_qc_filter] to restrict analyses and [battenberg_qc_report] to
#' compare QC outcomes between groups.
#'
#' @details The criteria, with the default (provisional) thresholds:
#' \itemize{
#'   \item `qc_fit_severe`: clonal fit score below `severe_fit_score` (0.10).
#'   \item `qc_fit_review`: clonal fit score from `severe_fit_score` up to
#'     `review_fit_score` (0.50).
#'   \item `qc_high_subclonal`: subclonal fraction above
#'     `max_subclonal_fraction` (0.90).
#'   \item `qc_low_assessable`: allele-specific assessable fraction below
#'     `min_assessable_fraction` (0.70).
#'   \item `qc_ploidy_discordant`: |final ploidy - clonal fit ploidy| above
#'     `max_ploidy_delta` (0.5).
#'   \item `qc_homdel_review`: largest homozygous deletion above
#'     `max_homdel_mb` (10 Mb), or homozygous deletion covering more than
#'     `max_homdel_fraction` (0.01) of the called genome.
#'   \item `qc_homdel_severe`: largest homozygous deletion above
#'     `severe_homdel_mb` (50 Mb).
#'   \item `qc_missing_fit_info`: purity, final ploidy, clonal fit ploidy or
#'     (when it is used) the clonal fit score is missing.
#' }
#' The clonal fit score (FRAC_GENOME distance) is read as the fraction of the
#' genome compatible with clonal states only under Battenberg's default
#' distance metric (`clonality_dist_metric = 0`); for any other or unknown
#' (NA) metric the fit-score criteria are not applied, and the reasons say so.
#'
#' `allelic_qc_status` is `severe_review` if the fit or the homozygous
#' deletion criterion is severe, otherwise
#' `unknown` if fitting information is missing, otherwise `review` if any
#' other criterion is met, otherwise `pass`. `allelic_qc_reasons` lists every
#' criterion met (separated by "; "). These thresholds are starting points for
#' inspecting cohort-wide distributions, not established Battenberg
#' standards, and a low fit score can reflect genuine heterogeneity rather
#' than a technical failure: the flags identify cases for review and
#' sensitivity analyses.
#'
#' @param battenberg Output of [summarize_battenberg].
#' @param qc_thresholds Optional named list overriding any of the defaults:
#' `severe_fit_score` (0.10), `review_fit_score` (0.50),
#' `max_subclonal_fraction` (0.90), `min_assessable_fraction` (0.70),
#' `max_ploidy_delta` (0.5), `max_homdel_mb` (10), `max_homdel_fraction`
#' (0.01), `severe_homdel_mb` (50).
#' @param clonality_dist_metric The Battenberg `clonality_dist_metric` used to
#' produce the results: 0 (default, Battenberg's default) or NA if unknown.
#'
#' @return `battenberg` with the QC columns added to (or replaced in)
#' `summary`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bb = summarize_battenberg(battenberg_dir = "battenberg_check")
#' table(bb$summary$allelic_qc_status)
#' # stricter subclonality threshold
#' bb = battenberg_qc(bb, qc_thresholds = list(max_subclonal_fraction = 0.7))
#' }
battenberg_qc = function(battenberg,
                         qc_thresholds = NULL,
                         clonality_dist_metric = 0){
  defaults = list(severe_fit_score = 0.10, review_fit_score = 0.50,
                  max_subclonal_fraction = 0.90, min_assessable_fraction = 0.70,
                  max_ploidy_delta = 0.5, max_homdel_mb = 10,
                  max_homdel_fraction = 0.01, severe_homdel_mb = 50)
  unknown = setdiff(names(qc_thresholds), names(defaults))
  if(length(unknown) > 0){
    stop(paste0("unknown qc_thresholds: ", paste(unknown, collapse = ", "),
                "; valid: ", paste(names(defaults), collapse = ", ")))
  }
  th = utils::modifyList(defaults, as.list(qc_thresholds))
  s = battenberg$summary
  if(is.null(s) || nrow(s) == 0){
    return(battenberg)
  }
  use_score = !is.na(clonality_dist_metric) && clonality_dist_metric == 0
  score = s$clonal_fit_score
  flag = function(cond) !is.na(cond) & cond
  s$qc_fit_severe = if(use_score) flag(score < th$severe_fit_score) else FALSE
  s$qc_fit_review = if(use_score) flag(score >= th$severe_fit_score &
                                         score < th$review_fit_score) else FALSE
  s$qc_high_subclonal = flag(s$subclonal_fraction > th$max_subclonal_fraction)
  s$qc_low_assessable = flag(s$assessable_fraction < th$min_assessable_fraction) |
    is.na(s$assessable_fraction)
  s$qc_ploidy_discordant = flag(abs(s$ploidy_delta_clonal_final) > th$max_ploidy_delta)
  s$qc_homdel_severe = flag(s$largest_homdel_mb > th$severe_homdel_mb)
  s$qc_homdel_review = flag(s$largest_homdel_mb > th$max_homdel_mb) |
    flag(s$homdel_fraction > th$max_homdel_fraction)
  s$qc_missing_fit_info = is.na(s$purity) | is.na(s$final_profile_ploidy) |
    is.na(s$clonal_fit_ploidy) | (use_score & is.na(score))

  labels = c(qc_fit_severe = paste0("clonal fit score < ", th$severe_fit_score),
             qc_fit_review = paste0("clonal fit score < ", th$review_fit_score),
             qc_high_subclonal = paste0("subclonal fraction > ", th$max_subclonal_fraction),
             qc_low_assessable = paste0("assessable fraction < ", th$min_assessable_fraction),
             qc_ploidy_discordant = paste0("|final - clonal fit ploidy| > ", th$max_ploidy_delta),
             qc_homdel_severe = paste0("homozygous deletion > ", th$severe_homdel_mb, " Mb"),
             qc_homdel_review = paste0("homozygous deletion > ", th$max_homdel_mb,
                                       " Mb or > ", th$max_homdel_fraction, " of the genome"),
             qc_missing_fit_info = "missing fitting information")
  s$allelic_qc_reasons = vapply(seq_len(nrow(s)), function(i){
    hits = names(labels)[vapply(names(labels), function(f) isTRUE(s[[f]][i]), logical(1))]
    reasons = unname(labels[hits])
    if(!use_score){
      reasons = c(reasons, "fit score not used (clonality_dist_metric not 0)")
    }
    if(length(reasons) == 0) "" else paste(reasons, collapse = "; ")
  }, character(1))
  s$allelic_qc_status = ifelse(s$qc_fit_severe | s$qc_homdel_severe, "severe_review",
                          ifelse(s$qc_missing_fit_info, "unknown",
                          ifelse(s$qc_fit_review | s$qc_high_subclonal |
                                   s$qc_low_assessable | s$qc_ploidy_discordant |
                                   s$qc_homdel_review,
                                 "review", "pass")))
  battenberg$summary = s
  battenberg
}

#' @title Restrict a Battenberg summary to samples passing allele-specific QC
#'
#' @description Return a copy of a [summarize_battenberg] result (all its
#' tables) limited to samples with the given QC statuses, e.g. only passing
#' samples for a primary analysis, or passing plus review for a sensitivity
#' analysis. The original object is not modified.
#'
#' @param battenberg Output of [summarize_battenberg] (with QC columns, as
#' added by [battenberg_qc]).
#' @param keep QC statuses to keep. Default "pass". Use
#' `c("pass", "review")` for a sensitivity analysis.
#'
#' @return The same list, with every table restricted to those samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' primary = battenberg_qc_filter(bb)
#' sensitivity = battenberg_qc_filter(bb, keep = c("pass", "review"))
#' }
battenberg_qc_filter = function(battenberg, keep = "pass"){
  if(!"allelic_qc_status" %in% colnames(battenberg$summary)){
    stop("no QC columns; run battenberg_qc() first")
  }
  valid = c("pass", "review", "severe_review", "unknown")
  if(!all(keep %in% valid)){
    stop(paste("keep must be among:", paste(valid, collapse = ", ")))
  }
  samples = battenberg$summary$sample[battenberg$summary$allelic_qc_status %in% keep]
  lapply(battenberg, function(tab){
    if(is.data.frame(tab) && "sample" %in% colnames(tab)) tab[tab$sample %in% samples, , drop = FALSE] else tab
  })
}

#' @title Compare Battenberg QC outcomes between groups
#'
#' @description Summarize allele-specific QC by group (e.g. NMF component)
#' and whole-genome doubling status, so that differences in QC failure
#' between groups, which could themselves create apparent differences in LOH
#' or allelic imbalance, are visible before comparing groups.
#'
#' @param battenberg Output of [summarize_battenberg] (with QC columns).
#' @param groups Named vector of group labels, named by sample (as in
#' `battenberg$summary$sample`). Samples not in it are labelled NA.
#'
#' @return A list: `status` (per group and doubling status: samples, counts
#' per QC status and percent passing), `reasons` (per group: samples flagged
#' for each criterion) and `distributions` (per group and doubling status:
#' median and interquartile range of the clonal fit score, subclonal fraction,
#' purity, assessable fraction and |final - clonal fit ploidy|).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report = battenberg_qc_report(bb, groups = setNames(nmf_group, sample_ids))
#' report$status
#' }
battenberg_qc_report = function(battenberg, groups){
  s = battenberg$summary
  if(!"allelic_qc_status" %in% colnames(s)){
    stop("no QC columns; run battenberg_qc() first")
  }
  s$group = unname(groups[s$sample])
  s$wgd_heuristic = ifelse(s$doubled, "doubled", "not doubled")
  statuses = c("pass", "review", "severe_review", "unknown")
  cells = split(s, list(s$group, s$wgd_heuristic), drop = TRUE)
  status = do.call(rbind, lapply(cells, function(x){
    counts = table(factor(x$allelic_qc_status, levels = statuses))
    data.frame(group = x$group[1], doubled = x$wgd_heuristic[1], n = nrow(x),
               t(as.matrix(counts)),
               percent_pass = round(100 * counts[["pass"]] / nrow(x), 1),
               check.names = FALSE)
  }))
  flags = grep("^qc_", colnames(s), value = TRUE)
  reasons = do.call(rbind, lapply(split(s, s$group), function(x){
    data.frame(group = x$group[1], n = nrow(x), t(colSums(x[, flags, drop = FALSE])))
  }))
  s$abs_ploidy_delta = abs(s$ploidy_delta_clonal_final)
  metrics = c("clonal_fit_score", "subclonal_fraction", "purity",
              "assessable_fraction", "abs_ploidy_delta")
  distributions = do.call(rbind, lapply(cells, function(x){
    do.call(rbind, lapply(metrics, function(m){
      q = stats::quantile(x[[m]], c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
      data.frame(group = x$group[1], doubled = x$wgd_heuristic[1], metric = m,
                 n = sum(!is.na(x[[m]])), q25 = q[1], median = q[2], q75 = q[3])
    }))
  }))
  rownames(status) = rownames(reasons) = rownames(distributions) = NULL
  list(status = status, reasons = reasons, distributions = distributions)
}
