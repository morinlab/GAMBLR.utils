#' @title Test co-occurrence and mutual exclusivity between pairs of events
#'
#' @description Test every pair of events (e.g. CNV features, mutated genes,
#' or both) in a sample-by-event status matrix for co-occurrence or mutual
#' exclusivity. Returns one row per pair with the 2x2 counts, an odds ratio, p
#' and adjusted q values, and flags for pairs whose events share a chromosome
#' arm or overlap, which are dependent by construction.
#'
#' @details
#' \strong{Status matrix:} rows are samples and columns are events. Values
#' greater than zero count as an event (so CNV levels or mutation counts can
#' be used as is); NA is treated as unknown, and each pair is tested on the
#' samples known for both events. Events present in fewer than `min_samples`
#' samples, or absent from fewer than `min_samples` samples, are not tested.
#'
#' \strong{Methods:}
#' \itemize{
#'   \item "fisher" (default): two-sided Fisher's exact test on the 2x2 table.
#'     The odds ratio is the conditional maximum likelihood estimate.
#'   \item "logistic": logistic regression of the second event on the first
#'     plus the per-sample `burden` covariate (e.g. percent genome altered or
#'     the number of events). The odds ratio is the first event's coefficient,
#'     adjusted for burden, and the p value is a likelihood ratio test of that
#'     coefficient (robust to pairs that never co-occur, where the
#'     coefficient itself is not well defined; the odds ratio is then
#'     reported as 0 or Inf). The model is not symmetric: event_b is the
#'     outcome and event_a the predictor. This corrects the main bias of
#'     unadjusted tests on CNV data: samples with unstable genomes carry many
#'     events, which makes most pairs look co-occurring and hides exclusivity.
#' }
#'
#' \strong{Leave-out burden:} a genome-wide burden such as PGA includes the
#' two events being tested, so a pair of large events (e.g. 1p loss and 1q
#' gain) raises the burden of the samples that carry both, and adjusting for
#' it partly adjusts away their association. This over-adjustment is largest
#' for arm-level events and negligible for focal ones. To avoid it, pass
#' per-chromosome PGA as `burden` (`calculate_pga(..., per_chromosome =
#' TRUE)`): each pair is then adjusted for the PGA of the rest of the genome,
#' leaving out the chromosomes of both events (weighted by chromosome length,
#' using the chromosome arms of `genome_build`). Events without coordinates
#' (e.g. mutations) leave nothing out.
#'
#' \strong{Gain and loss burden:} tumours differ in whether they mostly gain
#' or mostly lose DNA, and a single total burden can't account for that: the
#' leftover tendency makes gains look co-occurring with gains and exclusive
#' with losses (and the same for losses), which shows up as separate gain and
#' loss blocks when clustering. When `burden` has gain and loss columns (as
#' [calculate_pga] returns), both are used as covariates by default.
#'
#' \strong{Dependent pairs:} with `feature_info`, pairs are flagged when their
#' regions overlap (`overlapping`), such as an arm-level loss and a focal loss
#' inside it, or a gain and a loss of the same region; when they are on the
#' same chromosome arm (`same_arm`, which includes overlapping pairs); and
#' when they are on the same chromosome (`same_chromosome`). Overlapping and
#' same-arm pairs are partly or wholly dependent by construction, whereas
#' pairs on opposite arms of a chromosome (e.g. 1p loss and 1q gain, as from
#' an isochromosome 1q) can be biologically meaningful. `exclude` removes
#' flagged pairs before testing, so they also don't count towards the
#' multiple testing correction.
#' Events missing from `feature_info` (e.g. mutated genes when only CNV
#' features are given) get NA flags and are never excluded.
#'
#' @param status Sample-by-event matrix or data frame, e.g. from
#' [segmented_data_to_cnv_status] or [get_coding_ssm_status], or both joined
#' by sample. Sample IDs are taken from `sample_id_column` if it is a column,
#' and from rownames otherwise.
#' @param these_samples_metadata Required metadata whose sample_id column
#' defines the samples used, e.g. to run on a subset such as genomes only.
#' Samples in `status` but not in the metadata are dropped (and the burden is
#' subset to match); metadata samples missing from `status` are reported.
#' Also the source of `covariates`.
#' @param covariates Optional metadata columns to adjust for with
#' `method = "logistic"`, e.g. `c("seq_type")` or `c("seq_type", "purity")`.
#' Numeric columns are used as they are; other columns are treated as
#' categories (one indicator per level beyond the first). Samples with a
#' missing value are dropped. Use this for technical factors that affect
#' which events are detected, such as the assay, so that events that are
#' simply easier to detect in one assay don't look co-occurring.
#' @param feature_info Optional data frame with columns name, chrom, start and
#' end describing the events (e.g. [get_cnv_features], or the distinct name,
#' chrom, start and end columns of the long [segmented_data_to_cnv_status]
#' output). Used to flag and optionally exclude dependent pairs. Events are
#' matched to name, or to gene-labelled names such as "19p13.3~CD70_loss"
#' (from `gene_labels = TRUE`) when feature_info has representative_gene and
#' direction columns, as [get_cnv_features] does. Chromosome names are
#' compared without any "chr" prefix.
#' @param burden Per-sample burden covariate(s) for `method = "logistic"`: a
#' named numeric vector (names are sample IDs), a data frame with a
#' `sample_id_column` column and numeric column(s) chosen by `burden_column`
#' (e.g. the output of [calculate_pga]), or per-chromosome PGA (the output of
#' `calculate_pga(..., per_chromosome = TRUE)`) for a leave-out burden (see
#' Details; recommended). Ideally computed with the same centring as the
#' events, e.g. `centering = "median", cutoff = 0.1` for the
#' [segmented_data_to_cnv_status] defaults. Samples without a burden value
#' are dropped.
#' @param burden_column Column(s) of `burden` to use as covariates when it is
#' a data frame. Default NULL uses PGA_gain and PGA_loss together when
#' present (as from [calculate_pga]; see Details), otherwise the only numeric
#' column. For per-chromosome PGA, choose among "PGA", "PGA_gain" and
#' "PGA_loss" (e.g. "PGA" for the total only).
#' @param method "fisher" (default) or "logistic" (needs `burden` and/or
#' `covariates`). See Details.
#' @param events Optional character vector of event names (columns of
#' `status`) to test. Default all.
#' @param min_samples Minimum number of samples with, and without, an event
#' for it to be tested. Default 3.
#' @param exclude Pairs to remove before testing: "none", "overlapping",
#' "same_arm" (default; includes overlapping pairs) or "same_chromosome"
#' (includes both). Needs `feature_info`.
#' @param genome_build Genome build of `feature_info`, used for chromosome
#' arms (same_arm) and chromosome lengths (leave-out burden): "grch37"
#' (default) or "hg38".
#' @param p_adjust_method Method passed to [stats::p.adjust] for the q values.
#' Default "BH".
#' @param sample_id_column Name of the sample ID column of `status` (and of
#' `burden` when it is a data frame). Default "sample_id".
#'
#' @return A data frame with one row per tested pair, sorted by p value:
#' event_a, event_b, n_samples, n_both, n_a_only, n_b_only, n_neither,
#' odds_ratio, log_odds_ratio, p_value, q_value, association
#' ("co-occurring" or "mutually exclusive", from the odds ratio), and, with
#' `feature_info`, same_chromosome, same_arm, overlapping and each event's
#' coordinates
#' (chrom_a, start_a, end_a, chrom_b, start_b, end_b; NA for events not in
#' `feature_info`).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' segs = GAMBLR.open::get_cn_segments(these = dlbcl_meta)
#' cnv_status = segmented_data_to_cnv_status(segs, genome_build = "grch37",
#'                                           these_samples_metadata = dlbcl_meta)
#'
#' # unadjusted Fisher tests, without pairs of overlapping features
#' assoc = pairwise_event_association(cnv_status, dlbcl_meta,
#'                                    feature_info = get_cnv_features("grch37"))
#'
#' # adjusted for genomic instability: per-chromosome PGA, centred like the
#' # CNV calls, so each pair is adjusted for the rest of the genome
#' pga = calculate_pga(segs, centering = "median", cutoff = 0.1,
#'                     per_chromosome = TRUE)
#' assoc = pairwise_event_association(cnv_status, dlbcl_meta,
#'                                    feature_info = get_cnv_features("grch37"),
#'                                    burden = pga,
#'                                    method = "logistic")
#'
#' # several assays together, adjusting for the assay as well
#' assoc = pairwise_event_association(cnv_status, dlbcl_meta,
#'                                    feature_info = get_cnv_features("grch37"),
#'                                    burden = pga,
#'                                    covariates = "seq_type",
#'                                    method = "logistic")
#'
#' # the same on a subset of samples, e.g. genomes only
#' genome_meta = dplyr::filter(dlbcl_meta, seq_type == "genome")
#' assoc_genome = pairwise_event_association(cnv_status,
#'                                           these_samples_metadata = genome_meta,
#'                                           feature_info = get_cnv_features("grch37"),
#'                                           burden = pga,
#'                                           method = "logistic")
#'
#' # CNVs and mutations together
#' mutations = get_coding_ssm_status(gene_symbols = c("TP53", "EZH2", "MYD88"),
#'                                   these = dlbcl_meta, maf_data = maf)
#' both = dplyr::inner_join(tibble::rownames_to_column(cnv_status, "sample_id"),
#'                          mutations, by = "sample_id")
#' assoc = pairwise_event_association(both, dlbcl_meta,
#'                                    feature_info = get_cnv_features("grch37"))
#' }
pairwise_event_association = function(status,
                                      these_samples_metadata,
                                      feature_info = NULL,
                                      burden = NULL,
                                      covariates = NULL,
                                      burden_column = NULL,
                                      method = "fisher",
                                      events = NULL,
                                      min_samples = 3,
                                      exclude = "same_arm",
                                      p_adjust_method = "BH",
                                      sample_id_column = "sample_id",
                                      genome_build = "grch37"){
  if(!method %in% c("fisher", "logistic")){
    stop("method must be \"fisher\" or \"logistic\"")
  }
  if(!exclude %in% c("none", "overlapping", "same_arm", "same_chromosome")){
    stop("exclude must be \"none\", \"overlapping\", \"same_arm\" or \"same_chromosome\"")
  }
  arms = switch(genome_build,
                grch37 = , hg19 = GAMBLR.data::chromosome_arms_grch37,
                hg38 = , grch38 = GAMBLR.data::chromosome_arms_hg38,
                stop("genome_build must be \"grch37\" or \"hg38\""))
  arms = as.data.frame(arms)
  arms$chromosome = gsub("^chr", "", as.character(arms$chromosome))
  if(missing(these_samples_metadata) || is.null(these_samples_metadata)){
    stop("these_samples_metadata is required: its sample_id column defines the samples used")
  }
  if(method == "logistic" && is.null(burden) && is.null(covariates)){
    stop("method = \"logistic\" needs a burden and/or covariates")
  }
  if(method == "fisher" && !is.null(covariates)){
    stop("covariates need method = \"logistic\"")
  }
  these_samples_metadata = as.data.frame(these_samples_metadata)
  unknown_covariates = setdiff(covariates, colnames(these_samples_metadata))
  if(length(unknown_covariates) > 0){
    stop(paste("covariates not in these_samples_metadata:",
               paste(unknown_covariates, collapse = ", ")))
  }

  # status as a numeric matrix with sample IDs as rownames
  status = as.data.frame(status, check.names = FALSE)
  if(!is.null(sample_id_column) && sample_id_column %in% colnames(status)){
    rownames(status) = as.character(status[[sample_id_column]])
    status[[sample_id_column]] = NULL
  }
  if(!is.null(events)){
    unknown = events[!events %in% colnames(status)]
    if(length(unknown) > 0){
      stop(paste("events not in status:", paste(unknown, collapse = ", ")))
    }
    status = status[, events, drop = FALSE]
  }
  non_numeric = colnames(status)[!vapply(status, function(x)
    is.numeric(x) || is.logical(x), logical(1))]
  if(length(non_numeric) > 0){
    stop(paste("status columns must be numeric or logical:",
               paste(non_numeric, collapse = ", ")))
  }
  status = as.matrix(status) > 0
  storage.mode(status) = "integer"

  # samples: those in the metadata
  if(!"sample_id" %in% colnames(these_samples_metadata)){
    stop("these_samples_metadata needs a sample_id column")
  }
  wanted = unique(as.character(these_samples_metadata$sample_id))
  in_meta = rownames(status) %in% wanted
  if(sum(in_meta) == 0){
    stop("no samples in status are in these_samples_metadata")
  }
  missing_status = sum(!wanted %in% rownames(status))
  message("using ", sum(in_meta), " sample(s) in these_samples_metadata",
          if(any(!in_meta)) paste0(" (", sum(!in_meta),
                                   " other sample(s) in status dropped)")
          else "",
          if(missing_status > 0) paste0("; ", missing_status,
                                        " metadata sample(s) not in status")
          else "")
  status = status[in_meta, , drop = FALSE]
  meta_row = match(rownames(status), as.character(these_samples_metadata$sample_id))
  if(!is.null(covariates)){
    has_values = stats::complete.cases(
      these_samples_metadata[meta_row, covariates, drop = FALSE])
    if(sum(has_values) == 0){
      stop("no samples have values for all covariates")
    }
    if(any(!has_values)){
      message(sum(!has_values), " sample(s) with missing covariate values were dropped")
    }
    status = status[has_values, , drop = FALSE]
  }

  # burden covariates aligned to the samples. Each covariate is either
  # genome-wide (values) or per chromosome (altered length per chromosome,
  # so each pair can leave out its own chromosomes)
  burden_covs = NULL
  if(!is.null(burden)){
    burden_covs = .pairwise_burden_covariates(burden, burden_column,
                                              sample_id_column, arms,
                                              genome_build)
    keep = rownames(status) %in% burden_covs$samples
    if(sum(keep) == 0){
      stop("no samples in status have a burden value")
    }
    if(any(!keep)){
      message(sum(!keep), " sample(s) without a burden value were dropped")
    }
    status = status[keep, , drop = FALSE]
    burden_covs$covariates = lapply(burden_covs$covariates, function(cv){
      if(is.null(cv$altered)){
        cv$values = as.numeric(cv$values[rownames(status)])
      }else{
        cv$altered = cv$altered[rownames(status), , drop = FALSE]
        cv$total_altered = rowSums(cv$altered)
      }
      cv
    })
    message("adjusting for burden: ",
            paste(names(burden_covs$covariates), collapse = ", "),
            if(any(vapply(burden_covs$covariates,
                          function(cv) !is.null(cv$altered), logical(1))))
              " (leaving out each pair's chromosomes)" else "")
  }

  # metadata covariates as a design matrix aligned to the samples
  covariate_matrix = NULL
  if(!is.null(covariates)){
    meta_row = match(rownames(status),
                     as.character(these_samples_metadata$sample_id))
    parts = lapply(covariates, function(cv){
      values = these_samples_metadata[[cv]][meta_row]
      if(is.numeric(values)){
        return(stats::setNames(data.frame(values), cv))
      }
      values = factor(as.character(values))
      if(nlevels(values) < 2){
        message("covariate ", cv, " has a single value in these samples and was ignored")
        return(NULL)
      }
      indicators = stats::model.matrix(~ values)[, -1, drop = FALSE]
      colnames(indicators) = paste0(cv, "=", levels(values)[-1])
      as.data.frame(indicators, check.names = FALSE)
    })
    parts = parts[!vapply(parts, is.null, logical(1))]
    if(length(parts) > 0){
      covariate_matrix = as.matrix(do.call(cbind, parts))
      message("adjusting for covariates: ",
              paste(colnames(covariate_matrix), collapse = ", "))
    }
  }
  if(method == "logistic" && is.null(burden_covs) && is.null(covariate_matrix)){
    stop("nothing to adjust for: the covariates have a single value in these samples")
  }

  # testable events
  n_with = colSums(status == 1, na.rm = TRUE)
  n_without = colSums(status == 0, na.rm = TRUE)
  testable = n_with >= min_samples & n_without >= min_samples
  if(sum(testable) < 2){
    stop("fewer than two events occur in at least min_samples samples (and are absent from min_samples)")
  }
  status = status[, testable, drop = FALSE]
  event_names = colnames(status)

  pairs = utils::combn(length(event_names), 2)
  result = data.frame(event_a = event_names[pairs[1, ]],
                      event_b = event_names[pairs[2, ]],
                      stringsAsFactors = FALSE)

  # flag pairs that are dependent by construction
  if(!is.null(feature_info)){
    feature_info = as.data.frame(feature_info)
    required = c("name", "chrom", "start", "end")
    missing_cols = required[!required %in% colnames(feature_info)]
    if(length(missing_cols) > 0){
      stop(paste("feature_info is missing columns:",
                 paste(missing_cols, collapse = ", ")))
    }
    feature_info = feature_info[!duplicated(feature_info$name), , drop = FALSE]
    # events may use gene-labelled names (gene_labels = TRUE, e.g.
    # "19p13.3~CD70_loss"), so match those as well as plain names
    lookup_names = feature_info$name
    lookup_rows = seq_len(nrow(feature_info))
    if(all(c("representative_gene", "direction") %in% colnames(feature_info))){
      separators = unique(sub("^[A-Za-z0-9._-]*([^A-Za-z0-9._-]+).*$", "\\1",
                               grep("[^A-Za-z0-9._-]", event_names,
                                    value = TRUE)))
      for(sep in c("~", separators)){
        labelled = tryCatch(.label_cnv_feature_names(feature_info, TRUE, sep),
                            error = function(e) NULL)
        if(!is.null(labelled)){
          lookup_names = c(lookup_names, labelled)
          lookup_rows = c(lookup_rows, seq_len(nrow(feature_info)))
        }
      }
    }
    lookup_rows = lookup_rows[!duplicated(lookup_names)]
    lookup_names = lookup_names[!duplicated(lookup_names)]
    unmatched = event_names[!event_names %in% lookup_names]
    if(length(unmatched) == length(event_names)){
      warning("no events match feature_info names, so no pairs are flagged or excluded")
    }else if(length(unmatched) > 0){
      message(length(unmatched), " event(s) not in feature_info get no ",
              "chromosome flags: ",
              paste(utils::head(unmatched, 10), collapse = ", "),
              if(length(unmatched) > 10) ", ..." else "")
    }
    ia = lookup_rows[match(result$event_a, lookup_names)]
    ib = lookup_rows[match(result$event_b, lookup_names)]
    chrom = gsub("^chr", "", as.character(feature_info$chrom))
    result$same_chromosome = chrom[ia] == chrom[ib]
    # features that only touch end to end (e.g. 1p ending where 1q starts)
    # don't overlap
    result$overlapping = result$same_chromosome &
      feature_info$start[ia] < feature_info$end[ib] &
      feature_info$start[ib] < feature_info$end[ia]
    # chromosome arms each feature touches (a feature spanning the
    # centromere touches both)
    feature_chrom = gsub("^chr", "", as.character(feature_info$chrom))
    feature_arms = lapply(seq_len(nrow(feature_info)), function(k){
      on_chrom = arms[arms$chromosome == feature_chrom[k], , drop = FALSE]
      hit = on_chrom$start < feature_info$end[k] &
        feature_info$start[k] < on_chrom$end
      paste0(on_chrom$chromosome[hit], on_chrom$arm[hit])
    })
    shares_arm = mapply(function(x, y){
      if(is.na(x) || is.na(y)) return(NA)
      length(intersect(feature_arms[[x]], feature_arms[[y]])) > 0
    }, ia, ib)
    result$same_arm = result$overlapping | shares_arm
    # event coordinates, e.g. for ordering events along the genome in plots
    result$chrom_a = chrom[ia]
    result$start_a = feature_info$start[ia]
    result$end_a = feature_info$end[ia]
    result$chrom_b = chrom[ib]
    result$start_b = feature_info$start[ib]
    result$end_b = feature_info$end[ib]
    drop_pair = switch(exclude,
                       none = rep(FALSE, nrow(result)),
                       overlapping = result$overlapping,
                       same_arm = result$same_arm,
                       same_chromosome = result$same_chromosome)
    drop_pair[is.na(drop_pair)] = FALSE
    if(any(drop_pair)){
      message(sum(drop_pair), " pair(s) excluded (", exclude, ")")
    }
    result = result[!drop_pair, , drop = FALSE]
  }else if(exclude != "none"){
    message("exclude needs feature_info; no pairs were excluded")
  }
  if(nrow(result) == 0){
    stop("no pairs left to test")
  }

  pair_chrom_a = if("chrom_a" %in% colnames(result)) result$chrom_a else
    rep(NA_character_, nrow(result))
  pair_chrom_b = if("chrom_b" %in% colnames(result)) result$chrom_b else
    rep(NA_character_, nrow(result))
  if(!is.null(burden_covs) && is.null(feature_info) &&
     any(vapply(burden_covs$covariates,
                function(cv) !is.null(cv$altered), logical(1)))){
    message("per-chromosome burden without feature_info: no chromosomes can be left out, so the genome-wide burden is used")
  }

  test_pair = function(a, b, chrom_a, chrom_b){
    x = status[, a]
    y = status[, b]
    known = !is.na(x) & !is.na(y)
    x = x[known]
    y = y[known]
    counts = c(n_samples = length(x),
               n_both = sum(x == 1 & y == 1),
               n_a_only = sum(x == 1 & y == 0),
               n_b_only = sum(x == 0 & y == 1),
               n_neither = sum(x == 0 & y == 0))
    log_or = NA_real_
    p = NA_real_
    if(method == "fisher"){
      ft = tryCatch(stats::fisher.test(matrix(counts[c("n_both", "n_b_only",
                                                       "n_a_only", "n_neither")],
                                              nrow = 2)),
                    error = function(e) NULL)
      if(!is.null(ft)){
        log_or = log(unname(ft$estimate))
        p = ft$p.value
      }
    }else{
      design_null = matrix(1, nrow = sum(known), ncol = 1)
      if(!is.null(burden_covs)){
        leave_out = unique(stats::na.omit(c(chrom_a, chrom_b)))
        z = vapply(burden_covs$covariates, function(cv){
          if(is.null(cv$altered)){
            return(cv$values)
          }
          out = intersect(leave_out, colnames(cv$altered))
          (cv$total_altered -
             rowSums(cv$altered[, out, drop = FALSE])) /
            (sum(cv$chrom_len) - sum(cv$chrom_len[out]))
        }, numeric(nrow(status)))
        design_null = cbind(design_null, z[known, , drop = FALSE])
      }
      if(!is.null(covariate_matrix)){
        design_null = cbind(design_null,
                            covariate_matrix[known, , drop = FALSE])
      }
      design_full = cbind(design_null, x)
      fit_null = suppressWarnings(stats::glm.fit(design_null, y,
                                                 family = stats::binomial()))
      fit_full = suppressWarnings(stats::glm.fit(design_full, y,
                                                 family = stats::binomial()))
      log_or = unname(fit_full$coefficients[ncol(design_full)])
      # an empty cell separates the data, so the coefficient diverges; report
      # the limit, as Fisher's estimate does
      if(counts[["n_both"]] == 0 || counts[["n_neither"]] == 0){
        log_or = -Inf
      }else if(counts[["n_a_only"]] == 0 || counts[["n_b_only"]] == 0){
        log_or = Inf
      }
      p =stats::pchisq(max(fit_null$deviance - fit_full$deviance, 0),
                        df = 1, lower.tail = FALSE)
    }
    c(counts, log_odds_ratio = log_or, p_value = p)
  }

  tested = t(mapply(test_pair, result$event_a, result$event_b,
                    pair_chrom_a, pair_chrom_b))
  tested = as.data.frame(tested)
  rownames(tested) = NULL
  flags = result[, setdiff(colnames(result), c("event_a", "event_b")),
                 drop = FALSE]
  result = data.frame(result[, c("event_a", "event_b")], tested,
                      stringsAsFactors = FALSE)
  for(n in c("n_samples", "n_both", "n_a_only", "n_b_only", "n_neither")){
    result[[n]] = as.integer(result[[n]])
  }
  result$odds_ratio = exp(result$log_odds_ratio)
  result$q_value = stats::p.adjust(result$p_value, method = p_adjust_method)
  result$association = ifelse(is.na(result$log_odds_ratio), NA_character_,
                              ifelse(result$log_odds_ratio > 0,
                                     "co-occurring", "mutually exclusive"))
  result = cbind(result[, c("event_a", "event_b", "n_samples", "n_both",
                            "n_a_only", "n_b_only", "n_neither", "odds_ratio",
                            "log_odds_ratio", "p_value", "q_value",
                            "association")],
                 flags)
  result = result[order(result$p_value), , drop = FALSE]
  rownames(result) = NULL
  result
}

# Burden covariates for pairwise_event_association. Accepts a named vector,
# a data frame of per-sample values (e.g. calculate_pga: PGA, PGA_gain,
# PGA_loss) or per-chromosome PGA (chr<N>_pga[_gain|_loss] columns). Without
# burden_column, gain and loss burdens are used as two covariates when
# available (so tumours that tend to gain and tumours that tend to lose are
# both accounted for), otherwise the single total.
.pairwise_burden_covariates = function(burden, burden_column, sample_id_column,
                                       arms, genome_build){
  if(!is.data.frame(burden)){
    if(is.null(names(burden))){
      stop("burden must be named by sample ID, or be a data frame")
    }
    burden = burden[!is.na(burden)]
    return(list(samples = names(burden),
                covariates = list(burden = list(values = burden))))
  }
  burden = as.data.frame(burden)
  if(!sample_id_column %in% colnames(burden)){
    stop(paste0("burden needs a ", sample_id_column, " column"))
  }
  ids = as.character(burden[[sample_id_column]])
  per_chrom_cols = grep("^chr.+_pga(_gain|_loss)?$", colnames(burden),
                        value = TRUE)
  if(length(per_chrom_cols) > 1){
    # per chromosome: choose total or gain + loss
    suffixes = c(PGA = "_pga", PGA_gain = "_pga_gain", PGA_loss = "_pga_loss")
    available = names(suffixes)[vapply(suffixes, function(sfx)
      any(endsWith(per_chrom_cols, sfx) &
            !(sfx == "_pga" & grepl("_pga_(gain|loss)$", per_chrom_cols))),
      logical(1))]
    if(is.null(burden_column)){
      burden_column = if(all(c("PGA_gain", "PGA_loss") %in% available))
        c("PGA_gain", "PGA_loss") else "PGA"
    }
    unknown = burden_column[!burden_column %in% available]
    if(length(unknown) > 0){
      stop(paste0("for per-chromosome burden, burden_column must be among: ",
                  paste(available, collapse = ", ")))
    }
    covariates = lapply(stats::setNames(burden_column, burden_column),
                        function(bc){
      sfx = suffixes[[bc]]
      cols = grep(paste0("^chr.+", sfx, "$"), colnames(burden), value = TRUE)
      if(sfx == "_pga"){
        cols = cols[!grepl("_pga_(gain|loss)$", cols)]
      }
      chrom = gsub("^chr", "", sub(paste0(sfx, "$"), "", cols))
      chrom_len = vapply(chrom, function(ch){
        sum(arms$end[arms$chromosome == ch] - arms$start[arms$chromosome == ch])
      }, numeric(1))
      unknown_chrom = chrom[chrom_len == 0]
      if(length(unknown_chrom) > 0){
        stop(paste("burden has chromosomes not in the", genome_build,
                   "arms:", paste(unknown_chrom, collapse = ", ")))
      }
      altered = sweep(as.matrix(burden[, cols, drop = FALSE]), 2, chrom_len,
                      "*")
      colnames(altered) = chrom
      rownames(altered) = ids
      list(altered = altered, chrom_len = chrom_len)
    })
    complete = Reduce(intersect, lapply(covariates, function(cv)
      rownames(cv$altered)[stats::complete.cases(cv$altered)]))
    return(list(samples = complete, covariates = covariates))
  }
  # one value per sample
  numeric_cols = setdiff(colnames(burden)[vapply(burden, is.numeric,
                                                 logical(1))],
                         sample_id_column)
  if(is.null(burden_column)){
    if(all(c("PGA_gain", "PGA_loss") %in% numeric_cols)){
      burden_column = c("PGA_gain", "PGA_loss")
    }else if(length(numeric_cols) == 1){
      burden_column = numeric_cols
    }else{
      stop("burden has more than one numeric column; choose with burden_column")
    }
  }
  unknown = burden_column[!burden_column %in% numeric_cols]
  if(length(unknown) > 0){
    stop(paste("burden_column not a numeric column of burden:",
               paste(unknown, collapse = ", ")))
  }
  covariates = lapply(stats::setNames(burden_column, burden_column),
                      function(bc) list(values = stats::setNames(burden[[bc]],
                                                                 ids)))
  complete = ids[stats::complete.cases(burden[, burden_column, drop = FALSE])]
  list(samples = complete, covariates = covariates)
}
