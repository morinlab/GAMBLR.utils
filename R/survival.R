#' @title Summarize survival by group
#'
#' @description Kaplan-Meier summary of a time-to-event outcome, overall or
#' per group: number of patients and events, median follow-up (reverse
#' Kaplan-Meier), median survival with its 95% confidence interval, and
#' survival at fixed time points, with a log-rank test between groups.
#'
#' @details Rows with a missing time or event are dropped (counted in
#' `n_dropped`). Survival at a time point is NA for a group whose follow-up
#' ends before it. With `strata_col`, the log-rank test is stratified (e.g.
#' by study), so groups are compared within strata; the per-group summaries
#' still pool the strata.
#'
#' @param data A data frame with one row per patient.
#' @param time_col Column with the follow-up time (e.g. "OS_YEARS").
#' @param event_col Column with the event code: 1 = event, 0 = censored
#'   (e.g. "CODE_OS").
#' @param group_col Optional column to summarize and compare by.
#' @param strata_col Optional column to stratify the log-rank test by.
#' @param timepoints Times at which to report survival, in the units of
#'   `time_col`. Default: 2 and 5.
#'
#' @return A list with `summary` (a data frame with one row per group, or a
#'   single row "all": group, n, events, median_followup, median_survival,
#'   median_lower, median_upper, and surv_<t>, surv_<t>_lower,
#'   surv_<t>_upper for each time point), `logrank` (a list with chisq, df
#'   and p, or NULL without groups) and `n_dropped`.
#'
#' @export
#'
#' @examples
#' d = data.frame(OS_YEARS = c(1, 2, 3, 4, 5, 6, 7, 8),
#'                CODE_OS = c(1, 0, 1, 1, 0, 1, 0, 0),
#'                arm = rep(c("A", "B"), 4))
#' summarize_survival(d, "OS_YEARS", "CODE_OS", group_col = "arm")
summarize_survival = function(data, time_col = "OS_YEARS", event_col = "CODE_OS", group_col = NULL,
                              strata_col = NULL, timepoints = c(2, 5)){
  check_survival_columns(data, c(time_col, event_col, group_col, strata_col))
  d = as.data.frame(data)
  keep = !is.na(d[[time_col]]) & !is.na(d[[event_col]]) &
    (if(is.null(group_col)) TRUE else !is.na(d[[group_col]])) &
    (if(is.null(strata_col)) TRUE else !is.na(d[[strata_col]]))
  n_dropped = sum(!keep)
  d = d[keep, , drop = FALSE]
  if(nrow(d) == 0) stop("no rows with both a time and an event", call. = FALSE)
  d$.time = as.numeric(d[[time_col]])
  d$.event = as.numeric(d[[event_col]])
  d$.group = if(is.null(group_col)) "all" else as.character(d[[group_col]])

  one_group = function(g){
    x = d[d$.group == g, , drop = FALSE]
    fit = survival::survfit(survival::Surv(.time, .event) ~ 1, data = x)
    med = summary(fit)$table
    # median follow-up: reverse Kaplan-Meier (censoring as the event)
    fu = summary(survival::survfit(survival::Surv(.time, 1 - .event) ~ 1, data = x))$table
    row = data.frame(group = g, n = nrow(x), events = sum(x$.event == 1),
                     median_followup = unname(fu["median"]),
                     median_survival = unname(med["median"]),
                     median_lower = unname(med["0.95LCL"]), median_upper = unname(med["0.95UCL"]),
                     stringsAsFactors = FALSE)
    for(t in timepoints){
      s = summary(fit, times = t, extend = TRUE)
      beyond = t > max(x$.time)
      tag = paste0("surv_", format(t, trim = TRUE))
      row[[tag]] = if(beyond) NA_real_ else s$surv
      row[[paste0(tag, "_lower")]] = if(beyond) NA_real_ else s$lower
      row[[paste0(tag, "_upper")]] = if(beyond) NA_real_ else s$upper
    }
    row
  }
  groups = sort(unique(d$.group))
  summary_df = do.call(rbind, lapply(groups, one_group))
  rownames(summary_df) = NULL

  logrank = NULL
  if(!is.null(group_col) && length(groups) > 1){
    f = if(is.null(strata_col)) survival::Surv(.time, .event) ~ .group else
      stats::as.formula(paste0("survival::Surv(.time, .event) ~ .group + strata(", strata_col, ")"))
    sd = survival::survdiff(f, data = d)
    df = length(groups) - 1
    logrank = list(chisq = unname(sd$chisq), df = df,
                   p = stats::pchisq(sd$chisq, df = df, lower.tail = FALSE),
                   stratified_by = strata_col)
  }
  list(summary = summary_df, logrank = logrank, n_dropped = n_dropped)
}

#' @title Cox proportional hazards model
#'
#' @description Fit a Cox model and return hazard ratios with 95% confidence
#' intervals and Wald p-values, the numbers of patients and events used, and
#' a test of the proportional hazards assumption for each term.
#'
#' @details Rows with a missing time, event, covariate or stratum are dropped
#' (counted in `n_dropped`). Character covariates become factors whose
#' reference level is the one named in `reference`, else the most common
#' level. With `strata_col`, each stratum (e.g. study) has its own baseline
#' hazard and no hazard ratio of its own.
#'
#' @param data A data frame with one row per patient.
#' @param time_col Column with the follow-up time.
#' @param event_col Column with the event code (1 = event, 0 = censored).
#' @param covariates Columns to include as covariates.
#' @param strata_col Optional column to stratify the baseline hazard by.
#' @param reference Optional named list giving the reference level of
#'   categorical covariates, e.g. `list(COO_consensus = "GCB")`.
#'
#' @return A list with `coefficients` (a data frame: term, hr, hr_lower,
#'   hr_upper, p, ph_p), `n`, `events`, `ph_global_p`, `stratified_by` and
#'   `n_dropped`.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' d = data.frame(OS_YEARS = rexp(60), CODE_OS = rbinom(60, 1, 0.7),
#'                mutated = rep(c("yes", "no"), 30), study = rep(c("a", "b", "c"), 20))
#' survival_cox(d, "OS_YEARS", "CODE_OS", covariates = "mutated", strata_col = "study",
#'              reference = list(mutated = "no"))
survival_cox = function(data, time_col = "OS_YEARS", event_col = "CODE_OS", covariates, strata_col = NULL,
                        reference = list()){
  check_survival_columns(data, c(time_col, event_col, covariates, strata_col))
  d = as.data.frame(data)
  keep = stats::complete.cases(d[, c(time_col, event_col, covariates, strata_col), drop = FALSE])
  n_dropped = sum(!keep)
  d = d[keep, , drop = FALSE]
  if(nrow(d) == 0) stop("no complete rows for the model", call. = FALSE)
  for(cv in covariates){
    if(is.character(d[[cv]]) || is.logical(d[[cv]])){
      v = as.character(d[[cv]])
      ref = if(!is.null(reference[[cv]])) reference[[cv]] else names(sort(table(v), decreasing = TRUE))[1]
      if(!ref %in% v) stop("reference level ", ref, " isn't a value of ", cv, call. = FALSE)
      d[[cv]] = stats::relevel(factor(v), ref = ref)
    }
  }
  rhs = paste(sprintf("`%s`", covariates), collapse = " + ")
  if(!is.null(strata_col)) rhs = paste0(rhs, " + strata(`", strata_col, "`)")
  f = stats::as.formula(paste0("survival::Surv(`", time_col, "`, `", event_col, "`) ~ ", rhs))
  fit = survival::coxph(f, data = d)
  s = summary(fit)
  ci = s$conf.int
  coefs = data.frame(term = rownames(ci), hr = ci[, "exp(coef)"], hr_lower = ci[, "lower .95"],
                     hr_upper = ci[, "upper .95"], p = s$coefficients[, "Pr(>|z|)"],
                     stringsAsFactors = FALSE)
  zph = tryCatch(survival::cox.zph(fit)$table, error = function(e) NULL)
  coefs$ph_p = NA_real_
  ph_global = NA_real_
  if(!is.null(zph)){
    # cox.zph reports one test per covariate (not per factor level)
    for(cv in covariates){
      rows = startsWith(coefs$term, sprintf("`%s`", cv)) | startsWith(coefs$term, cv)
      if(cv %in% rownames(zph)) coefs$ph_p[rows] = zph[cv, "p"]
    }
    if("GLOBAL" %in% rownames(zph)) ph_global = zph["GLOBAL", "p"]
  }
  coefs$term = gsub("`", "", coefs$term)
  rownames(coefs) = NULL
  list(coefficients = coefs, n = fit$n, events = fit$nevent, ph_global_p = ph_global,
       stratified_by = strata_col, n_dropped = n_dropped)
}

#' @title Table 1: patient characteristics
#'
#' @description Summarize patient characteristics the way a clinical paper's
#' Table 1 does, overall or by group: counts and percentages for categorical
#' variables, median and interquartile range for numeric ones, the number
#' missing, and (with groups) a p-value per variable.
#'
#' @details Percentages are of the patients with a value (missing patients
#' are counted separately). P-values compare groups with Fisher's exact test
#' (simulated for tables larger than 2 x 2) for categorical variables and the
#' Kruskal-Wallis test (the Wilcoxon rank-sum test for two groups) for
#' numeric ones. Rows with no group value are dropped.
#'
#' @param data A data frame with one row per patient.
#' @param variables Columns to summarize.
#' @param group_col Optional column to split the table by.
#' @param test Whether to add p-values when grouped. Default: TRUE.
#'
#' @return A data frame with one row per variable level (categorical) or per
#'   numeric variable: variable, level, overall, one column per group, and p.
#'   Cells are "n (pct%)" or "median [Q1, Q3]"; the level "(missing)" counts
#'   patients without a value.
#'
#' @export
#'
#' @examples
#' d = data.frame(age = c(45, 60, 72, 55, 66, 38), sex = c("F", "M", "M", "F", "M", NA),
#'                group = c("A", "A", "B", "B", "B", "A"))
#' table_one(d, c("age", "sex"), group_col = "group")
table_one = function(data, variables, group_col = NULL, test = TRUE){
  check_survival_columns(data, c(variables, group_col))
  d = as.data.frame(data)
  if(!is.null(group_col)) d = d[!is.na(d[[group_col]]), , drop = FALSE]
  groups = if(is.null(group_col)) character(0) else sort(unique(as.character(d[[group_col]])))
  pct = function(k, n) if(n == 0) "0 (0%)" else sprintf("%d (%.0f%%)", k, 100 * k / n)
  med = function(x){
    x = x[!is.na(x)]
    if(length(x) == 0) return("")
    q = stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
    sprintf("%s [%s, %s]", signif(q[2], 3), signif(q[1], 3), signif(q[3], 3))
  }
  header = data.frame(variable = "N", level = "", overall = as.character(nrow(d)), stringsAsFactors = FALSE)
  for(g in groups) header[[g]] = as.character(sum(d[[group_col]] == g))
  if(length(groups) > 0) header$p = NA_real_
  out = list(header)
  for(v in variables){
    x = d[[v]]
    sub = function(g) if(is.null(g)) x else x[d[[group_col]] == g]
    if(is.numeric(x)){
      rows = data.frame(variable = v, level = "median [IQR]", overall = med(x), stringsAsFactors = FALSE)
      for(g in groups) rows[[g]] = med(sub(g))
      p = NA_real_
      if(test && length(groups) > 1 && sum(!is.na(x)) > 0){
        p = tryCatch(if(length(groups) == 2) stats::wilcox.test(x ~ d[[group_col]], exact = FALSE)$p.value else
          stats::kruskal.test(x ~ factor(d[[group_col]]))$p.value, error = function(e) NA_real_)
      }
    }else{
      x = as.character(x)
      levels_v = sort(unique(stats::na.omit(x)))
      rows = data.frame(variable = v, level = levels_v,
                        overall = vapply(levels_v, function(l) pct(sum(x %in% l), sum(!is.na(x))), ""),
                        stringsAsFactors = FALSE)
      for(g in groups){
        xg = sub(g)
        rows[[g]] = vapply(levels_v, function(l) pct(sum(xg %in% l), sum(!is.na(xg))), "")
      }
      p = NA_real_
      if(test && length(groups) > 1 && length(levels_v) > 1){
        tab = table(x, d[[group_col]])
        p = tryCatch(if(all(dim(tab) == 2)) stats::fisher.test(tab)$p.value else
          stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value, error = function(e) NA_real_)
      }
    }
    n_missing = sum(is.na(x))
    if(n_missing > 0){
      miss = data.frame(variable = v, level = "(missing)", overall = as.character(n_missing),
                        stringsAsFactors = FALSE)
      for(g in groups) miss[[g]] = as.character(sum(is.na(sub(g))))
      rows = rbind(rows, miss)
    }
    if(length(groups) > 0){
      rows$p = NA_real_
      rows$p[1] = p
    }
    out[[length(out) + 1]] = rows
  }
  res = do.call(rbind, out)
  rownames(res) = NULL
  res
}

# Stop with the names of any columns missing from data
check_survival_columns = function(data, columns){
  missing = setdiff(columns, colnames(data))
  if(length(missing) > 0){
    stop("column(s) not in data: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

# columns survfit/survdiff formulas refer to inside summarize_survival()
utils::globalVariables(c(".time", ".event", ".group"))
