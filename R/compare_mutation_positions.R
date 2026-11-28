#' @title Compare where in a protein mutations fall between two groups
#'
#' @description Tests whether the positions of a gene's mutations differ
#' between two groups of patients, three ways, each with a permutation test of
#' the patients' group labels (so ties at hotspots don't matter and a
#' patient's mutations always move together):
#' \itemize{
#'   \item overall: the Kolmogorov-Smirnov distance between the two groups'
#'     position distributions;
#'   \item per domain (when `domains` is given): the difference between the
#'     groups in the share of mutations inside the domain, with
#'     Benjamini-Hochberg q-values across domains;
#'   \item hotspots: whether the recurrent residues (positions mutated in at
#'     least `min_recurrence` patients, pooled) are shared differently, by a
#'     chi-square statistic on the residue-by-group table (other positions
#'     pooled into one row), with each residue also tested with Fisher's
#'     exact test (BH across residues).
#' }
#'
#' @details Each row of `data` counts once, so pass one row per patient per
#' residue (e.g. with `dplyr::distinct()`) for patient-level counting. Rows
#' with a missing group, position or patient are dropped.
#'
#' @param data A data frame with one row per mutation (or per patient and
#'   residue).
#' @param group_col Column with the two groups.
#' @param position_col Column with the amino acid position (numeric).
#' @param patient_col Column identifying the patient (the unit permuted).
#' @param domains Optional data frame of protein domains with columns `name`,
#'   `start` and `end`.
#' @param n_perm Number of permutations. Default 10000.
#' @param min_recurrence Residues mutated in at least this many rows (pooled)
#'   are hotspots. Default 3.
#' @param seed Random seed for the permutations, so results reproduce.
#'   Default 1.
#'
#' @return A list with `groups`, `overall` (n per group, KS distance D,
#'   permutation p), `domains` (per domain: mutations inside per group, share
#'   per group, difference, permutation p, q), `hotspots` (the chi-square test
#'   and per-residue table with Fisher p and q) and `n_perm`.
#'
#' @export
#'
#' @examples
#' set.seed(3)
#' d = data.frame(patient = paste0("p", 1:60), group = rep(c("A", "B"), each = 30),
#'                position = c(sample(c(20:30, 22, 22, 22), 30, TRUE), sample(100:300, 30)))
#' r = compare_mutation_positions(d, "group", "position", "patient", n_perm = 2000,
#'                                domains = data.frame(name = "N-term", start = 1, end = 50))
#' r$overall
compare_mutation_positions = function(data, group_col, position_col, patient_col, domains = NULL, n_perm = 10000,
                                      min_recurrence = 3, seed = 1){
  check_survival_columns(data, c(group_col, position_col, patient_col))
  d = as.data.frame(data)[, c(group_col, position_col, patient_col)]
  names(d) = c("group", "position", "patient")
  d = d[stats::complete.cases(d), , drop = FALSE]
  d$group = as.character(d$group)
  d$position = as.numeric(d$position)
  groups = unique(d$group)
  if(length(groups) != 2) stop("compare_mutation_positions needs exactly two groups; got ", length(groups))
  if(any(table(d$group) < 2)) stop("each group needs at least two mutations")
  # a patient belongs to one group; their rows move together when labels are permuted
  pg = unique(d[, c("patient", "group")])
  if(any(duplicated(pg$patient))) stop("some patients are in both groups")
  row_patient = match(d$patient, pg$patient)

  ks_d = function(lab){
    a = d$position[lab == groups[1]]; b = d$position[lab == groups[2]]
    grid = sort(unique(d$position))
    max(abs(stats::ecdf(a)(grid) - stats::ecdf(b)(grid)))
  }
  share_diff = function(lab, inside) mean(inside[lab == groups[1]]) - mean(inside[lab == groups[2]])

  hot = names(which(table(d$position) >= min_recurrence))
  cat_of = ifelse(as.character(d$position) %in% hot, as.character(d$position), "other")
  chisq_stat = function(lab){
    tab = table(factor(cat_of, levels = c(hot, "other")), factor(lab, levels = groups))
    tab = tab[rowSums(tab) > 0, , drop = FALSE]
    if(nrow(tab) < 2) return(0)
    e = outer(rowSums(tab), colSums(tab)) / sum(tab)
    sum((tab - e)^2 / e)
  }
  domain_inside = if(!is.null(domains) && nrow(domains) > 0){
    lapply(seq_len(nrow(domains)), function(i) d$position >= domains$start[i] & d$position <= domains$end[i])
  } else list()

  observed_ks = ks_d(d$group)
  observed_chi = if(length(hot) > 0) chisq_stat(d$group) else NA_real_
  observed_dom = vapply(domain_inside, function(ins) share_diff(d$group, ins), numeric(1))

  set.seed(seed)
  perm_ks = numeric(n_perm); perm_chi = numeric(n_perm)
  perm_dom = matrix(0, n_perm, length(domain_inside))
  for(k in seq_len(n_perm)){
    lab = sample(pg$group)[row_patient]
    perm_ks[k] = ks_d(lab)
    if(length(hot) > 0) perm_chi[k] = chisq_stat(lab)
    for(j in seq_along(domain_inside)) perm_dom[k, j] = share_diff(lab, domain_inside[[j]])
  }
  perm_p = function(obs, null) (1 + sum(null >= obs - 1e-12)) / (n_perm + 1)

  n = table(factor(d$group, levels = groups))
  overall = data.frame(group_1 = groups[1], group_2 = groups[2], n_1 = as.integer(n[1]), n_2 = as.integer(n[2]),
                       ks_d = observed_ks, p_value = perm_p(observed_ks, perm_ks), stringsAsFactors = FALSE)

  dom = NULL
  if(length(domain_inside) > 0){
    dom = do.call(rbind, lapply(seq_along(domain_inside), function(j){
      ins = domain_inside[[j]]
      data.frame(domain = domains$name[j], start = domains$start[j], end = domains$end[j],
                 inside_1 = sum(ins & d$group == groups[1]), inside_2 = sum(ins & d$group == groups[2]),
                 share_1 = mean(ins[d$group == groups[1]]), share_2 = mean(ins[d$group == groups[2]]),
                 difference = observed_dom[j], p_value = perm_p(abs(observed_dom[j]), abs(perm_dom[, j])),
                 stringsAsFactors = FALSE)
    }))
    dom$q_value = stats::p.adjust(dom$p_value, "BH")
  }

  hotspots = list(test = NULL, residues = NULL)
  if(length(hot) > 0){
    hotspots$test = data.frame(n_hotspots = length(hot), chisq = observed_chi, p_value = perm_p(observed_chi, perm_chi))
    res = do.call(rbind, lapply(hot, function(h){
      a = sum(d$position == as.numeric(h) & d$group == groups[1]); b = sum(d$position == as.numeric(h) & d$group == groups[2])
      ft = stats::fisher.test(matrix(c(a, n[1] - a, b, n[2] - b), 2))
      data.frame(position = as.numeric(h), n_1 = a, n_2 = b, share_1 = a / n[1], share_2 = b / n[2],
                 odds_ratio = unname(ft$estimate), p_value = ft$p.value)
    }))
    res$q_value = stats::p.adjust(res$p_value, "BH")
    rownames(res) = NULL
    hotspots$residues = res[order(res$p_value), , drop = FALSE]
  }
  list(groups = groups, overall = overall, domains = dom, hotspots = hotspots, n_perm = n_perm)
}
