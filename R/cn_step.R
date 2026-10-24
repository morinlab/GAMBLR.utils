#' @title Copy number step at a position
#'
#' @description For each sample, the copy number on either side of a genomic
#' position and the difference between them (the "step"), e.g. to find
#' copy number transitions at a recurrent breakpoint, such as those left by
#' an unbalanced translocation, whether or not the rearrangement itself was
#' detected.
#'
#' @details Each side is the length-weighted mean copy number over `flank`
#' bases, leaving out `exclude` bases next to the position (where the exact
#' breakpoint varies between samples):
#' \itemize{
#'   \item left: from `pos - exclude - flank` to `pos - exclude`;
#'   \item right: from `pos + exclude` to `pos + exclude + flank`.
#' }
#' Sides stop at the chromosome ends, so a side can be shorter than `flank`
#' (e.g. the telomeric side of IGH, about 1 Mb from the 14q telomere).
#' `step` is right minus left, so on a q arm it is the telomeric side minus
#' the centromeric side. A side covered by segments for less than
#' `min_coverage` of its length is NA, except in samples from modules whose
#' gaps are copy neutral (controlfreec; see [add_segment_baseline]), whose
#' gaps count as 2. Placeholder rows (`dummy_segment == 1`) are ignored.
#'
#' @param seg_data Segments with ID, chrom, start, end and CN (centred or
#' not; see [center_segments]).
#' @param chrom,pos The position (chromosome with or without "chr").
#' @param flank Bases on each side to average over. Default 2e6.
#' @param exclude Bases next to the position left out on each side. Default
#' 1e5.
#' @param min_coverage Fraction of a side that must be covered by segments.
#' Default 0.5.
#' @param these_samples_metadata Optional metadata: only its samples, all of
#' them (those without segments get NA).
#' @param genome_build "grch37" (default) or "hg38", for the chromosome
#' lengths.
#'
#' @return A data frame with one row per sample: sample_id, cn_left,
#' cn_right and step.
#'
#' @export
#'
#' @examples
#' seg = data.frame(ID = c("s1", "s1", "s2"), chrom = "18",
#'                  start = c(5e7, 60.8e6, 5e7), end = c(60.8e6, 7e7, 7e7),
#'                  CN = c(3, 2, 2))
#' cn_step(seg, chrom = "18", pos = 60.8e6)
cn_step = function(seg_data, chrom, pos, flank = 2e6, exclude = 1e5, min_coverage = 0.5,
                   these_samples_metadata = NULL, genome_build = "grch37"){
  s = as.data.frame(seg_data)
  if(!all(c("ID", "chrom", "start", "end", "CN") %in% colnames(s))){
    stop("seg_data needs columns ID, chrom, start, end and CN")
  }
  ids = if(is.null(these_samples_metadata)) unique(as.character(s$ID))
        else unique(as.character(these_samples_metadata$sample_id))
  neutral_gaps = names(.fixed_baselines(s))
  s = s[.real_segments(s) & as.character(s$ID) %in% ids & !is.na(s$CN) &
          gsub("^chr", "", as.character(s$chrom)) == gsub("^chr", "", as.character(chrom)), ,
        drop = FALSE]
  side = function(from, to){
    ov = pmin(s$end, to) - pmax(s$start, from)
    keep = ov > 0
    id = factor(as.character(s$ID[keep]), levels = ids)
    covered = tapply(ov[keep], id, sum)
    total = tapply(ov[keep] * as.numeric(s$CN[keep]), id, sum)
    covered[is.na(covered)] = 0; total[is.na(total)] = 0
    width = to - from
    # copy-neutral gaps (controlfreec) count as 2
    gap = ids %in% neutral_gaps
    total[gap] = total[gap] + 2 * (width - covered[gap])
    covered[gap] = width
    value = as.numeric(total / covered)
    value[as.numeric(covered) < min_coverage * width] = NA
    value
  }
  # sides stop at the chromosome ends (e.g. IGH is about 1 Mb from 14qter)
  chrom_len = .chrom_lengths(genome_build)[gsub("^chr", "", as.character(chrom))]
  if(is.na(chrom_len)) chrom_len = Inf
  left = side(max(0, pos - exclude - flank), pos - exclude)
  right = side(pos + exclude, min(chrom_len, pos + exclude + flank))
  data.frame(sample_id = ids, cn_left = left, cn_right = right, step = right - left,
             stringsAsFactors = FALSE)
}

#' @title Concordance of copy number steps at two positions
#'
#' @description Measure, per sample, the copy number steps at two positions
#' (see [cn_step]) and whether they go in the directions expected if the
#' two positions are joined by a rearrangement, by group (e.g. samples with
#' and without a known translocation).
#'
#' @details An unbalanced reciprocal translocation leaves copy number steps
#' at both breakpoints: a gained derivative is higher on the centromeric
#' side of one and the telomeric side of the other. Steps are right minus
#' left (higher minus lower coordinate), which is telomeric minus
#' centromeric on a q arm but the reverse on a p arm, so the expected signs
#' are opposite when both loci are on q arms (or both on p arms), e.g.
#' t(14;18) (IGH and BCL2) or t(8;14) (MYC and IGH), and the same when one
#' is on a p arm, e.g. t(2;8) (IGK on 2p and MYC). This is worked out for
#' each partner from the chromosome arm coordinates of `genome_build`, unless
#' given (`expect`, or a partner's own `expect`).
#'
#' \strong{Several partners.} `b` can be a list of partner loci (e.g. IGH,
#' IGK and IGL for MYC). Each sample is then assigned its best-supported
#' partner (`b_partner`): a transition at both loci in the expected
#' directions first, then transitions at both, then a transition at the
#' partner, larger steps breaking ties; the `b_*` columns are that
#' partner's. `partner_steps` has every partner's step for every sample and
#' `partner_summary` the counts per partner and group.
#'
#' For each group: the number of samples with both steps measured, the
#' Spearman correlation of a's step and the partner's (oriented so the
#' expected combination is always negative; `b_step_oriented`),
#' and, among samples where both steps are at least `min_step` copies, how
#' many go the expected way (`n_expected`) or not (`n_other`), with a
#' binomial test against an even split. With exactly two groups, a Fisher
#' test compares their shares of expected steps.
#'
#' @param seg_data Segments with ID, chrom, start, end and CN.
#' @param a,b The two loci (`b` can also be a list of partner loci; see
#' Details), each a list with chrom, optionally a `name` (used in
#' `step_type`, default "a" and "b") and `expect`, and either `pos` (the
#' step at that position, as in [cn_step]) or `from` and `to` (each
#' sample's largest transition in that range, with its position, as in
#' [cn_transition]; optionally with its own `flank`, default 5e5), e.g.
#' `list(chrom = "18", from = 60.5e6, to = 61.5e6)`.
#' @param these_samples_metadata Metadata with sample_id (and `group_by`).
#' @param group_by Optional metadata column to compare groups by, e.g.
#' "bcl2_ba".
#' @param expect "opposite" or "same" signs for every partner. Default
#' NULL: worked out from the loci's chromosome arms (see Details).
#' @param min_step Smallest step (in copies) counted as a transition.
#' Default 0.5.
#' @param flank,exclude,min_coverage,genome_build As in [cn_step] (`flank`
#' and `exclude` for loci given as a position).
#'
#' @return A list: `steps` (one row per sample: sample_id, group, a's and
#' b's pos (the given position, or where in the range the transition is),
#' cn_left, cn_right and step, `b_partner` and `b_step_oriented`,
#' `pattern`: "expected", "other",
#' "one" (a step at only one position) or "none", and `step_type`: the
#' directions, e.g. "BCL2 down, IGH up", "BCL2 up only" or "none", a factor
#' with the expected combinations first, e.g. for `split_by` in
#' GAMBLR.viz::plot_multi_region_CN), `summary` (one row per
#' group), `partner_summary` and `partner_steps` (per partner), `fisher`
#' (the comparison of two groups, or NULL), and the `loci`, `expect` (per
#' partner) and `min_step` used.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # t(14;18): BCL2 (major breakpoint region) and IGH-J, grch37
#' conc = cn_step_concordance(seg,
#'                            a = list(chrom = "18", pos = 60793500),
#'                            b = list(chrom = "14", pos = 106329000),
#'                            these_samples_metadata = meta, group_by = "bcl2_ba")
#' conc$summary
#'
#' # MYC with any of the three immunoglobulin loci as its partner
#' conc = cn_step_concordance(seg,
#'          a = list(name = "MYC", chrom = "8", from = 128.2e6, to = 129.6e6),
#'          b = list(list(name = "IGH", chrom = "14", from = 105.5e6, to = 107.3e6),
#'                   list(name = "IGK", chrom = "2", from = 88.8e6, to = 90.3e6),
#'                   list(name = "IGL", chrom = "22", from = 22.3e6, to = 23.4e6)),
#'          these_samples_metadata = meta, group_by = "myc_ba")
#' conc$partner_summary
#' }
cn_step_concordance = function(seg_data, a, b, these_samples_metadata, group_by = NULL,
                               expect = NULL, min_step = 0.5, flank = 2e6,
                               exclude = 1e5, min_coverage = 0.5, genome_build = "grch37"){
  if(!is.null(expect)) expect = match.arg(expect, c("opposite", "same"))
  if(missing(these_samples_metadata)) stop("these_samples_metadata is required")
  meta = as.data.frame(these_samples_metadata)
  meta = meta[!duplicated(meta$sample_id), , drop = FALSE]
  # b: one locus or a list of partner loci
  partners = if(!is.null(b[["chrom"]])) list(b) else b
  for(k in seq_along(partners)){
    if(is.null(partners[[k]][["name"]])){
      partners[[k]]$name = if(length(partners) == 1) "b" else paste0("b", k)
    }
  }
  partner_names = vapply(partners, function(p) as.character(p[["name"]]), "")
  if(any(duplicated(partner_names))) stop("partner names in b must be unique")
  a_name = if(is.null(a[["name"]])) "a" else a[["name"]]
  step_at = function(p){
    if(!is.null(p[["from"]])){
      # a range: each sample's largest transition in it, with its position
      return(cn_transition(seg_data, chrom = p[["chrom"]], from = as.numeric(p[["from"]]),
                           to = as.numeric(p[["to"]]),
                           flank = if(is.null(p[["flank"]])) 5e5 else as.numeric(p[["flank"]]),
                           min_coverage = min_coverage, these_samples_metadata = meta,
                           genome_build = genome_build))
    }
    out = cn_step(seg_data, chrom = p[["chrom"]], pos = as.numeric(p[["pos"]]),
                  flank = flank, exclude = exclude, min_coverage = min_coverage,
                  these_samples_metadata = meta, genome_build = genome_build)
    out$pos = as.numeric(p[["pos"]])
    out
  }
  # the expected direction for each partner: opposite signs when both loci
  # are on the same kind of arm (both q, or both p), the same sign otherwise
  a_arm = .locus_arm(a, genome_build)
  partner_expect = vapply(partners, function(p){
    if(!is.null(p[["expect"]])) return(match.arg(p[["expect"]], c("opposite", "same")))
    if(!is.null(expect)) return(expect)
    if(.locus_arm(p, genome_build) == a_arm) "opposite" else "same"
  }, "")
  names(partner_expect) = partner_names
  sa = step_at(a)
  big_a = abs(sa$step) >= min_step
  # every partner's step, and whether it goes with a's as expected
  per_partner = lapply(seq_along(partners), function(k){
    sb = step_at(partners[[k]])
    big_b = abs(sb$step) >= min_step
    same_sign = sign(sa$step) == sign(sb$step)
    expected = if(partner_expect[k] == "opposite") !same_sign else same_sign
    data.frame(sample_id = meta$sample_id, partner = partner_names[k], expect = unname(partner_expect[k]),
               pos = sb$pos, cn_left = sb$cn_left, cn_right = sb$cn_right, step = sb$step,
               both = big_a & big_b, expected = expected, stringsAsFactors = FALSE)
  })
  # each sample's best-supported partner: steps at both loci in the expected
  # directions first, then at both, then a step at the partner; larger
  # steps break ties
  score = sapply(per_partner, function(x){
    sc = ifelse(is.na(x$step), -Inf,
         ifelse(x$both %in% TRUE & x$expected %in% TRUE, 3,
         ifelse(x$both %in% TRUE, 2, ifelse(abs(x$step) >= min_step, 1, 0))))
    # unmeasured partners rank last
    ifelse(is.na(x$step), -1e9,
           sc + pmin(abs(x$step), 20) / 100 + ifelse(is.na(sa$step), 0, pmin(abs(sa$step), 20) / 1000))
  })
  score = matrix(score, ncol = length(partners))
  best = max.col(score, ties.method = "first")
  chosen = do.call(rbind, lapply(seq_along(best), function(i) per_partner[[best[i]]][i, ]))
  steps = data.frame(sample_id = meta$sample_id,
                     group = if(is.null(group_by)) "all" else as.character(meta[[group_by]]),
                     a_pos = sa$pos, a_cn_left = sa$cn_left, a_cn_right = sa$cn_right, a_step = sa$step,
                     b_partner = chosen$partner,
                     b_pos = chosen$pos, b_cn_left = chosen$cn_left, b_cn_right = chosen$cn_right,
                     b_step = chosen$step, stringsAsFactors = FALSE)
  big_b = abs(steps$b_step) >= min_step
  expected = chosen$expected
  steps$pattern = ifelse(is.na(steps$a_step) | is.na(steps$b_step), NA_character_,
                  ifelse(big_a & big_b, ifelse(expected, "expected", "other"),
                  ifelse(big_a | big_b, "one", "none")))
  # the partner's step oriented so that the expected combination is always
  # opposite in sign to a's (for correlations across partners)
  steps$b_step_oriented = ifelse(chosen$expect == "same", -1, 1) * steps$b_step
  # the directions of the steps (tel - cen on q arms), named after the
  # loci; the expected combinations first
  dir = function(v) ifelse(v > 0, "up", "down")
  combos = lapply(seq_along(partners), function(k){
    bn = partner_names[k]
    cc = c(paste0(a_name, " down, ", bn, " up"), paste0(a_name, " up, ", bn, " down"),
           paste0(a_name, " up, ", bn, " up"), paste0(a_name, " down, ", bn, " down"))
    if(partner_expect[k] == "same") cc[c(3, 4, 1, 2)] else cc
  })
  types = c(unlist(lapply(combos, `[`, 1:2)), unlist(lapply(combos, `[`, 3:4)),
            paste0(a_name, c(" down", " up"), " only"),
            unlist(lapply(partner_names, function(bn) paste0(bn, c(" down", " up"), " only"))), "none")
  big_a = big_a %in% TRUE; big_b = big_b %in% TRUE
  step_type = ifelse(is.na(steps$pattern), NA_character_,
              ifelse(big_a & big_b, paste0(a_name, " ", dir(steps$a_step), ", ", steps$b_partner, " ",
                                           dir(steps$b_step)),
              ifelse(big_a, paste0(a_name, " ", dir(steps$a_step), " only"),
              ifelse(big_b, paste0(steps$b_partner, " ", dir(steps$b_step), " only"), "none"))))
  steps$step_type = factor(step_type, levels = types)
  groups = unique(steps$group)
  summary = do.call(rbind, lapply(groups, function(g){
    x = steps[steps$group %in% g & !is.na(steps$pattern), , drop = FALSE]
    rho = if(nrow(x) >= 5) suppressWarnings(stats::cor.test(x$a_step, x$b_step_oriented,
                                                            method = "spearman"))
          else NULL
    n_exp = sum(x$pattern == "expected"); n_oth = sum(x$pattern == "other")
    data.frame(group = g, n_measured = nrow(x),
               n_step_a = sum(abs(x$a_step) >= min_step),
               n_step_b = sum(abs(x$b_step) >= min_step),
               n_expected = n_exp, n_other = n_oth,
               binomial_p = if(n_exp + n_oth > 0) stats::binom.test(n_exp, n_exp + n_oth)$p.value else NA,
               spearman_rho = if(is.null(rho)) NA else unname(rho$estimate),
               spearman_p = if(is.null(rho)) NA else rho$p.value,
               stringsAsFactors = FALSE)
  }))
  # each partner on its own, per group
  partner_steps = do.call(rbind, per_partner)
  partner_steps$group = steps$group[match(partner_steps$sample_id, steps$sample_id)]
  partner_summary = do.call(rbind, lapply(partner_names, function(bn){
    do.call(rbind, lapply(groups, function(g){
      x = partner_steps[partner_steps$partner == bn & partner_steps$group %in% g &
                          !is.na(partner_steps$step) & !is.na(sa$step[match(partner_steps$sample_id, sa$sample_id)]), ,
                        drop = FALSE]
      both = x$both %in% TRUE
      data.frame(partner = bn, group = g, expect = partner_expect[match(bn, partner_names)],
                 n_measured = nrow(x), n_step = sum(abs(x$step) >= min_step),
                 n_expected = sum(both & x$expected), n_other = sum(both & !x$expected),
                 stringsAsFactors = FALSE)
    }))
  }))
  fisher = NULL
  if(length(groups) == 2){
    tab = rbind(summary$n_expected, summary$n_other)
    dimnames(tab) = list(c("expected", "other"), summary$group)
    if(all(colSums(tab) > 0)) fisher = list(table = tab, p.value = stats::fisher.test(tab)$p.value,
                                            odds_ratio = unname(stats::fisher.test(tab)$estimate))
  }
  list(steps = steps, summary = summary, partner_summary = partner_summary,
       partner_steps = partner_steps[, c("sample_id", "group", "partner", "expect", "pos", "cn_left",
                                         "cn_right", "step", "both", "expected")],
       fisher = fisher, loci = list(a = a, b = partners), expect = partner_expect,
       min_step = min_step)
}

# The arm ("p" or "q") a locus (a list with chrom and pos, or from and to)
# lies on, from the chromosome arm coordinates of the build
.locus_arm = function(locus, genome_build){
  arms = as.data.frame(if(genome_build %in% c("hg38", "grch38")) GAMBLR.data::chromosome_arms_hg38
                       else GAMBLR.data::chromosome_arms_grch37)
  arms$chromosome = gsub("^chr", "", as.character(arms$chromosome))
  chrom = gsub("^chr", "", as.character(locus[["chrom"]]))
  at = if(!is.null(locus[["pos"]])) as.numeric(locus[["pos"]])
       else mean(as.numeric(c(locus[["from"]], locus[["to"]])))
  q_start = arms$start[arms$chromosome == chrom & arms$arm == "q"]
  if(length(q_start) == 0) return("q")   # acrocentric or unknown: q
  if(at >= q_start[1]) "q" else "p"
}

#' @title The largest copy number transition in a range
#'
#' @description For each sample, find where in a range its copy number
#' changes most, and by how much: the segment boundary in the range with
#' the largest step, measured as in [cn_step] on either side of it. Useful
#' where a breakpoint isn't at one fixed position, e.g. to locate the copy
#' number transitions an unbalanced translocation leaves across a locus.
#'
#' @details Candidate positions are the sample's segment boundaries within
#' `from`-`to` (including the edges of gaps in controlfreec samples, whose
#' gaps are copy neutral). At each, the step is the length-weighted mean
#' copy number over `flank` bases to the right minus that to the left
#' (sides covered less than `min_coverage` don't count); the boundary with
#' the largest absolute step is kept. The position is the segment boundary,
#' so it locates the copy number change, not necessarily the breakpoint
#' itself.
#'
#' @param seg_data Segments with ID, chrom, start, end and CN.
#' @param chrom,from,to The range searched.
#' @param flank Bases on each side of a boundary to average over. Default
#' 5e5.
#' @param min_coverage,these_samples_metadata,genome_build As in [cn_step].
#'
#' @return A data frame with one row per sample: sample_id, pos, cn_left,
#' cn_right and step. A sample whose segments cover the range without a
#' boundary in it has step 0 and pos NA; one without enough coverage has
#' NA throughout.
#'
#' @export
#'
#' @examples
#' seg = data.frame(ID = "s1", chrom = "18", start = c(5e7, 60.9e6), end = c(60.9e6, 7e7),
#'                  CN = c(3, 2))
#' cn_transition(seg, chrom = "18", from = 60e6, to = 62e6)
cn_transition = function(seg_data, chrom, from, to, flank = 5e5, min_coverage = 0.5,
                         these_samples_metadata = NULL, genome_build = "grch37"){
  s = as.data.frame(seg_data)
  if(!all(c("ID", "chrom", "start", "end", "CN") %in% colnames(s))){
    stop("seg_data needs columns ID, chrom, start, end and CN")
  }
  ids = if(is.null(these_samples_metadata)) unique(as.character(s$ID))
        else unique(as.character(these_samples_metadata$sample_id))
  neutral_gaps = names(.fixed_baselines(s))
  chrom = gsub("^chr", "", as.character(chrom))
  chrom_len = .chrom_lengths(genome_build)[chrom]
  if(is.na(chrom_len)) chrom_len = Inf
  s = s[.real_segments(s) & as.character(s$ID) %in% ids & !is.na(s$CN) &
          gsub("^chr", "", as.character(s$chrom)) == chrom &
          s$end > from - flank & s$start < to + flank, , drop = FALSE]
  by_sample = split(s, as.character(s$ID))
  rows = lapply(ids, function(id){
    x = by_sample[[id]]
    none = data.frame(sample_id = id, pos = NA_real_, cn_left = NA_real_, cn_right = NA_real_,
                      step = NA_real_)
    if(is.null(x) || nrow(x) == 0) return(none)
    gap_cn = if(id %in% neutral_gaps) 2 else NA
    mean_cn = function(a, b){
      a = max(0, a); b = min(chrom_len, b)
      ov = pmin(x$end, b) - pmax(x$start, a)
      keep = ov > 0
      covered = sum(ov[keep])
      total = sum(ov[keep] * x$CN[keep])
      if(!is.na(gap_cn)){
        total = total + gap_cn * (b - a - covered); covered = b - a
      }
      if(covered < min_coverage * (b - a)) NA else total / covered
    }
    cand = sort(unique(c(x$start, x$end)))
    cand = cand[cand >= from & cand <= to]
    if(length(cand) == 0){
      # covered without a boundary: measured, with no transition
      level = mean_cn(from, to)
      if(is.na(level)) return(none)
      return(data.frame(sample_id = id, pos = NA_real_, cn_left = level, cn_right = level, step = 0))
    }
    left = vapply(cand, function(p) mean_cn(p - flank, p), numeric(1))
    right = vapply(cand, function(p) mean_cn(p, p + flank), numeric(1))
    step = right - left
    if(all(is.na(step))) return(none)
    k = which.max(abs(step))
    data.frame(sample_id = id, pos = cand[k], cn_left = left[k], cn_right = right[k],
               step = step[k])
  })
  do.call(rbind, rows)
}
