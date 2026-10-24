#' @title Check that segments cover each sample's genome
#'
#' @description Check, for each sample, that its copy number segments cover
#' enough of the autosomal genome for genome-wide estimates (such as a
#' baseline copy number) to be meaningful. Data that were filtered to a
#' region or a few chromosomes fail.
#'
#' @details Coverage is the union of each sample's autosomal segments
#' (overlaps are not double counted; placeholder rows with
#' `dummy_segment == 1` are left out) over the autosomal genome size of the
#' build (from `GAMBLR.data` chromosome arm coordinates). Complete GAMBL
#' samples typically cover 0.7-0.95 (callers leave centromeres, gaps and
#' placeholder regions uncalled); a region covers well under 0.01.
#'
#' @param seg_data Segments with ID, chrom, start and end.
#' @param genome_build "grch37" or "hg38". Default: the genome build of
#' `seg_data` if set, otherwise "grch37".
#' @param min_genome_coverage Minimum fraction of the autosomal genome a
#' sample's segments must cover. Default 0.5.
#' @param min_autosomes Minimum number of autosomes with segments. Default 20.
#'
#' @return A data frame with one row per sample: ID, coverage, n_autosomes
#' and ok (meeting both minimums).
#'
#' @export
#'
#' @examples
#' seg = data.frame(ID = "s1", chrom = as.character(1:22), start = 1, end = 5e7, CN = 2)
#' check_segment_coverage(seg)
check_segment_coverage = function(seg_data, genome_build = NULL, min_genome_coverage = 0.5,
                                  min_autosomes = 20){
  s = as.data.frame(seg_data)
  if(!all(c("ID", "chrom", "start", "end") %in% colnames(s))){
    stop("seg_data needs columns ID, chrom, start and end")
  }
  if(is.null(genome_build)) genome_build = get_genome_build(seg_data)
  if(is.null(genome_build)) genome_build = "grch37"
  autosome_bp = .autosome_length(genome_build)
  s$chrom = gsub("^chr", "", as.character(s$chrom))
  s = s[.is_autosome(s$chrom) & .real_segments(s), c("ID", "chrom", "start", "end"), drop = FALSE]
  ids = unique(as.character(as.data.frame(seg_data)$ID))
  if(nrow(s) == 0){
    return(data.frame(ID = ids, coverage = 0, n_autosomes = 0L, ok = FALSE))
  }
  s = s[order(s$ID, s$chrom, s$start), , drop = FALSE]
  n = nrow(s)
  g = as.integer(factor(paste(s$ID, s$chrom, sep = "\r")))
  g = cumsum(c(TRUE, g[-1] != g[-n]))
  # union length: each segment adds what lies beyond the furthest end so far
  offset = (g - 1) * (max(s$end) + 1)
  reach = c(-Inf, cummax(s$end + offset)[-n]) - offset
  reach[c(TRUE, g[-1] != g[-n])] = -Inf
  added = pmax(0, s$end - pmax(s$start, reach))
  coverage = tapply(added, s$ID, sum) / autosome_bp
  n_auto = tapply(s$chrom, s$ID, function(x) length(unique(x)))
  out = data.frame(ID = ids, coverage = as.numeric(coverage[ids]),
                   n_autosomes = as.integer(n_auto[ids]))
  out$coverage[is.na(out$coverage)] = 0
  out$n_autosomes[is.na(out$n_autosomes)] = 0L
  out$ok = out$coverage >= min_genome_coverage & out$n_autosomes >= min_autosomes
  out
}

#' @title Add each sample's baseline copy number to segments
#'
#' @description Add `baseline_cn` and `baseline_method` columns to copy
#' number segments: each sample's baseline (genome-wide) copy number, used
#' to call gains and losses relative to the sample rather than to diploid.
#'
#' @details Segments that already have a `baseline_cn` from the requested
#' method are returned unchanged (e.g. region segments from
#' `GAMBLR.results::get_cn_segments(region = ...)`, whose baselines were
#' computed genome-wide); a "mode" baseline is also converted to
#' "mode_rounded" without recomputing. Otherwise the baseline is computed from the
#' segments, which is only meaningful for whole-genome data: each sample is
#' first checked with [check_segment_coverage], and if any falls short (e.g.
#' the data were filtered to a region) this stops rather than return a
#' misleading baseline.
#'
#' Methods, over autosomal segments weighted by length (placeholder rows
#' with `dummy_segment == 1` are always left out). These are the same in
#' every GAMBLR function that takes a baseline or ploidy method (e.g.
#' [segmented_data_to_ploidy_estimate], [center_segments],
#' `GAMBLR.results::get_cn_segments`):
#' \itemize{
#'   \item "mode": the mean copy number within the sample's dominant copy
#'     number state (the rounded state covering the most of the genome).
#'     Robust to a large burden of gains or losses, and continuous, so a
#'     sample whose neutral level sits at 1.9 gets 1.9. This is the
#'     baseline `adjust_for_ploidy` in [segmented_data_to_cn_matrix] uses.
#'   \item "mode_rounded": that dominant state itself, a whole number
#'     (`round()` of "mode").
#'   \item "mean": the mean copy number over the genome, which reflects the
#'     burden of gains and losses (conventional ploidy).
#' }
#' Samples whose segments all come from controlfreec (a `module` column of
#' "controlfreec") get a baseline of 2 under every method, without the
#' coverage check: controlfreec doesn't report copy-neutral segments, so
#' most of the genome is missing and an estimate from the rest would be
#' wrong.
#'
#' Ties between dominant states go to the lower state. Only autosomes
#' (1-22) count, and copy number above 20 counts as 20, in every function
#' that estimates a baseline (including the cached baselines in
#' GAMBLR.results), so they all agree.
#' Copy number is the CN column, or `2 * 2^log.ratio` without one.
#'
#' @param seg_data Segments with ID, chrom, start, end and CN or log.ratio.
#' @param method "mode" (default), "mode_rounded" or "mean". See Details.
#' @param genome_build,min_genome_coverage,min_autosomes As in
#' [check_segment_coverage].
#'
#' @return `seg_data` with `baseline_cn` and `baseline_method` columns.
#'
#' @export
#'
#' @examples
#' seg = data.frame(ID = "s1", chrom = as.character(1:22), start = 1, end = 5e7,
#'                  CN = c(3, rep(2, 21)))
#' add_segment_baseline(seg)
add_segment_baseline = function(seg_data, method = "mode", genome_build = NULL,
                                min_genome_coverage = 0.5, min_autosomes = 20){
  method = match.arg(method, baseline_methods())
  existing = .existing_baseline(seg_data, method)
  if(!is.null(existing)){
    seg_data$baseline_cn = existing
    seg_data$baseline_method = method
    return(seg_data)
  }
  if("baseline_method" %in% colnames(seg_data)){
    message("recomputing baseline_cn with method \"", method, "\" (was \"",
            paste(unique(stats::na.omit(seg_data$baseline_method)), collapse = ", "), "\")")
  }
  fixed = .fixed_baselines(as.data.frame(seg_data))
  cov = check_segment_coverage(seg_data, genome_build, min_genome_coverage, min_autosomes)
  # samples with a fixed baseline don't need to cover the genome
  cov$ok[cov$ID %in% names(fixed)] = TRUE
  if(any(!cov$ok)){
    bad = cov[!cov$ok, , drop = FALSE]
    stop(nrow(bad), " sample(s) have segments covering too little of the genome to estimate ",
         "a baseline (less than ", min_genome_coverage, " of the autosomal genome or fewer than ",
         min_autosomes, " autosomes):\n  ", .list_ids(bad$ID), "\n",
         "If most samples fail, the data were probably filtered to a region: use ",
         "whole-genome segments, or segments with a baseline_cn column (e.g. from ",
         "GAMBLR.results::get_cn_segments(region = ...)). Otherwise, ",
         .remove_failing_hint(check = "coverage_ok"), call. = FALSE)
  }
  s = as.data.frame(seg_data)
  cn = if("CN" %in% colnames(s)) as.numeric(s$CN) else 2 * 2^as.numeric(s$log.ratio)
  auto = .is_autosome(s$chrom) & !is.na(cn) & .real_segments(s)
  base = .baseline_from(s$ID[auto], cn[auto], as.numeric(s$end[auto]) - as.numeric(s$start[auto]),
                        method)
  base[names(fixed)] = fixed
  seg_data$baseline_cn = unname(base[as.character(s$ID)])
  if(anyNA(seg_data$baseline_cn)) .stop_missing_baseline(seg_data$ID[is.na(seg_data$baseline_cn)])
  seg_data$baseline_method = method
  seg_data
}

# Autosomes (1-22, with or without a chr prefix): the chromosomes every
# baseline, ploidy and centring calculation uses (sex chromosomes and other
# contigs, e.g. unplaced or random ones, never count)
.is_autosome = function(chrom){
  gsub("^chr", "", as.character(chrom)) %in% as.character(1:22)
}

# Modules (callers) whose samples get a fixed baseline instead of an
# estimate: controlfreec doesn't report copy-neutral segments, so most of
# its genome is missing and an estimate from what's left is meaningless
.fixed_baseline_modules = c(controlfreec = 2)

# Fixed baselines, named by sample, for samples whose (real) segments all
# come from a module in .fixed_baseline_modules
.fixed_baselines = function(seg){
  if(!"module" %in% colnames(seg) || nrow(seg) == 0) return(stats::setNames(numeric(0), character(0)))
  real = .real_segments(seg)
  ids = as.character(seg$ID[real]); module = as.character(seg$module[real])
  only = tapply(module, ids, function(m) if(length(unique(m)) == 1) m[1] else NA_character_)
  only = only[!is.na(only) & only %in% names(.fixed_baseline_modules)]
  stats::setNames(as.numeric(.fixed_baseline_modules[only]), names(only))
}

# Highest copy number counted by baseline and ploidy estimates (higher
# values count as this), so the same sample gets the same baseline whatever
# cap a function applies to the values it returns
.baseline_max_cn = 20

# Rows that are real segments rather than placeholders (dummy_segment == 1),
# which never count towards a sample's baseline, ploidy or coverage
.real_segments = function(seg){
  if(!"dummy_segment" %in% colnames(seg)) return(rep(TRUE, nrow(seg)))
  dummy = suppressWarnings(as.numeric(seg$dummy_segment))
  is.na(dummy) | dummy != 1
}

#' @title Baseline (ploidy) methods
#'
#' @description The names of the methods for a sample's baseline copy
#' number used throughout GAMBLR: "mode", "mode_rounded" and "mean". See
#' [add_segment_baseline] for what each one is.
#'
#' @return A character vector.
#'
#' @export
#'
#' @examples
#' baseline_methods()
baseline_methods = function(){
  c("mode", "mode_rounded", "mean")
}

# The baseline_cn of the requested method from a baseline already on the
# segments: as is when the method matches, rounded when "mode_rounded" is
# wanted from "mode"; NULL when it would have to be recomputed
.existing_baseline = function(seg_data, method){
  if(!all(c("baseline_cn", "baseline_method") %in% colnames(seg_data))) return(NULL)
  have = unique(stats::na.omit(seg_data$baseline_method))
  base = as.numeric(seg_data$baseline_cn)
  if(length(have) != 1) return(NULL)
  if(!(have == method || (have == "mode" && method == "mode_rounded"))) return(NULL)
  # samples given no baseline (e.g. by GAMBLR.results::get_cn_segments, for
  # segments too sparse or without CN) can't be used
  if(anyNA(base)){
    .stop_missing_baseline(seg_data$ID[is.na(base)])
  }
  if(have == method) base else round(base)
}

# Stop on samples without a baseline, naming them and how to remove them
.stop_missing_baseline = function(ids){
  ids = unique(as.character(ids))
  stop(length(ids), " sample(s) have no baseline_cn (their segments are too sparse, run past ",
       "the chromosome ends, or have no CN):\n  ", .list_ids(ids), "\n",
       "Remove them first, e.g. with GAMBLR.results::annotate_cn_availability(metadata) and ",
       "dplyr::filter(metadata, cn_available) before getting the segments, or ",
       "seg_data = dplyr::filter(seg_data, !is.na(baseline_cn))", call. = FALSE)
}

# Per-sample baseline from segment copy numbers and lengths (methods as in
# add_segment_baseline), with copy number capped at max_cn
.baseline_from = function(id, cn, len, method, max_cn = .baseline_max_cn){
  cn = pmin(cn, max_cn)
  plain = function(x) stats::setNames(as.numeric(x), names(x))
  if(method == "mean"){
    return(plain(tapply(cn * len, id, sum) / tapply(len, id, sum)))
  }
  state = round(cn)
  tot = tapply(len, list(id, state), sum)
  tot[is.na(tot)] = 0
  dominant = stats::setNames(as.numeric(colnames(tot))[max.col(tot, ties.method = "first")],
                             rownames(tot))
  if(method == "mode_rounded") return(dominant)
  # mean copy number within the dominant state
  in_state = state == dominant[as.character(id)]
  plain(tapply(cn[in_state] * len[in_state], id[in_state], sum) /
          tapply(len[in_state], id[in_state], sum))
}

# Autosomal genome length (bp) for a build, from chromosome arm coordinates
.autosome_length = function(genome_build){
  arms = as.data.frame(if(genome_build %in% c("hg38", "grch38")) GAMBLR.data::chromosome_arms_hg38
                       else GAMBLR.data::chromosome_arms_grch37)
  arms$chromosome = gsub("^chr", "", as.character(arms$chromosome))
  arms = arms[arms$chromosome %in% as.character(1:22), , drop = FALSE]
  sum(as.numeric(tapply(arms$end, arms$chromosome, max)))
}

#' @title Check segmented data, sample by sample
#'
#' @description Run the checks that GAMBLR functions apply to copy number
#' segments and report, for each sample, whether it passes, so that failing
#' samples can be found and removed in one step: genome coverage (enough of
#' the genome for genome-wide estimates such as a baseline), genome build
#' (no segments beyond the chromosome ends of `genome_build`, which suggests
#' coordinates from another build) and overlapping segments. Samples whose
#' segments all come from controlfreec pass the coverage check, since they
#' get a fixed baseline (see [add_segment_baseline]).
#'
#' @param seg_data Segments with ID, chrom, start and end.
#' @param genome_build "grch37" or "hg38". Default: the genome build of
#' `seg_data` if set, otherwise "grch37".
#' @param min_genome_coverage,min_autosomes As in [check_segment_coverage].
#' @param chrom_end_tolerance How far (bp) a segment may run past a
#' chromosome end before it counts. Default 1e5.
#'
#' @return A data frame with one row per sample: ID, n_segments, coverage,
#' n_autosomes, n_beyond_chrom_end, n_overlapping, the checks coverage_ok,
#' build_ok and overlap_ok, and ok (all three).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qc = check_segments(seg, genome_build = "grch37")
#' qc[!qc$ok, ]
#' seg = dplyr::filter(seg, !ID %in% qc$ID[!qc$ok])
#' }
check_segments = function(seg_data, genome_build = NULL, min_genome_coverage = 0.5,
                          min_autosomes = 20, chrom_end_tolerance = 1e5){
  if(is.null(genome_build)) genome_build = get_genome_build(seg_data)
  if(is.null(genome_build)) genome_build = "grch37"
  s = as.data.frame(seg_data)
  out = check_segment_coverage(s, genome_build, min_genome_coverage, min_autosomes)
  out$n_segments = as.integer(table(factor(s$ID, levels = out$ID)))
  chrom_len = .chrom_lengths(genome_build)[gsub("^chr", "", as.character(s$chrom))]
  beyond = !is.na(chrom_len) & s$end > chrom_len + chrom_end_tolerance
  out$n_beyond_chrom_end = as.integer(table(factor(s$ID[beyond], levels = out$ID)))
  ov = check_segment_overlaps(s)
  out$n_overlapping = as.integer(table(factor(ov$ID, levels = out$ID)))
  # samples with a fixed baseline (controlfreec) don't need to cover the genome
  out$coverage_ok = out$ok | out$ID %in% names(.fixed_baselines(s))
  out$build_ok = out$n_beyond_chrom_end == 0
  out$overlap_ok = out$n_overlapping == 0
  out$ok = out$coverage_ok & out$build_ok & out$overlap_ok
  out[, c("ID", "n_segments", "coverage", "n_autosomes", "n_beyond_chrom_end",
          "n_overlapping", "coverage_ok", "build_ok", "overlap_ok", "ok")]
}

# Full chromosome lengths for a build, named without "chr" (the arm tables
# stop early on chrY, so they can't be used for this)
.chrom_lengths = function(genome_build){
  if(genome_build %in% c("hg38", "grch38")){
    return(c("1" = 248956422, "2" = 242193529, "3" = 198295559, "4" = 190214555,
             "5" = 181538259, "6" = 170805979, "7" = 159345973, "8" = 145138636,
             "9" = 138394717, "10" = 133797422, "11" = 135086622, "12" = 133275309,
             "13" = 114364328, "14" = 107043718, "15" = 101991189, "16" = 90338345,
             "17" = 83257441, "18" = 80373285, "19" = 58617616, "20" = 64444167,
             "21" = 46709983, "22" = 50818468, "X" = 156040895, "Y" = 57227415))
  }
  stats::setNames(as.numeric(.grch37_chrom_lengths), gsub("^chr", "", names(.grch37_chrom_lengths)))
}

# The lines to show users for removing samples that fail a check
.remove_failing_hint = function(object = "seg_data", check = "ok"){
  paste0("Find and remove them with:\n",
         "  qc = GAMBLR.utils::check_segments(", object, ")\n",
         "  ", object, " = dplyr::filter(", object, ", !ID %in% qc$ID[!qc$", check, "])")
}

# A list of IDs for an error message: all up to 30, then truncated
.list_ids = function(ids){
  ids = unique(ids)
  paste(if(length(ids) <= 30) ids else c(utils::head(ids, 30), "..."), collapse = ", ")
}
