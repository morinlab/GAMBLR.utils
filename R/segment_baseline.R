#' @title Check that segments cover each sample's genome
#'
#' @description Check, for each sample, that its copy number segments cover
#' enough of the autosomal genome for genome-wide estimates (such as a
#' baseline copy number) to be meaningful. Data that were filtered to a
#' region or a few chromosomes fail.
#'
#' @details Coverage is the union of each sample's autosomal segments
#' (overlaps are not double counted) over the autosomal genome size of the
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
  s = s[s$chrom %in% as.character(1:22), c("ID", "chrom", "start", "end"), drop = FALSE]
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
#' computed genome-wide). Otherwise the baseline is computed from the
#' segments, which is only meaningful for whole-genome data: each sample is
#' first checked with [check_segment_coverage], and if any falls short (e.g.
#' the data were filtered to a region) this stops rather than return a
#' misleading baseline.
#'
#' Methods, over autosomal segments weighted by length:
#' \itemize{
#'   \item "mode": the most common (rounded) copy number state.
#'   \item "mean": the mean copy number (can be fractional).
#' }
#' Copy number is the CN column, or `2 * 2^log.ratio` without one.
#'
#' @param seg_data Segments with ID, chrom, start, end and CN or log.ratio.
#' @param method "mode" (default) or "mean".
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
  method = match.arg(method, c("mode", "mean"))
  if(all(c("baseline_cn", "baseline_method") %in% colnames(seg_data)) &&
     all(seg_data$baseline_method == method, na.rm = TRUE) &&
     !anyNA(seg_data$baseline_cn)){
    return(seg_data)
  }
  if("baseline_method" %in% colnames(seg_data)){
    message("recomputing baseline_cn with method \"", method, "\" (was \"",
            paste(unique(stats::na.omit(seg_data$baseline_method)), collapse = ", "), "\")")
  }
  cov = check_segment_coverage(seg_data, genome_build, min_genome_coverage, min_autosomes)
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
  auto = gsub("^chr", "", as.character(s$chrom)) %in% as.character(1:22) & !is.na(cn)
  base = .baseline_from(s$ID[auto], cn[auto], as.numeric(s$end[auto]) - as.numeric(s$start[auto]),
                        method)
  seg_data$baseline_cn = unname(base[as.character(s$ID)])
  seg_data$baseline_method = method
  seg_data
}

# Per-sample baseline from segment copy numbers and lengths
.baseline_from = function(id, cn, len, method){
  if(method == "mean"){
    return(tapply(cn * len, id, sum) / tapply(len, id, sum))
  }
  tot = tapply(len, list(id, round(cn)), sum)
  tot[is.na(tot)] = 0
  stats::setNames(as.numeric(colnames(tot))[max.col(tot, ties.method = "first")], rownames(tot))
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
#' coordinates from another build) and overlapping segments.
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
  out$coverage_ok = out$ok
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
