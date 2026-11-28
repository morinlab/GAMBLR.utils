#' @title Score and choose among several segmentations of the same samples
#'
#' @description For each sample, score every segmentation (e.g. every
#' pipeline's segments from `GAMBLR.results::get_cn_segments(..., pipeline =
#' "all")`) with a QC penalty and mark the lowest-scoring one as `selected`.
#' The penalties are reported separately, so the choice can be checked, e.g.
#' with `GAMBLR.viz::plot_sample_cn_comparison`.
#'
#' Run it on the segments of the whole cohort or project at once, with
#' `these_samples_metadata` for all of those samples, not one sample at a
#' time: the over-segmentation penalty compares each segmentation with the
#' median over all samples of the same seq_type and genome version, and with
#' only a few samples that median is unreliable.
#'
#' @details Each segmentation is scored on its own QC metrics
#' ([segmented_data_qc], autosomes only), not by comparison with the others.
#' The score is the sum of these penalties (0 is best):
#' \itemize{
#'   \item `pen_coverage`: `weights$coverage` x how far the fraction of the
#'     autosomes covered by segments falls below `min_coverage` (0.9).
#'   \item `pen_oversegmentation`: `weights$oversegmentation` x
#'     log2(segments per 100 Mb / the reference median), when positive. The
#'     reference is the median segments per 100 Mb over every scored
#'     segmentation (all samples and sources) of the same seq_type and
#'     native genome version (grch37, hg19 and hs37d5 count as one version;
#'     hg38), so a segmentation's score doesn't depend on which other
#'     sources the sample has. Pass `reference` (e.g. the
#'     `oversegmentation_reference` attribute of an earlier result) to reuse
#'     fixed medians.
#'   \item `pen_short`: `weights$short` x the fraction of segments shorter
#'     than `short_segment` (1 Mb).
#'   \item `pen_deep_loss`: `weights$deep_loss` x the fraction of the covered
#'     genome below 0.5 copies (homozygous deletions are small, so a large
#'     fraction means a failed profile).
#'   \item `pen_contradiction`: `weights$contradiction` when deep losses cover
#'     more than 1% of the genome but their log ratio is above -0.5 (the copy
#'     number contradicts the data, e.g. a failed absolute fit).
#' }
#' With `eligible`, segmentations marked FALSE (e.g. ones that can't give a
#' genome-wide baseline, which functions such as
#' [segmented_data_to_cn_matrix] need) are ranked after the eligible ones and
#' never selected, so a sample whose segmentations are all ineligible has
#' none selected.
#'
#' Copy number is the CN column; a segmentation without CN (e.g. sequenza
#' capture segments) is scored on `2 * 2^log.ratio`. Ties are broken by
#' `prefer`, then alphabetically.
#'
#' A QC score doesn't check that a profile is right: a smooth but wrong
#' profile can score well. Look at the selected segmentation next to the
#' others before relying on it.
#'
#' @param seg_data Segments of the whole cohort or project from several
#' sources (see Description), with ID, chrom, start, end, CN and/or
#' log.ratio, and the column named in `by`.
#' @param these_samples_metadata Metadata with sample_id, seq_type and
#' genome_build (the native genome version), giving each sample's stratum for
#' the over-segmentation reference. Segments of samples not in it are
#' dropped. Default NULL: all samples form one stratum (with a message).
#' @param by Column naming each segment's source. Default "pipeline".
#' @param exclude Sources never scored or selected. Default the published
#' segmentations bundled in GAMBLR.data ("dlbclass", "schmitz_array"), which
#' are references, not GAMBL results.
#' @param prefer Optional order of sources used to break ties.
#' @param weights Named list of penalty weights: coverage (10),
#' oversegmentation (1), short (1), deep_loss (20) and contradiction (1).
#' Missing names keep their defaults.
#' @param min_coverage Fraction of the autosomes a segmentation should cover.
#' Default 0.9.
#' @param short_segment Length (bp) below which a segment counts as short.
#' Default 1e6.
#' @param genome_build Genome build of the coordinates, for the autosome
#' length ("grch37" or "hg38"). Default: from `seg_data`, else "grch37".
#' @param eligible Optional data frame with ID, the column named in `by`, and
#' eligible (logical): which segmentations may be selected (see Details).
#' Segmentations not listed are eligible.
#' @param reference Optional fixed over-segmentation medians: a data frame
#' with seq_type, genome and segments_per_100mb, as in the
#' `oversegmentation_reference` attribute of the result. Strata it lacks are
#' computed from `seg_data` (with a message).
#'
#' @return A data frame with one row per sample and source: ID and sample_id
#' (the same value; sample_id lets the table be joined to metadata with e.g.
#' `dplyr::left_join(metadata, scores)` without naming columns), seq_type and
#' genome (the stratum; NA without metadata), the source
#' (named as `by`), the QC metrics used (n_segments, covered_fraction,
#' segments_per_100mb, short_segment_fraction, deep_loss_fraction,
#' deep_loss_logratio), each penalty, score, eligible, rank (1 = best) and
#' selected (TRUE for each sample's chosen source). Excluded sources aren't included.
#' The medians used are in `attr(result, "oversegmentation_reference")`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # the whole cohort (or project), not one sample
#' meta = GAMBLR.results::get_gambl_metadata() %>% dplyr::filter(cohort == "dlbcl_chapuy")
#' segs = GAMBLR.results::get_cn_segments(meta, pipeline = "all")
#' scores = compare_segmentations(segs, these_samples_metadata = meta)
#' dplyr::filter(scores, selected)
#' GAMBLR.viz::plot_sample_cn_comparison(segs, "DLBCL10541T", selection = scores)
#' }
compare_segmentations = function(seg_data,
                                  these_samples_metadata = NULL,
                                  by = "pipeline",
                                  exclude = c("dlbclass", "schmitz_array"),
                                  prefer = NULL,
                                  weights = list(),
                                  min_coverage = 0.9,
                                  short_segment = 1e6,
                                  genome_build = NULL,
                                  eligible = NULL,
                                  reference = NULL){
  w = utils::modifyList(list(coverage = 10, oversegmentation = 1, short = 1, deep_loss = 20,
                             contradiction = 1), weights)
  if(is.null(genome_build)){
    genome_build = tryCatch(get_genome_build(seg_data), error = function(e) NULL)
    if(is.null(genome_build) || !genome_build %in% c("grch37", "hg38")) genome_build = "grch37"
  }
  autosome_mb = sum(as.numeric(.chrom_lengths(genome_build)[as.character(1:22)])) / 1e6
  segs = as.data.frame(seg_data)
  if(!by %in% colnames(segs)) stop("seg_data has no '", by, "' column")
  segs$.source = as.character(segs[[by]])
  segs = segs[!segs$.source %in% exclude, , drop = FALSE]
  if(nrow(segs) == 0) stop("no segments left after excluding: ", paste(exclude, collapse = ", "))
  if(!is.null(these_samples_metadata)){
    meta = as.data.frame(these_samples_metadata)
    missing_cols = setdiff(c("sample_id", "seq_type", "genome_build"), colnames(meta))
    if(length(missing_cols)) stop("these_samples_metadata lacks: ", paste(missing_cols, collapse = ", "))
    outside = setdiff(unique(segs$ID), meta$sample_id)
    if(length(outside)){
      message("dropping ", length(outside), " sample(s) not in these_samples_metadata")
      segs = segs[segs$ID %in% meta$sample_id, , drop = FALSE]
      if(nrow(segs) == 0) stop("no segments left for the samples in these_samples_metadata")
    }
    meta = meta[!duplicated(meta$sample_id), , drop = FALSE]
  }else{
    message("no these_samples_metadata: over-segmentation is measured against one stratum of all samples given")
    meta = NULL
  }

  qc = do.call(rbind, lapply(split(segs, segs$.source), function(s){
    # a source without CN is scored on its log ratios
    scale = if("CN" %in% colnames(s) && any(!is.na(s$CN))) "CN" else "log_ratio"
    q = segmented_data_qc(s, value_scale = scale, short_segment = short_segment)
    q$source = s$.source[1]
    q
  }))
  qc$covered_fraction = pmin(qc$covered_mb / autosome_mb, 1)
  qc$pen_coverage = w$coverage * pmax(min_coverage - qc$covered_fraction, 0)
  # over-segmentation reference: median per seq_type and native genome version
  if(is.null(meta)){
    qc$seq_type = NA_character_
    qc$genome = NA_character_
  }else{
    m = match(qc$sample_id, meta$sample_id)
    qc$seq_type = as.character(meta$seq_type[m])
    qc$genome = .native_genome_version(meta$genome_build[m])
  }
  stratum = paste(qc$seq_type, qc$genome)
  first = !duplicated(stratum)
  computed = data.frame(seq_type = qc$seq_type[first], genome = qc$genome[first],
                        segments_per_100mb = as.vector(tapply(qc$segments_per_100mb, stratum, stats::median)[stratum[first]]),
                        n_segmentations = as.vector(table(stratum)[stratum[first]]),
                        stringsAsFactors = FALSE)
  if(is.null(reference)){
    ref = computed
  }else{
    reference = as.data.frame(reference)
    missing_cols = setdiff(c("seq_type", "genome", "segments_per_100mb"), colnames(reference))
    if(length(missing_cols)) stop("reference lacks: ", paste(missing_cols, collapse = ", "))
    reference$genome = .native_genome_version(reference$genome)
    new = computed[!paste(computed$seq_type, computed$genome) %in% paste(reference$seq_type, reference$genome), , drop = FALSE]
    if(nrow(new)) message("reference lacks ", paste(new$seq_type, new$genome, collapse = ", "),
                          "; computing those medians from seg_data")
    if(!"n_segmentations" %in% colnames(reference)) reference$n_segmentations = NA_integer_
    ref = rbind(reference[, colnames(computed)], new)
  }
  ref_median = ref$segments_per_100mb[match(stratum, paste(ref$seq_type, ref$genome))]
  qc$pen_oversegmentation = w$oversegmentation * pmax(log2(qc$segments_per_100mb / ref_median), 0)
  qc$pen_short = w$short * qc$short_segment_fraction
  qc$pen_deep_loss = w$deep_loss * qc$deep_loss_fraction
  qc$pen_contradiction = w$contradiction * (qc$deep_loss_fraction > 0.01 &
                                              !is.na(qc$deep_loss_logratio) & qc$deep_loss_logratio > -0.5)
  pens = c("pen_coverage", "pen_oversegmentation", "pen_short", "pen_deep_loss", "pen_contradiction")
  qc$score = round(rowSums(qc[, pens]), 3)

  qc$eligible = TRUE
  if(!is.null(eligible)){
    e = as.data.frame(eligible)
    missing_cols = setdiff(c("ID", by, "eligible"), colnames(e))
    if(length(missing_cols)) stop("eligible lacks: ", paste(missing_cols, collapse = ", "))
    at = match(paste(qc$sample_id, qc$source), paste(e$ID, e[[by]]))
    qc$eligible = !(e$eligible[at] %in% FALSE)
  }

  # rank within each sample: eligible first, then score, prefer order, name
  pref = match(qc$source, prefer)
  pref[is.na(pref)] = length(prefer) + 1
  o = order(qc$sample_id, !qc$eligible, qc$score, pref, qc$source)
  qc = qc[o, , drop = FALSE]
  qc$rank = stats::ave(seq_len(nrow(qc)), qc$sample_id, FUN = seq_along)
  qc$selected = qc$rank == 1 & qc$eligible

  # sample_id duplicates ID so the table joins to metadata without naming columns
  out = data.frame(ID = qc$sample_id, sample_id = qc$sample_id, seq_type = qc$seq_type,
                   genome = qc$genome, source = qc$source,
                   n_segments = qc$n_segments,
                   covered_fraction = round(qc$covered_fraction, 3),
                   segments_per_100mb = round(qc$segments_per_100mb, 2),
                   short_segment_fraction = round(qc$short_segment_fraction, 3),
                   deep_loss_fraction = round(qc$deep_loss_fraction, 4),
                   deep_loss_logratio = round(qc$deep_loss_logratio, 3),
                   round(qc[, pens], 3), score = qc$score, eligible = qc$eligible, rank = qc$rank,
                   selected = qc$selected, stringsAsFactors = FALSE)
  names(out)[names(out) == "source"] = by
  rownames(out) = NULL
  rownames(ref) = NULL
  attr(out, "oversegmentation_reference") = ref
  out
}

# Native genome version of a genome_build value: the GRCh37 builds are one
# version, the GRCh38 builds another
.native_genome_version = function(genome_build){
  g = tolower(as.character(genome_build))
  ifelse(g %in% c("grch37", "hg19", "hs37d5"), "grch37",
         ifelse(g %in% c("hg38", "grch38"), "hg38", g))
}

