#' @title Find overlapping segments within samples
#'
#' @description Report copy number segments that overlap another segment of
#' the same sample (and seq_type) on the same chromosome. A sample's
#' segments should tile its genome without overlapping; overlaps usually
#' mean the data were duplicated (e.g. a table bound or joined to a copy of
#' itself) or that two profiles share a sample ID, and they cause double
#' counting and inconsistent results downstream. See
#' [resolve_segment_overlaps] to clean them up.
#'
#' @details Within each sample (by `by`) and chromosome, segments are sorted
#' by start and each is compared with the furthest end reached by the
#' segments before it, so a segment overlapping any earlier one is found.
#' Overlaps of up to `tolerance` bp are ignored (e.g. 1 for data where a
#' segment ends on the base the next one starts on).
#'
#' @param seg_data Segmented data with ID, chrom, start and end.
#' @param by Columns identifying one profile. Default: "ID", plus
#' "seg_seq_type" when present (as returned by `get_cn_segments` for samples
#' with more than one seq_type).
#' @param tolerance Overlaps of this many bp or fewer are ignored. Default 0.
#'
#' @return A data frame with one row per overlapping segment: the `by`
#' columns, chrom, its start and end, the sequence it shares with earlier
#' segments in bp (`overlap_bp`; its whole length if nested in one),
#' whether an earlier segment has exactly the same coordinates
#' (`same_coordinates`) and whether that earlier segment is an identical row
#' (`duplicate_row`). Zero rows when there are no overlaps.
#'
#' @export
#'
#' @examples
#' seg = data.frame(ID = "s1", chrom = "1", start = c(1, 101, 101),
#'                  end = c(100, 200, 200), CN = 2)
#' check_segment_overlaps(seg)
check_segment_overlaps = function(seg_data, by = NULL, tolerance = 0){
  s = as.data.frame(seg_data)
  if(!all(c("ID", "chrom", "start", "end") %in% colnames(s))){
    stop("seg_data needs columns ID, chrom, start and end")
  }
  if(is.null(by)){
    by = c("ID", intersect("seg_seq_type", colnames(s)))
  }
  if(!all(by %in% colnames(s))){
    stop("by columns not in seg_data: ", paste(setdiff(by, colnames(s)), collapse = ", "))
  }
  ov = .segment_overlap_rows(s, by, tolerance)
  if(length(ov$rows) == 0){
    out = s[0, c(by, "chrom", "start", "end"), drop = FALSE]
    out$overlap_bp = numeric(0)
    out$same_coordinates = logical(0)
    out$duplicate_row = logical(0)
    return(out)
  }
  out = s[ov$rows, c(by, "chrom", "start", "end"), drop = FALSE]
  out$overlap_bp = ov$overlap_bp
  out$same_coordinates = ov$same_coordinates
  out$duplicate_row = ov$duplicate_row
  rownames(out) = NULL
  out
}

#' @title Remove overlapping segments within samples
#'
#' @description Clean up segmented data in which segments of the same
#' sample overlap (see [check_segment_overlaps]), reporting what was
#' changed.
#'
#' @details Exact duplicate rows are always removed first (keeping one
#' copy); this fixes data that were bound or joined to a copy of
#' themselves. Remaining overlaps are handled by `method`:
#' \itemize{
#'   \item "trim" (default): segments with the same coordinates as an
#'     earlier one (but different values) are dropped, keeping the first,
#'     and a segment that partly overlaps earlier ones starts after the
#'     furthest end reached (keeping its own convention: one base after for
#'     1-based coordinates where segments don't share bases, i.e. when
#'     `one_based = TRUE`); segments left with no length are dropped.
#'   \item "drop_samples": samples (profiles) that still have overlapping
#'     segments are removed entirely.
#'   \item "error": stop if any overlaps remain.
#' }
#' Overlaps of up to `tolerance` bp are left alone. Attributes such as the
#' genome build of a `seg_data` object are kept.
#'
#' @param seg_data Segmented data with ID, chrom, start and end.
#' @param method "trim" (default), "drop_samples" or "error". See Details.
#' @param by,tolerance As in [check_segment_overlaps].
#' @param one_based Set to TRUE (default) when coordinates are 1-based and
#' inclusive (consecutive segments are end, end + 1), as in GAMBL seg files,
#' so trimmed segments start one base after the previous end; FALSE for
#' 0-based half-open coordinates.
#' @param verbose Report what was removed or changed. Default TRUE.
#'
#' @return `seg_data` without overlapping segments.
#'
#' @export
#'
#' @examples
#' seg = data.frame(ID = "s1", chrom = "1", start = c(1, 101, 101, 150),
#'                  end = c(100, 200, 200, 300), CN = c(2, 3, 3, 4))
#' resolve_segment_overlaps(seg)
resolve_segment_overlaps = function(seg_data, method = "trim", by = NULL,
                                    tolerance = 0, one_based = TRUE,
                                    verbose = TRUE){
  method = match.arg(method, c("trim", "drop_samples", "error"))
  s = as.data.frame(seg_data)
  if(is.null(by)){
    by = c("ID", intersect("seg_seq_type", colnames(s)))
  }
  n0 = nrow(s)
  # exact duplicate rows
  dup = duplicated(s)
  s = s[!dup, , drop = FALSE]
  if(verbose && any(dup)){
    message("removed ", sum(dup), " duplicate row(s)")
  }
  ov = .segment_overlap_rows(s, by, tolerance)
  if(length(ov$rows) > 0){
    profile = do.call(paste, c(s[by], sep = "\r"))
    affected = unique(profile[ov$rows])
    if(method == "error"){
      stop(length(ov$rows), " overlapping segment(s) in ", length(affected),
           " sample(s); see check_segment_overlaps()")
    }
    if(method == "drop_samples"){
      drop = profile %in% affected
      if(verbose){
        message("removed ", length(affected), " sample(s) with overlapping segments (",
                sum(drop), " rows)")
      }
      s = s[!drop, , drop = FALSE]
    }else{
      # trim, repeatedly, until no overlaps remain (trimming one segment can
      # leave it overlapping a later one)
      n_same = 0; n_trim = 0; n_empty = 0
      repeat{
        ov = .segment_overlap_rows(s, by, tolerance)
        if(length(ov$rows) == 0) break
        same = ov$rows[ov$same_coordinates]
        n_same = n_same + length(same)
        partial = setdiff(ov$rows, same)
        new_start = ov$reach[!ov$same_coordinates] + as.numeric(one_based)
        s$start[partial] = new_start
        n_trim = n_trim + length(partial)
        empty = partial[s$end[partial] < s$start[partial] + as.numeric(!one_based)]
        n_empty = n_empty + length(empty)
        s = s[-c(same, empty), , drop = FALSE]
      }
      if(verbose){
        message("dropped ", n_same, " segment(s) repeating an earlier segment's ",
                "coordinates, trimmed ", n_trim, " partly overlapping segment(s)",
                if(n_empty > 0) paste0(" (", n_empty, " left empty and dropped)") else "",
                " in ", length(affected), " sample(s)")
      }
    }
  }
  if(verbose && nrow(s) == n0){
    message("no overlapping segments")
  }
  rownames(s) = NULL
  # keep class and attributes (e.g. genome_build of seg_data)
  for(a in setdiff(names(attributes(seg_data)), c("names", "row.names"))){
    attr(s, a) = attr(seg_data, a)
  }
  s
}

# Rows (indices into s) of segments overlapping an earlier segment of the
# same profile and chromosome by more than tolerance bp, with the overlap,
# the furthest end reached before them and whether an earlier segment has
# the same coordinates or is an identical row.
.segment_overlap_rows = function(s, by, tolerance){
  n = nrow(s)
  if(n < 2) return(list(rows = integer(0)))
  key = do.call(paste, c(s[c(by, "chrom")], sep = "\r"))
  o = order(key, s$start, s$end)
  k = key[o]; st = as.numeric(s$start[o]); en = as.numeric(s$end[o])
  # furthest end reached by the earlier segments of the same profile/chrom
  first = c(TRUE, k[-1] != k[-n])
  reach = rep(-Inf, n)
  run_max = -Inf
  for(i in seq_len(n)){
    if(first[i]) run_max = -Inf
    reach[i] = run_max
    run_max = max(run_max, en[i])
  }
  # shared sequence: up to the earlier reach or this segment's end
  # (a nested segment overlaps along its whole length)
  overlap_bp = pmin(reach, en) - st
  hit = which(overlap_bp > tolerance)
  if(length(hit) == 0) return(list(rows = integer(0)))
  coord_key = paste(k, st, en, sep = "\r")
  same = duplicated(coord_key)[hit]
  row_key = do.call(paste, c(s[o, , drop = FALSE], sep = "\r"))
  dup_row = duplicated(row_key)[hit]
  list(rows = o[hit], overlap_bp = overlap_bp[hit], reach = reach[hit],
       same_coordinates = same, duplicate_row = dup_row)
}

#' @title Check loaded segments for overlaps
#'
#' @description Used by `get_cn_segments` (GAMBLR.results, GAMBLR.open) to
#' check that no sample (profile) has overlapping segments, and to warn,
#' stop or clean them up.
#'
#' @param seg_data Segmented data.
#' @param action "warn" (default): warn, naming the samples; "error": stop;
#' "resolve": clean up with [resolve_segment_overlaps] (method "trim");
#' "none": skip the check.
#' @param by Columns identifying one profile (see [check_segment_overlaps]).
#'
#' @return `seg_data`, cleaned up with `action = "resolve"`.
#'
#' @export
#' @keywords internal
handle_segment_overlaps = function(seg_data, action = "warn", by = NULL){
  action = match.arg(action, c("warn", "error", "resolve", "none"))
  if(action == "none" || is.null(seg_data) || nrow(seg_data) == 0){
    return(seg_data)
  }
  if(action == "resolve"){
    return(resolve_segment_overlaps(seg_data, method = "trim", by = by))
  }
  ov = check_segment_overlaps(seg_data, by = by)
  if(nrow(ov) > 0){
    ids = unique(ov$ID)
    text = paste0(nrow(ov), " overlapping segment(s) in ", length(ids),
                  " sample(s) (", sum(ov$duplicate_row), " duplicate rows): ",
                  paste(utils::head(ids, 10), collapse = ", "),
                  if(length(ids) > 10) ", ..." else "",
                  ". Segments of a sample should not overlap; see ",
                  "check_segment_overlaps() and resolve_segment_overlaps().")
    if(action == "error") stop(text, call. = FALSE)
    warning(text, call. = FALSE)
  }
  seg_data
}
