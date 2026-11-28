#' @title Keep each sample's best segmentation
#'
#' @description From segments of several sources (e.g. every pipeline's
#' segments from `GAMBLR.results::get_cn_segments(..., pipeline = "all")`),
#' keep only each sample's chosen source: the one [compare_segmentations]
#' selects, or the one named in a saved selection.
#'
#' @details Without `selection`, the sources are scored here with
#' [compare_segmentations], so give `seg_data` and `these_samples_metadata`
#' for the whole cohort or project (see there). The choice depends on the
#' sources available: a new pipeline can change it. To make it repeatable,
#' save the scores (`attr(result, "selection")`, e.g. with
#' `readr::write_tsv`) and pass them back as `selection` next time; each
#' sample then keeps the source chosen before.
#'
#' @param seg_data Segments with ID, chrom, start, end, and the column named
#' in `by`.
#' @param these_samples_metadata Metadata with sample_id (and seq_type and
#' genome_build, for scoring). Segments of samples not in it are dropped, with
#' a message. Default NULL: all samples in `seg_data`.
#' @param selection Optional saved choice: the table from
#' [compare_segmentations] (its `selected` rows are used), or any table with
#' sample_id (or ID) and the `by` column, and optionally seq_type, with one
#' row per sample. Default NULL: score `seg_data` now.
#' @param by Column naming each segment's source. Default "pipeline".
#' @param ... Further arguments for [compare_segmentations] (e.g. `weights`,
#' `prefer`, `reference`) when scoring here.
#'
#' @return `seg_data`'s rows for each sample's chosen source (samples without
#' a choice, or without segments from their chosen source, are dropped, with
#' a message). The selection used, one row per sample, is in
#' `attr(result, "selection")`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta = GAMBLR.results::get_gambl_metadata() %>% dplyr::filter(cohort == "dlbcl_chapuy")
#' segs = GAMBLR.results::get_cn_segments(meta, pipeline = "all")
#' best = select_cn_segments(segs, meta)
#'
#' # save the choice and reuse it, so later runs keep it
#' readr::write_tsv(attr(best, "selection"), "cn_selection.tsv")
#' best = select_cn_segments(segs, meta, selection = readr::read_tsv("cn_selection.tsv"))
#' }
select_cn_segments = function(seg_data,
                              these_samples_metadata = NULL,
                              selection = NULL,
                              by = "pipeline",
                              ...){
  genome_build = attr(seg_data, "genome_build")
  segs = as.data.frame(seg_data)
  if(!by %in% colnames(segs)) stop("seg_data has no '", by, "' column")
  if(!is.null(these_samples_metadata)){
    outside = !segs$ID %in% these_samples_metadata$sample_id
    if(any(outside)){
      message("dropping ", length(unique(segs$ID[outside])), " sample(s) not in these_samples_metadata")
      segs = segs[!outside, , drop = FALSE]
    }
  }
  if(is.null(selection)){
    selection = compare_segmentations(segs, these_samples_metadata, by = by, ...)
  }

  # one chosen source per sample (and seq_type, where the selection has it)
  sel = as.data.frame(selection)
  if(!"sample_id" %in% colnames(sel) && "ID" %in% colnames(sel)) sel$sample_id = sel$ID
  missing_cols = setdiff(c("sample_id", by), colnames(sel))
  if(length(missing_cols)) stop("selection lacks: ", paste(missing_cols, collapse = ", "))
  if("selected" %in% colnames(sel)) sel = sel[sel$selected %in% TRUE, , drop = FALSE]
  if(!"seq_type" %in% colnames(sel)) sel$seq_type = NA_character_
  sel$seq_type = as.character(sel$seq_type)
  if(any(duplicated(sel[, c("sample_id", "seq_type")]))){
    stop("selection has more than one ", by, " for: ",
         paste(utils::head(unique(sel$sample_id[duplicated(sel[, c("sample_id", "seq_type")])]), 10),
               collapse = ", "))
  }

  # match on sample_id and seq_type, or sample_id alone where the selection
  # or the segments have no seq_type
  seg_type = if("seg_seq_type" %in% colnames(segs)) as.character(segs$seg_seq_type) else rep(NA_character_, nrow(segs))
  typed = sel[!is.na(sel$seq_type), , drop = FALSE]
  untyped = sel[is.na(sel$seq_type), , drop = FALSE]
  want = typed[[by]][match(paste(segs$ID, seg_type), paste(typed$sample_id, typed$seq_type))]
  by_id = untyped[[by]][match(segs$ID, untyped$sample_id)]
  if(all(is.na(seg_type))) by_id = ifelse(is.na(by_id), typed[[by]][match(segs$ID, typed$sample_id)], by_id)
  want = ifelse(is.na(want), by_id, want)

  no_choice = setdiff(unique(segs$ID), sel$sample_id)
  if(length(no_choice)) message("dropping ", length(no_choice), " sample(s) with no choice in selection")
  out = segs[!is.na(want) & as.character(segs[[by]]) == want, , drop = FALSE]
  absent = setdiff(intersect(unique(segs$ID), sel$sample_id), out$ID)
  if(length(absent)){
    message(length(absent), " sample(s) have no segments from their chosen ", by, ": ",
            paste(utils::head(absent, 10), collapse = ", "))
  }
  rownames(out) = NULL
  if(!is.null(genome_build) && genome_build %in% c("grch37", "hg38")) out = create_seg_data(out, genome_build)
  attr(out, "selection") = sel[sel$sample_id %in% unique(segs$ID), , drop = FALSE]
  out
}
