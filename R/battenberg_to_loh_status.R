#' @title Call LOH at CNV features from Battenberg results
#'
#' @description Call copy-neutral LOH, LOH with copy loss (deletion), or any
#' complete LOH at each feature in a feature table, for every sample in a
#' [summarize_battenberg] result. Uses the same calling rules as
#' [segmented_data_to_cnv_status], so the result can be combined with gain
#' and loss calls, e.g. for [test_feature_associations] or
#' `GAMBLR.viz::stacked_CN_heatmap`.
#'
#' @details Each assessable segment (see [summarize_battenberg]) is given a
#' copy number of 1 if it belongs to an LOH event of the requested type and 2
#' otherwise; unassessable segments (e.g. subclonal segments whose states
#' disagree) are left out, so they reduce a feature's coverage rather than
#' counting as no LOH. These are passed to [segmented_data_to_cnv_status]
#' with `value_scale = "CN"`, `centering = "none"` and `background = "none"`
#' (measuring focal features against their arm would hide focal LOH inside
#' arm-level LOH). A feature is therefore called when the covered part of it
#' is at least `min_coverage` of its length and LOH of the requested type
#' makes up at least its `min_fraction` of the covered length (defaults 0.5
#' and 0.5). Features with too little coverage are 0, like features without
#' LOH (use `encoding = "value"` to see them as NA).
#'
#' LOH types, from the `loh_type` of [summarize_battenberg]'s
#' `complete_loh_events` (relative to `cnloh_reference`):
#' \itemize{
#'   \item "cnloh": copy-neutral LOH (`loh_type` neutral).
#'   \item "deletion": LOH with copy loss (`loh_type` loss).
#'   \item "any": any complete LOH, including LOH with gain.
#' }
#' Only events of at least `min_event_size` are used.
#'
#' The default features are the loss-direction features of
#' [get_cnv_features] plus every autosomal chromosome arm (except the
#' acrocentric p arms) that isn't already one of them, since LOH usually acts
#' at deletion loci and arm-level LOH is common. Feature names end in
#' `_cnloh`, `_deletion_loh` or `_loh` instead of `_loss`, so they can be
#' combined with gain and loss calls (`cbind` for the matrix format).
#'
#' Only grch37 is supported, as for [segmented_data_to_cnv_status].
#'
#' @param battenberg Output of [summarize_battenberg] (with
#' `complete_loh_events` and `segments`).
#' @param loh_type "cnloh" (default), "deletion" or "any". See Details.
#' @param cnv_features Optional feature table (see [get_cnv_features]).
#' Default: the loss features plus all arms (see Details). Its direction is
#' set to "loss" for every row.
#' @param exclude_whole_chromosome Set to TRUE to leave out LOH from
#' whole-chromosome events (e.g. whole-chromosome uniparental disomy), so
#' arm and focal features reflect segmental LOH only. Default FALSE.
#' @param min_event_size Minimum event size (bp). Default 1e6.
#' @param genome_build Genome build of the Battenberg results. Default
#' "grch37" (the only one supported).
#' @param these_samples_metadata Optional metadata; passed to
#' [segmented_data_to_cnv_status] (samples are restricted to and ordered by
#' its sample_id).
#' @param ... Other arguments passed to [segmented_data_to_cnv_status], e.g.
#' `encoding`, `output_format` ("matrix" or "long") or `gene_labels`.
#'
#' @return As [segmented_data_to_cnv_status], with renamed features.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bb = summarize_battenberg(battenberg_dir = "battenberg_results")
#' cnloh = battenberg_to_loh_status(bb, encoding = "binary")
#' deletion_loh = battenberg_to_loh_status(bb, loh_type = "deletion",
#'                                         encoding = "binary")
#' # segmental copy-neutral LOH only
#' cnloh_segmental = battenberg_to_loh_status(bb, exclude_whole_chromosome = TRUE,
#'                                            encoding = "binary")
#' }
battenberg_to_loh_status = function(battenberg,
                                    loh_type = "cnloh",
                                    cnv_features = NULL,
                                    exclude_whole_chromosome = FALSE,
                                    min_event_size = 1e6,
                                    genome_build = "grch37",
                                    these_samples_metadata = NULL,
                                    ...){
  loh_type = match.arg(loh_type, c("cnloh", "deletion", "any"))
  if(!is.list(battenberg) || is.null(battenberg$segments) ||
     !"complete_loh_events" %in% names(battenberg)){
    stop("battenberg must be the output of summarize_battenberg() (with complete_loh_events and segments)")
  }
  extra = list(...)
  fixed = c("value_scale", "centering", "background", "single_thresh",
            "high_thresh", "seg_data")
  if(any(names(extra) %in% fixed)){
    stop(paste("these are set by battenberg_to_loh_status:",
               paste(intersect(names(extra), fixed), collapse = ", ")))
  }
  if(is.null(cnv_features)){
    cnv_features = .loh_default_features(genome_build)
  }
  cnv_features = as.data.frame(cnv_features)
  cnv_features$direction = "loss"
  cnv_features$background_arm = NA_character_

  # events of the requested type
  events = battenberg$complete_loh_events
  if(!is.null(events) && nrow(events) > 0){
    keep = switch(loh_type,
                  cnloh = events$loh_type == "neutral",
                  deletion = events$loh_type == "loss",
                  any = rep(TRUE, nrow(events)))
    keep = keep & events$size >= min_event_size
    if(exclude_whole_chromosome){
      keep = keep & events$extent != "whole_chromosome"
    }
    events = events[keep, , drop = FALSE]
  }

  # assessable segments, CN 1 if in one of those events, otherwise 2
  seg = battenberg$segments
  seg = seg[seg$assessable, , drop = FALSE]
  in_event = rep(FALSE, nrow(seg))
  if(!is.null(events) && nrow(events) > 0){
    loh_seg = which(!is.na(seg$homologues) & seg$homologues == "one" &
                      !is.na(seg$total_cn))
    for(i in seq_len(nrow(events))){
      e = events[i, ]
      hit = loh_seg[seg$sample[loh_seg] == e$sample & seg$chrom[loh_seg] == e$chrom &
                      seg$start[loh_seg] >= e$start & seg$end[loh_seg] <= e$end &
                      seg$total_cn[loh_seg] == e$total_cn]
      in_event[hit] = TRUE
    }
  }
  seg_data = data.frame(ID = seg$sample, chrom = seg$chrom, start = seg$start,
                        end = seg$end, CN = ifelse(in_event, 1, 2),
                        stringsAsFactors = FALSE)
  args = c(list(seg_data = seg_data, cnv_features = cnv_features,
                genome_build = genome_build, value_scale = "CN",
                centering = "none", background = "none",
                single_thresh = 0.5, high_thresh = 1.5),
           extra)
  if(!is.null(these_samples_metadata)){
    args$these_samples_metadata = these_samples_metadata
  }
  out = do.call(segmented_data_to_cnv_status, args)

  suffix = switch(loh_type, cnloh = "_cnloh", deletion = "_deletion_loh", any = "_loh")
  rename = function(v) sub("_loss$", suffix, v)
  if(identical(extra$output_format, "long")){
    out$name = rename(out$name)
  }else{
    colnames(out) = rename(colnames(out))
  }
  out
}

# Loss features from get_cnv_features() plus every autosomal arm (except
# acrocentric p arms) that isn't already a loss feature
.loh_default_features = function(genome_build){
  features = as.data.frame(get_cnv_features(genome_build))
  features = features[features$direction == "loss", , drop = FALSE]
  arms = as.data.frame(GAMBLR.data::chromosome_arms_grch37)
  arms$chromosome = gsub("^chr", "", as.character(arms$chromosome))
  arms = arms[arms$chromosome %in% as.character(1:22) &
                !(arms$chromosome %in% c("13", "14", "15", "21", "22") &
                    arms$arm == "p"), ]
  arm_names = paste0(arms$chromosome, arms$arm)
  arms = arms[!arm_names %in% features$region_name[features$type == "arm"], ]
  arm_names = paste0(arms$chromosome, arms$arm)
  if(nrow(arms) > 0){
    new = features[rep(NA_integer_, nrow(arms)), , drop = FALSE]
    new$name = paste0(arm_names, "_loss")
    new$type = "arm"
    new$region_name = arm_names
    new$chrom = arms$chromosome
    new$start = arms$start
    new$end = arms$end
    new$direction = "loss"
    new$min_fraction = 0.5
    new$min_coverage = 0.5
    features = rbind(features, new)
  }
  features = features[order(as.integer(features$chrom), features$start), ]
  rownames(features) = NULL
  features
}
