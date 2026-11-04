#' @title Copy number type of each segmentation pipeline
#'
#' @description Whether each pipeline's CN column is absolute copy number
#' (purity- and ploidy-corrected, e.g. battenberg, PureCN, sequenza and the
#' published DLBClass/ABSOLUTE segments) or relative, i.e. only a rescaled
#' log ratio (e.g. CNVkit and gatk_cnv, whose CN is `round(2 *
#' 2^log.ratio)`). Copy number features are called on CN for absolute
#' pipelines and on log ratios for relative ones (see
#' [cnv_calling_settings]).
#'
#' @return A data frame: pipeline, cn_type ("absolute" or "relative").
#' Pipelines not listed are treated as relative.
#'
#' @export
#'
#' @examples
#' cn_pipeline_types()
cn_pipeline_types = function(){
  data.frame(pipeline = c("battenberg", "purecn", "sequenza", "controlfreec", "dlbclass",
                          "cnvkit", "gatk_cnv", "schmitz_array"),
             cn_type = c(rep("absolute", 5), rep("relative", 3)),
             stringsAsFactors = FALSE)
}

# Pipelines whose segments come from arrays (unfragmented, like the segments
# LymphGen's training data were built from)
.array_cn_pipelines = c("schmitz_array")

#' @title How to call copy number features from each pipeline's segments
#'
#' @description The settings for [segmented_data_to_cnv_status] that suit
#' each pipeline and feature set: absolute pipelines (see
#' [cn_pipeline_types]) are called on CN, centred on each sample's dominant
#' copy number ("mode"); relative ones on log ratios, centred on the sample's
#' median. DLBClass features are always called on log ratios, as DLBClass
#' defines them. Thresholds follow the feature set: LymphGen's copy number
#' categories are 0.5 and 1.5 copies on either scale (as
#' `GAMBLR.predict::lymphgen_input()` uses them); otherwise the defaults of the scale (0.1
#' and 0.9 for log ratios, as in DLBClass, and 0.5 and 1.5 for CN). LymphGen
#' sets also get, on array segments ("schmitz_array"), `gene_rule =
#' "lymphgen"` (LymphGen's own way of calling gene features, see
#' [segmented_data_to_cnv_status]), which reproduces its training data from
#' them; sequencing pipelines keep `gene_rule = "fraction"` (the length of a
#' change measured over its run of neighbouring segments), which matches
#' LymphGen's classifications best on exomes, whose segmentation splits broad
#' changes into pieces under 30 Mb. On the CN scale, LymphGen sets get level 2
#' relative to each sample's baseline (`homdel_max_cn` NA, i.e. NULL): a
#' change of at least `high_thresh` copies, as in LymphGen's training data,
#' where e.g. a loss from 4 to 2 copies in a tetraploid tumour is a HOMDEL.
#' Other sets keep `homdel_max_cn = 0.5` (a HOMDEL is CN 0.5 or less).
#'
#' @param pipeline Pipeline names, e.g. the `pipeline` column of
#' `GAMBLR.results::get_cn_segments(..., pipeline = "all")`.
#' @param feature_set Optional feature set ("LymphGen_model", "LymphGen",
#' "DLBClass", ...); both LymphGen sets get LymphGen's thresholds.
#' @param cnv_args Optional named list of [segmented_data_to_cnv_status]
#' arguments; its `value_scale`, `centering`, `single_thresh`,
#' `high_thresh`, `gene_rule` and `homdel_max_cn` override the settings for
#' every pipeline.
#' @param cn_type Optional CN type ("absolute" or "relative") for each
#' element of `pipeline`, overriding [cn_pipeline_types], e.g. the `cn_type`
#' column [calibrate_relative_cn] adds. `pipeline` is then not deduplicated.
#'
#' @return A data frame with one row per pipeline: pipeline, cn_type,
#' value_scale, centering, single_thresh, high_thresh, gene_rule,
#' homdel_max_cn (NA for relative level 2, i.e. NULL in
#' [segmented_data_to_cnv_status]).
#'
#' @export
#'
#' @examples
#' cnv_calling_settings(c("battenberg", "gatk_cnv"), feature_set = "LymphGen")
cnv_calling_settings = function(pipeline, feature_set = NULL, cnv_args = list(), cn_type = NULL){
  pipeline = as.character(pipeline)
  types = cn_pipeline_types()
  if(is.null(cn_type)){
    pipeline = unique(pipeline)
    cn_type = types$cn_type[match(pipeline, types$pipeline)]
  }else{
    if(length(cn_type) != length(pipeline)) stop("cn_type must have one value per pipeline")
    # NA entries fall back to the pipeline's type
    cn_type = as.character(cn_type)
    cn_type = ifelse(is.na(cn_type), types$cn_type[match(pipeline, types$pipeline)], cn_type)
  }
  cn_type[is.na(cn_type)] = "relative"
  value_scale = ifelse(cn_type == "absolute", "CN", "log_ratio")
  # DLBClass's copy number features are defined on log ratios (its seg2gsm
  # uses segment means), whatever the pipeline
  if(!is.null(feature_set) && identical(tolower(feature_set), "dlbclass")) value_scale[] = "log_ratio"
  if(!is.null(cnv_args$value_scale)) value_scale = rep(cnv_args$value_scale, length(pipeline))
  centering = ifelse(value_scale == "CN", "mode", "median")
  if(!is.null(cnv_args$centering)) centering = rep(cnv_args$centering, length(pipeline))
  # LymphGen's thresholds and gene rule for both of its sets ("LymphGen", "LymphGen_model")
  lymphgen = !is.null(feature_set) && grepl("^lymphgen", tolower(feature_set))
  thresholds = lapply(value_scale, function(scale){
    defaults = if(lymphgen) c(single = 0.5, high = 1.5) else .resolve_cnv_thresholds(scale, NULL, NULL)
    c(single = if(is.null(cnv_args$single_thresh)) defaults[["single"]] else cnv_args$single_thresh,
      high = if(is.null(cnv_args$high_thresh)) defaults[["high"]] else cnv_args$high_thresh)
  })
  data.frame(pipeline = pipeline, cn_type = cn_type, value_scale = value_scale,
             centering = centering,
             single_thresh = vapply(thresholds, `[[`, 0, "single"),
             high_thresh = vapply(thresholds, `[[`, 0, "high"),
             # LymphGen's own gene rule (each segment on its own length) for
             # array segments; for sequencing, whose segmentation splits broad
             # changes into pieces under 30 Mb, lengths over runs of segments
             gene_rule = if(!is.null(cnv_args$gene_rule)) cnv_args$gene_rule else
               ifelse(lymphgen & pipeline %in% .array_cn_pipelines, "lymphgen", "fraction"),
             # level 2 on the CN scale: relative to the baseline for LymphGen
             # (NA, i.e. homdel_max_cn = NULL), absolute CN <= 0.5 otherwise
             homdel_max_cn = if("homdel_max_cn" %in% names(cnv_args)) {
               if(is.null(cnv_args$homdel_max_cn)) NA_real_ else cnv_args$homdel_max_cn
             } else if(lymphgen) NA_real_ else 0.5,
             stringsAsFactors = FALSE)
}

#' @title Each sample's copy number baseline
#'
#' @description The baseline [segmented_data_to_cnv_status] measures copy
#' number changes from, per sample: the length-weighted median log ratio
#' ("median", `value_scale = "log_ratio"`), the dominant copy number ("mode",
#' CN only), the length-weighted median CN ("median", CN), or diploid
#' ("none").
#'
#' @param seg_data Segments with ID, chrom, start, end and log.ratio or CN.
#' @param value_scale "log_ratio" or "CN".
#' @param centering "median", "mode" (CN only) or "none".
#' @param autosomes_only Use autosomes only. Default TRUE.
#'
#' @return A named vector of baselines (names = sample IDs), on the scale of
#' `value_scale`.
#'
#' @export
#'
#' @examples
#' segs = data.frame(ID = "s1", chrom = c("1", "2"), start = 1, end = c(2e8, 1e8),
#'                   CN = c(2, 1), log.ratio = c(0, -1))
#' cnv_sample_baseline(segs, "CN", "mode")
cnv_sample_baseline = function(seg_data, value_scale = "log_ratio", centering = "median",
                               autosomes_only = TRUE){
  segs = as.data.frame(seg_data)
  value_col = if(value_scale == "CN") "CN" else "log.ratio"
  if(!value_col %in% colnames(segs)) stop("seg_data has no ", value_col, " column")
  segs$chrom = gsub("^chr", "", as.character(segs$chrom))
  segs = segs[!is.na(segs[[value_col]]), , drop = FALSE]
  if(nrow(segs) == 0) return(stats::setNames(numeric(0), character(0)))
  centre = .cnv_sample_centre(segs, segs[[value_col]], segs$end - segs$start, value_scale,
                              centering, autosomes_only)
  first = !duplicated(segs$ID)
  stats::setNames(centre[first], segs$ID[first])
}
