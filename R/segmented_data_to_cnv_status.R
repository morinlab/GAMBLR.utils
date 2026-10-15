# CNV feature calling (arm-level and focal) from segmented data.
#
# TODO(hg38): only grch37/hg19 is supported. The bundled feature table is
# grch37-only, and the guard in .check_cnv_genome_build() refuses anything
# else. Proper hg38 support needs lifted-over focal features plus hg38 GISTIC
# arm boundaries (see data-raw/cnv_features.R).

# grch37 chromosome lengths, used to catch hg38 segments passed as grch37
.grch37_chrom_lengths = c(
  "1" = 249250621, "2" = 243199373, "3" = 198022430, "4" = 191154276,
  "5" = 180915260, "6" = 171115067, "7" = 159138663, "8" = 146364022,
  "9" = 141213431, "10" = 135534747, "11" = 135006516, "12" = 133851895,
  "13" = 115169878, "14" = 107349540, "15" = 102531392, "16" = 90354753,
  "17" = 81195210, "18" = 78077248, "19" = 59128983, "20" = 63025520,
  "21" = 48129895, "22" = 51304566, "X" = 155270560, "Y" = 59373566
)

.check_cnv_genome_build = function(seg_data, genome_build){
  attr_build = get_genome_build(seg_data)
  if(missing(genome_build) || is.null(genome_build)){
    if(is.null(attr_build)){
      stop("genome_build must be specified (only \"grch37\" is currently supported)")
    }
    genome_build = attr_build
  }
  if(!genome_build %in% c("grch37", "hg19")){
    stop(paste0(
      "CNV feature calling currently only supports grch37/hg19 (got \"",
      genome_build, "\"). The bundled feature coordinates are grch37-only; ",
      "hg38 support has not been implemented yet."
    ))
  }
  if(!is.null(attr_build) && !attr_build %in% c("grch37", "hg19")){
    stop(paste0("seg_data is labelled as \"", attr_build,
                "\" but genome_build is \"", genome_build, "\""))
  }
  # Segments running past a grch37 chromosome end are a strong sign of hg38
  # data (several hg38 chromosomes are longer). This can't catch every case,
  # since hg38 data within grch37 bounds looks identical.
  chrom_len = .grch37_chrom_lengths[seg_data$chrom]
  beyond = which(!is.na(chrom_len) & seg_data$end > chrom_len + 1e5)
  if(length(beyond) > 0){
    example = seg_data[beyond[1], ]
    affected = unique(seg_data$ID[beyond])
    stop(paste0(
      length(beyond), " segment(s) from ", length(affected), " sample(s) ",
      "extend beyond the grch37 chromosome end (e.g. ", example$ID, " chr",
      example$chrom, ":", format(example$start, scientific = FALSE), "-",
      format(example$end, scientific = FALSE), " > ",
      format(chrom_len[beyond[1]], scientific = FALSE),
      "). These look like hg38 coordinates, which are ",
      "not supported yet. Affected samples: ",
      paste(utils::head(affected, 10), collapse = ", "),
      ifelse(length(affected) > 10, ", ...", "")
    ))
  }
  invisible(genome_build)
}

# Thresholds default to values suited to the scale: DLBClass's 0.1/0.9 for
# copy number change from log ratios (small changes in impure tumours), and
# 0.5/1.5 for CN, where a whole-copy change is level 1 and two or more
# copies is level 2.
.resolve_cnv_thresholds = function(value_scale, single_thresh, high_thresh){
  defaults = switch(value_scale,
                    log_ratio = c(single = 0.1, high = 0.9),
                    CN = c(single = 0.5, high = 1.5))
  if(is.null(single_thresh)){
    single_thresh = defaults[["single"]]
  }
  if(is.null(high_thresh)){
    high_thresh = defaults[["high"]]
  }
  if(!is.numeric(single_thresh) || length(single_thresh) != 1 ||
     !is.numeric(high_thresh) || length(high_thresh) != 1){
    stop("single_thresh and high_thresh must be single numbers")
  }
  if(single_thresh < 0 || high_thresh < single_thresh){
    stop("thresholds must satisfy 0 <= single_thresh <= high_thresh")
  }
  c(single = single_thresh, high = high_thresh)
}

# DLBClass's calc_region_median generalised to any length fraction: the value
# reached once `fraction` of the total length has been accumulated, starting
# from the most extreme values in the event's direction (largest for gains,
# smallest for losses). fraction = 0.5 is the length-weighted median.
.weighted_region_value = function(values, lengths, fraction, from_top){
  # numeric, since genome-wide integer lengths overflow cumsum
  lengths = as.numeric(lengths)
  keep = lengths > 0
  values = values[keep]
  lengths = lengths[keep]
  if(length(values) == 0){
    return(NA_real_)
  }
  o = order(values, decreasing = from_top)
  cumulative = cumsum(lengths[o])
  values[o][which(cumulative >= fraction * sum(lengths))[1]]
}

.validate_cnv_features = function(cnv_features){
  required = c("name", "type", "region_name", "chrom", "start", "end",
               "direction", "min_fraction", "min_coverage")
  missing_cols = required[!required %in% colnames(cnv_features)]
  if(length(missing_cols) > 0){
    stop(paste("cnv_features is missing columns:",
               paste(missing_cols, collapse = ", ")))
  }
  if(!all(cnv_features$direction %in% c("gain", "loss"))){
    stop("cnv_features$direction must be \"gain\" or \"loss\"")
  }
  if(any(duplicated(cnv_features$name))){
    stop("cnv_features$name must be unique")
  }
  if(!"background_arm" %in% colnames(cnv_features)){
    cnv_features$background_arm = NA_character_
  }
  cnv_features$background_arm[cnv_features$background_arm %in% ""] = NA
  arms = cnv_features$region_name[cnv_features$type == "arm"]
  bg = unique(stats::na.omit(cnv_features$background_arm))
  bg = bg[bg != "none"]
  if(any(!bg %in% arms)){
    stop(paste("background_arm values must be \"none\" or match an arm-type",
               "region_name in cnv_features:",
               paste(bg[!bg %in% arms], collapse = ", ")))
  }
  cnv_features$chrom = gsub("^chr", "", as.character(cnv_features$chrom))
  as.data.frame(cnv_features)
}

# Feature names with the representative gene inserted before the direction,
# e.g. "2p16.1_gain" -> "2p16.1~REL_gain". Features without a gene, or whose
# name already starts with that gene (e.g. "CDKN2A_loss"), are unchanged.
.label_cnv_feature_names = function(cnv_features, gene_labels, gene_separator){
  if(!gene_labels){
    return(cnv_features$name)
  }
  if(!is.character(gene_separator) || length(gene_separator) != 1 ||
     gene_separator == ""){
    stop("gene_separator must be a single non-empty string")
  }
  gene = cnv_features$representative_gene
  if(is.null(gene)){
    return(cnv_features$name)
  }
  label = !is.na(gene) & gene != "" &
    !startsWith(cnv_features$name, paste0(gene, "_"))
  labelled = cnv_features$name
  has_direction = grepl("_(gain|loss)$", labelled)
  labelled[label & has_direction] = sub(
    "_(gain|loss)$", "", labelled[label & has_direction]
  )
  labelled[label] = paste0(labelled[label], gene_separator, gene[label],
                           ifelse(has_direction[label],
                                  paste0("_", cnv_features$direction[label]), ""))
  if(any(duplicated(labelled))){
    stop("gene labels produce duplicate feature names")
  }
  labelled
}

# The arm each feature is measured relative to (NA = the genome baseline).
# "auto" applies the DLBClass rule to the table in use: a non-arm feature is
# measured relative to an arm feature in the same table that contains it and
# has the same direction. A background_arm value overrides the rule for that
# row ("none" = always genome-relative). "column" uses only background_arm.
.resolve_cnv_backgrounds = function(cnv_features, background){
  if(!background %in% c("auto", "column", "none")){
    stop("background must be \"auto\", \"column\" or \"none\"")
  }
  resolved = rep(NA_character_, nrow(cnv_features))
  if(background == "none"){
    return(resolved)
  }
  if(background == "auto"){
    arm_rows = which(cnv_features$type == "arm")
    midpoint = (cnv_features$start + cnv_features$end) / 2
    for(i in which(cnv_features$type != "arm")){
      parent = arm_rows[cnv_features$chrom[arm_rows] == cnv_features$chrom[i] &
                          cnv_features$direction[arm_rows] == cnv_features$direction[i] &
                          cnv_features$start[arm_rows] <= midpoint[i] &
                          cnv_features$end[arm_rows] >= midpoint[i]]
      if(length(parent) > 0){
        # the most specific arm feature if several contain it
        parent_len = cnv_features$end[parent] - cnv_features$start[parent]
        resolved[i] = cnv_features$region_name[parent[which.min(parent_len)]]
      }
    }
  }
  explicit = !is.na(cnv_features$background_arm)
  resolved[explicit] = cnv_features$background_arm[explicit]
  resolved[resolved %in% "none"] = NA
  resolved
}

# Trim or drop segments overlapping blacklisted regions, considering only
# blacklist entries that overlap a feature by more than
# min_feature_overlap of that feature's length. Entries are applied in order,
# as in DLBClass.
.apply_cnv_blacklist = function(segs, blacklist, cnv_features,
                                min_feature_overlap, drop_fraction){
  blacklist = as.data.frame(blacklist)
  # accept common spellings, e.g. the GISTIC/DLBClass Chromosome/Start/End
  aliases = list(chrom = c("chrom", "chr", "chromosome"),
                 start = "start", end = "end")
  for(col in names(aliases)){
    if(!col %in% colnames(blacklist)){
      alt = which(tolower(colnames(blacklist)) %in% aliases[[col]])
      if(length(alt) == 0){
        stop("blacklist must have chrom, start and end columns")
      }
      colnames(blacklist)[alt[1]] = col
    }
  }
  blacklist$chrom = gsub("^chr", "", as.character(blacklist$chrom))
  blacklist$chrom[blacklist$chrom == "23"] = "X"
  blacklist$chrom[blacklist$chrom == "24"] = "Y"
  bl_len = blacklist$end - blacklist$start
  relevant = rep(FALSE, nrow(blacklist))
  for(i in seq_len(nrow(cnv_features))){
    f = cnv_features[i, ]
    relevant = relevant |
      (blacklist$chrom == f$chrom &
         blacklist$start <= f$end & blacklist$end >= f$start &
         bl_len > (f$end - f$start) * min_feature_overlap)
  }
  blacklist = blacklist[relevant, , drop = FALSE]

  removed = rep(FALSE, nrow(segs))
  for(i in seq_len(nrow(blacklist))){
    b_start = blacklist$start[i]
    b_end = blacklist$end[i]
    hit = which(!removed & segs$chrom == blacklist$chrom[i] &
                  segs$start <= b_end & segs$end >= b_start)
    for(j in hit){
      s_start = segs$start[j]
      s_end = segs$end[j]
      if(s_start < b_start && s_end > b_end){
        # blacklisted region inside the segment: only drop the segment when
        # the region makes up most of it
        if((b_end - b_start) / (s_end - s_start) > drop_fraction){
          removed[j] = TRUE
        }
      }else if(s_start < b_start && s_end > b_start){
        segs$end[j] = b_start
      }else if(s_start < b_end && s_end > b_end){
        segs$start[j] = b_end
      }else if(s_start >= b_start && s_end <= b_end){
        removed[j] = TRUE
      }
    }
  }
  segs[!removed, , drop = FALSE]
}

# Shared preprocessing for annotate_cnv and segmented_data_to_cnv_status.
# Returns the (blacklist-trimmed) segments with .value, the centred copy
# number change relative to the sample baseline, plus per-sample arm levels
# for features called relative to their arm.
.prepare_cnv_segments = function(seg_data, cnv_features, genome_build,
                                 these_samples_metadata, value_scale,
                                 centering, centering_autosomes_only,
                                 blacklist, blacklist_min_feature_overlap,
                                 blacklist_drop_fraction, background){
  if(!value_scale %in% c("log_ratio", "CN")){
    stop("value_scale must be \"log_ratio\" or \"CN\"")
  }
  if(!centering %in% c("median", "mode", "none")){
    stop("centering must be \"median\", \"mode\" or \"none\"")
  }
  if(centering == "mode" && value_scale != "CN"){
    stop("centering = \"mode\" requires value_scale = \"CN\"")
  }
  value_col = ifelse(value_scale == "log_ratio", "log.ratio", "CN")
  required_cols = c("ID", "chrom", "start", "end", value_col)
  if(!all(required_cols %in% colnames(seg_data))){
    stop(paste("seg_data must have columns:",
               paste(required_cols, collapse = ", ")))
  }
  segs = as.data.frame(seg_data)
  segs$chrom = gsub("^chr", "", as.character(segs$chrom))
  genome_build = .check_cnv_genome_build(
    structure(segs, genome_build = get_genome_build(seg_data)), genome_build
  )
  if(is.null(cnv_features)){
    cnv_features = get_cnv_features(genome_build)
  }
  cnv_features = .validate_cnv_features(cnv_features)
  cnv_features$.background = .resolve_cnv_backgrounds(cnv_features, background)

  if(!missing(these_samples_metadata) && !is.null(these_samples_metadata)){
    all_samples = unique(these_samples_metadata$sample_id)
    segs = segs[segs$ID %in% all_samples, , drop = FALSE]
  }else{
    all_samples = unique(segs$ID)
  }
  all_samples = sort(all_samples)
  absent = all_samples[!all_samples %in% segs$ID]
  if(length(absent) > 0){
    warning(paste0(
      length(absent), " sample(s) have no segments and will have no CNV ",
      "events called: ", paste(utils::head(absent, 10), collapse = ", "),
      ifelse(length(absent) > 10, ", ...", "")
    ))
  }
  segs = segs[!is.na(segs[[value_col]]), , drop = FALSE]

  if(!is.null(blacklist)){
    segs = .apply_cnv_blacklist(segs, blacklist, cnv_features,
                                blacklist_min_feature_overlap,
                                blacklist_drop_fraction)
  }
  seg_len = segs$end - segs$start
  raw = segs[[value_col]]

  if(centering == "mode"){
    baseline = segmented_data_to_ploidy_estimate(
      segs, method = "mode", autosomes_only = centering_autosomes_only
    )
    centre = baseline$ploidy[match(segs$ID, baseline$sample_id)]
    centre[is.na(centre)] = 2
  }else if(centering == "median"){
    in_baseline = rep(TRUE, nrow(segs))
    if(centering_autosomes_only){
      in_baseline = !grepl("^[XY]$", segs$chrom)
    }
    medians = vapply(split(seq_len(nrow(segs))[in_baseline], segs$ID[in_baseline]),
                     function(ix) .weighted_region_value(raw[ix], seg_len[ix],
                                                         0.5, FALSE),
                     numeric(1))
    centre = medians[segs$ID]
    centre[is.na(centre)] = ifelse(value_scale == "log_ratio", 0, 2)
  }else{
    centre = rep(ifelse(value_scale == "log_ratio", 0, 2), nrow(segs))
  }

  if(value_scale == "log_ratio"){
    # copy number change relative to diploid, as in DLBClass
    segs$.value = 2^(raw - centre + 1) - 2
  }else{
    segs$.value = raw - centre
  }

  # overlap of every segment with every feature on its chromosome
  overlaps = dplyr::inner_join(
    dplyr::mutate(segs, .seg_row = dplyr::row_number()),
    dplyr::select(cnv_features, name, type, region_name, chrom,
                  f_start = start, f_end = end, direction, min_fraction,
                  min_coverage, .background),
    by = "chrom", relationship = "many-to-many"
  ) %>%
    dplyr::mutate(.overlap = pmin(end, f_end) - pmax(start, f_start)) %>%
    dplyr::filter(.overlap > 0)

  overlaps$.event_value = overlaps$.value
  if(any(!is.na(cnv_features$.background))){
    arm_levels = overlaps %>%
      dplyr::filter(type == "arm",
                    region_name %in% cnv_features$.background) %>%
      dplyr::distinct(ID, region_name, .seg_row, .keep_all = TRUE) %>%
      dplyr::group_by(ID, .background = region_name) %>%
      dplyr::summarise(.arm_level = .weighted_region_value(.value, .overlap,
                                                           0.5, FALSE),
                       .groups = "drop")
    overlaps = dplyr::left_join(overlaps, arm_levels,
                                by = c("ID", ".background")) %>%
      dplyr::mutate(
        .event_value = ifelse(is.na(.arm_level), .value, .value - .arm_level)
      ) %>%
      dplyr::select(-.arm_level)
  }

  list(segs = segs,
       overlaps = overlaps,
       cnv_features = cnv_features,
       all_samples = all_samples)
}

#' @title Call arm-level and focal CNV events per sample
#'
#' @description Call copy number events defined in a feature table (by
#' default the curated features from [get_cnv_features]) for each sample in
#' a set of segments. With the DLBClass feature set, the DLBClass blacklist
#' (see `blacklist`) and `output_format = "gsm"`, this produces the copy
#' number portion of the DLBClass GSM input, following what the DLBClass
#' scripts are intended to do (see examples).
#'
#' @details Processing steps:
#' \enumerate{
#'   \item Optionally trim or drop segments overlapping blacklisted regions.
#'     Only blacklist entries overlapping a feature by more than
#'     `blacklist_min_feature_overlap` of that feature's length are used.
#'     A segment entirely inside a blacklisted region is dropped, one
#'     partially overlapping it is trimmed, and one containing it is dropped
#'     only if the region makes up more than `blacklist_drop_fraction` of it.
#'   \item Centre each sample on its baseline (`centering`) and convert to a
#'     copy number change: for `value_scale = "log_ratio"` this is
#'     `2^(log.ratio - median + 1) - 2`, for `"CN"` it is `CN - baseline`.
#'   \item For each feature, take the segments overlapping it (clipped to the
#'     feature). Features with a background arm (see `background`) are
#'     measured relative to the sample's level on that arm (its
#'     length-weighted median).
#'   \item If the overlapping segments cover less than `min_coverage` of the
#'     feature, no event is called. Otherwise the event value is the level
#'     reached by at least `min_fraction` of the covered length, counting
#'     from the most extreme values in the feature's direction (0.5 is the
#'     length-weighted median).
#'   \item The event level is 1 if the value passes `single_thresh` and 2 if
#'     it passes `high_thresh` (above for gains, below the negative for
#'     losses).
#' }
#' Overlaps are inclusive and coverage uses segment lengths clipped to each
#' feature.
#'
#' The thresholds are on the copy number change scale, and their defaults
#' follow `value_scale`: 0.1 and 0.9 for "log_ratio" (as in DLBClass, where
#' impure tumours give small changes), and 0.5 and 1.5 for "CN" (a one-copy
#' change is level 1, two or more copies is level 2). Values you supply are
#' used as given.
#'
#' Only grch37/hg19 is currently supported: `genome_build` must be given (or
#' carried by a `seg_data` object), and segments that run past grch37
#' chromosome ends are rejected as likely hg38.
#'
#' @param seg_data Segmented data with columns ID, chrom, start, end and
#' log.ratio (or CN when `value_scale = "CN"`).
#' @param cnv_features Feature table (see [cnv_features_grch37] for the
#' format). Default is all curated features with include = TRUE, from
#' [get_cnv_features].
#' @param these_samples_metadata Optional metadata whose sample_id column
#' defines the samples returned. Samples without segments get no events.
#' @param genome_build Genome build of `seg_data`. Only "grch37" (or "hg19")
#' is currently supported.
#' @param value_scale "log_ratio" (default) or "CN".
#' @param centering How each sample's genome-wide baseline is set before
#' calling. All are computed internally from `seg_data` (after blacklist
#' trimming), so nothing needs to be run beforehand:
#' \itemize{
#'   \item "median" (default): the length-weighted median of the sample's
#'     values (log.ratio or CN, depending on `value_scale`).
#'   \item "mode": the sample's dominant CN level, from
#'     [segmented_data_to_ploidy_estimate] with `method = "mode"`. Requires
#'     `value_scale = "CN"`. Samples without an estimate use 2.
#'   \item "none": no sample-specific baseline; values are relative to
#'     diploid (log.ratio 0 or CN 2).
#' }
#' @param centering_autosomes_only Exclude chrX/chrY when computing the
#' sample baseline (passed as `autosomes_only` for "mode"). Default FALSE, as
#' in DLBClass.
#' @param single_thresh Minimum copy number change for a level 1 event.
#' Default NULL: 0.1 for `value_scale = "log_ratio"`, 0.5 for "CN".
#' @param high_thresh Minimum copy number change for a level 2 event.
#' Default NULL: 0.9 for `value_scale = "log_ratio"`, 1.5 for "CN".
#' @param background How to choose the arm, if any, that each non-arm
#' feature is measured relative to:
#' \itemize{
#'   \item "auto" (default): the DLBClass rule, applied to the `cnv_features`
#'     table in use. A feature is measured relative to an arm feature in the
#'     same table that contains it and has the same direction (e.g. a focal
#'     18q gain relative to the 18q gain arm feature). Adding an arm feature
#'     to a table therefore makes the focal features on that arm, in the same
#'     direction, relative to it. A non-empty `background_arm` overrides the
#'     rule for that row: an arm-type `region_name` forces that arm, and
#'     "none" forces the genome baseline.
#'   \item "column": only the arms named in `background_arm`.
#'   \item "none": every feature is measured against the genome baseline.
#' }
#' @param blacklist Optional data frame of regions to trim from segments
#' before calling, such as known germline CNVs, with columns chrom, start and
#' end (Chromosome, Start and End are also accepted; chromosomes 23 and 24
#' are read as X and Y). Default NULL (none). No blacklist is bundled:
#' DLBClass used the GISTIC germline CNV file `CNV.hg19.bypos.111213.txt`
#' (from the GISTIC hg19 support files) with one added region, distributed in
#' the DLBCL-Classifier repository as
#' `CNV.hg19.bypos.111213.CR1_event_added.bed`.
#' @param blacklist_min_feature_overlap Default 0.01. See Details.
#' @param blacklist_drop_fraction Default 0.8. See Details.
#' @param encoding "level" (0/1/2, default), "binary" (0/1) or "value" (the
#' event's copy number change, NA when the feature isn't covered).
#' @param output_format "matrix" (default; samples as rows, features named by
#' `name` as columns) or "gsm" (DLBClass layout: a classifier_name column
#' using `dlbclass_name`, then one column per sample).
#' @param gene_labels Set to TRUE to add each feature's representative_gene
#' to its name, e.g. "2p16.1_gain" becomes "2p16.1~REL_gain". Features
#' without a representative gene keep their name. Ignored for
#' `output_format = "gsm"`, which always uses the DLBClass names. Default
#' FALSE.
#' @param gene_separator Separator between the region and the gene when
#' `gene_labels = TRUE`. Default "~".
#'
#' @return A data frame in the requested format.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' segs = GAMBLR.open::get_cn_segments(these = dlbcl_meta)
#' cnv_status = segmented_data_to_cnv_status(segs, genome_build = "grch37")
#'
#' # CN portion of the DLBClass GSM, using the blacklist from the
#' # DLBCL-Classifier repository
#' dlbclass_blacklist = read.delim("CNV.hg19.bypos.111213.CR1_event_added.bed")
#' gsm = segmented_data_to_cnv_status(
#'   segs,
#'   cnv_features = get_cnv_features("grch37", feature_set = "DLBClass"),
#'   genome_build = "grch37",
#'   these_samples_metadata = dlbcl_meta,
#'   blacklist = dlbclass_blacklist,
#'   output_format = "gsm"
#' )
#' }
segmented_data_to_cnv_status = function(seg_data,
                                        cnv_features = NULL,
                                        these_samples_metadata,
                                        genome_build,
                                        value_scale = "log_ratio",
                                        centering = "median",
                                        centering_autosomes_only = FALSE,
                                        single_thresh = NULL,
                                        high_thresh = NULL,
                                        background = "auto",
                                        blacklist = NULL,
                                        blacklist_min_feature_overlap = 0.01,
                                        blacklist_drop_fraction = 0.8,
                                        encoding = "level",
                                        output_format = "matrix",
                                        gene_labels = FALSE,
                                        gene_separator = "~"){
  if(!encoding %in% c("level", "binary", "value")){
    stop("encoding must be \"level\", \"binary\" or \"value\"")
  }
  if(!output_format %in% c("matrix", "gsm")){
    stop("output_format must be \"matrix\" or \"gsm\"")
  }
  if(missing(these_samples_metadata)){
    these_samples_metadata = NULL
  }
  if(missing(genome_build)){
    genome_build = NULL
  }
  prepared = .prepare_cnv_segments(
    seg_data, cnv_features, genome_build, these_samples_metadata,
    value_scale, centering, centering_autosomes_only, blacklist,
    blacklist_min_feature_overlap, blacklist_drop_fraction,
    background
  )
  cnv_features = prepared$cnv_features
  all_samples = prepared$all_samples
  thresholds = .resolve_cnv_thresholds(value_scale, single_thresh, high_thresh)
  single_thresh = thresholds[["single"]]
  high_thresh = thresholds[["high"]]

  calls = prepared$overlaps %>%
    dplyr::group_by(ID, name, direction, f_start, f_end,
                    min_fraction, min_coverage) %>%
    dplyr::summarise(
      .covered = sum(.overlap),
      .event_value = .weighted_region_value(.event_value, .overlap,
                                            min_fraction[1],
                                            direction[1] == "gain"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      .event_value = ifelse(.covered < min_coverage * (f_end - f_start),
                            NA_real_, .event_value),
      .signed = ifelse(direction == "gain", .event_value, -.event_value),
      level = dplyr::case_when(
        is.na(.signed) ~ 0,
        .signed > high_thresh ~ 2,
        .signed > single_thresh ~ 1,
        TRUE ~ 0
      )
    )

  result = matrix(ifelse(encoding == "value", NA_real_, 0),
                  nrow = length(all_samples), ncol = nrow(cnv_features),
                  dimnames = list(all_samples, cnv_features$name))
  cell = cbind(match(calls$ID, all_samples), match(calls$name, cnv_features$name))
  result[cell] = switch(encoding,
                        level = calls$level,
                        binary = as.numeric(calls$level > 0),
                        value = calls$.event_value)

  if(output_format == "gsm"){
    # with background = "auto", extra arm features change how DLBClass focal
    # features are measured, so a mixed table doesn't give DLBClass calls
    in_dlbclass = grepl("(^|,)\\s*DLBClass\\s*(,|$)",
                        if(is.null(cnv_features$feature_set)) ""
                        else ifelse(is.na(cnv_features$feature_set), "",
                                    cnv_features$feature_set))
    if(!all(in_dlbclass)){
      warning(paste0(
        sum(!in_dlbclass), " feature(s) are not in the DLBClass feature set, ",
        "so this output is not a DLBClass GSM (with background = \"auto\" ",
        "they can also change how DLBClass features are called). Use ",
        "get_cnv_features(genome_build, feature_set = \"DLBClass\") for ",
        "DLBClass input."
      ))
    }
    if("dlbclass_name" %in% colnames(cnv_features)){
      classifier_name = ifelse(is.na(cnv_features$dlbclass_name),
                               toupper(cnv_features$name),
                               cnv_features$dlbclass_name)
    }else{
      classifier_name = toupper(cnv_features$name)
    }
    gsm = data.frame(classifier_name = classifier_name,
                     t(result), check.names = FALSE)
    rownames(gsm) = NULL
    return(gsm)
  }
  colnames(result) = .label_cnv_feature_names(cnv_features, gene_labels,
                                              gene_separator)
  as.data.frame(result)
}

#' @title Annotate segments with the CNV events they support
#'
#' @description Return one row per segment, restricted to segments that
#' support at least one event in a feature table, with one column per
#' feature type (e.g. `arm_events`, `focal_events`) listing the supported
#' events separated by ";".
#'
#' @details A segment supports an event when it overlaps the feature and its
#' copy number change (after the same blacklist trimming, centering and arm
#' adjustment as [segmented_data_to_cnv_status]) passes `single_thresh` in
#' the feature's direction. This is per-segment evidence: whether a sample
#' has the event also depends on coverage and the fraction of the feature
#' affected, which [segmented_data_to_cnv_status] applies.
#'
#' Coordinates are returned after blacklist trimming, so they can differ
#' from the input. With `clip_to_features = TRUE` each supporting segment is
#' instead cut to the parts covered by the features it supports: overlapping
#' features (e.g. an arm and a gene on it) are merged into one piece, so
#' rows never overlap, and a segment supporting separate features gives one
#' row per piece, each listing the events it overlaps. This is useful for
#' showing only the driver regions, e.g. via [segmented_data_to_cn_matrix].
#'
#' Only grch37/hg19 is currently supported.
#'
#' @inheritParams segmented_data_to_cnv_status
#' @param keep_all Set to TRUE to also return segments that support no
#' events. Default FALSE. Can't be combined with `clip_to_features`.
#' @param clip_to_features Set to TRUE to trim each supporting segment to
#' the features it supports (see Details). Default FALSE.
#' @param gene_labels Set to TRUE to add each feature's representative_gene
#' to the event names, e.g. "2p16.1~REL_gain". Default FALSE.
#' @param gene_separator Separator between the region and the gene when
#' `gene_labels = TRUE`. Default "~".
#'
#' @return The segments (or, with `clip_to_features = TRUE`, the clipped
#' pieces of them) with added columns cnv_value (the centred copy number
#' change) and `<type>_events` for each feature type.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' segs = GAMBLR.open::get_cn_segments(these = dlbcl_meta)
#' annotated = annotate_cnv(segs, genome_build = "grch37")
#'
#' # heatmap of only the parts of segments that support curated events
#' driver_segs = annotate_cnv(segs, genome_build = "grch37",
#'                            clip_to_features = TRUE)
#' cn_mat = segmented_data_to_cn_matrix(driver_segs,
#'                                      these_samples_metadata = dlbcl_meta,
#'                                      fill_missing_with = "diploid",
#'                                      genome_build = "grch37")
#' }
annotate_cnv = function(seg_data,
                        cnv_features = NULL,
                        these_samples_metadata,
                        genome_build,
                        value_scale = "log_ratio",
                        centering = "median",
                        centering_autosomes_only = FALSE,
                        single_thresh = NULL,
                        background = "auto",
                        blacklist = NULL,
                        blacklist_min_feature_overlap = 0.01,
                        blacklist_drop_fraction = 0.8,
                        keep_all = FALSE,
                        clip_to_features = FALSE,
                        gene_labels = FALSE,
                        gene_separator = "~"){
  if(keep_all && clip_to_features){
    stop("keep_all and clip_to_features can't both be TRUE")
  }
  if(missing(these_samples_metadata)){
    these_samples_metadata = NULL
  }
  if(missing(genome_build)){
    genome_build = NULL
  }
  prepared = .prepare_cnv_segments(
    seg_data, cnv_features, genome_build, these_samples_metadata,
    value_scale, centering, centering_autosomes_only, blacklist,
    blacklist_min_feature_overlap, blacklist_drop_fraction,
    background
  )
  # annotate_cnv has no high_thresh, so don't check against its default
  single_thresh = .resolve_cnv_thresholds(
    value_scale, single_thresh,
    if(is.null(single_thresh)) NULL else single_thresh
  )[["single"]]
  segs = dplyr::mutate(prepared$segs, .seg_row = dplyr::row_number())
  feature_types = unique(prepared$cnv_features$type)
  event_names = stats::setNames(
    .label_cnv_feature_names(prepared$cnv_features, gene_labels, gene_separator),
    prepared$cnv_features$name
  )

  supported = prepared$overlaps %>%
    dplyr::filter(ifelse(direction == "gain",
                         .event_value > single_thresh,
                         .event_value < -single_thresh))

  if(clip_to_features){
    # each supported feature's part of the segment, merged where they
    # overlap so the returned pieces of a segment are disjoint
    parts = supported %>%
      dplyr::transmute(.seg_row, name, type,
                       .part_start = pmax(start, f_start),
                       .part_end = pmin(end, f_end)) %>%
      dplyr::arrange(.seg_row, .part_start) %>%
      dplyr::group_by(.seg_row) %>%
      dplyr::mutate(.piece = cumsum(
        .part_start > dplyr::lag(cummax(.part_end), default = -Inf)
      )) %>%
      dplyr::ungroup()
    pieces = parts %>%
      dplyr::group_by(.seg_row, .piece) %>%
      dplyr::summarise(.start = min(.part_start), .end = max(.part_end),
                       .groups = "drop")
    for(t in feature_types){
      events = parts %>%
        dplyr::filter(type == t) %>%
        dplyr::group_by(.seg_row, .piece) %>%
        dplyr::summarise(.events = paste(event_names[name], collapse = ";"),
                         .groups = "drop")
      pieces[[paste0(t, "_events")]] = events$.events[
        match(paste(pieces$.seg_row, pieces$.piece),
              paste(events$.seg_row, events$.piece))
      ]
    }
    segs = dplyr::inner_join(segs, pieces, by = ".seg_row") %>%
      dplyr::mutate(start = .start, end = .end) %>%
      dplyr::arrange(.seg_row, start) %>%
      dplyr::select(-.start, -.end, -.piece)
  }else{
    for(t in feature_types){
      events = supported %>%
        dplyr::filter(type == t) %>%
        dplyr::group_by(.seg_row) %>%
        dplyr::summarise(.events = paste(event_names[name], collapse = ";"),
                         .groups = "drop")
      segs[[paste0(t, "_events")]] = events$.events[match(segs$.seg_row,
                                                          events$.seg_row)]
    }
  }
  event_cols = paste0(feature_types, "_events")
  if(!keep_all){
    any_event = rowSums(!is.na(segs[, event_cols, drop = FALSE])) > 0
    segs = segs[any_event, , drop = FALSE]
  }
  segs = dplyr::rename(segs, cnv_value = .value) %>%
    dplyr::select(-.seg_row)
  rownames(segs) = NULL
  segs
}
