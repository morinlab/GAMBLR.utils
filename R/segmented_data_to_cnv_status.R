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
      "not supported yet. Affected samples:\n  ", .list_ids(affected), "\n",
      .remove_failing_hint(check = "build_ok")
    ), call. = FALSE)
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
  # optional focality limit (see .drop_broad_cnv_calls)
  if(!"max_event_length" %in% colnames(cnv_features)){
    cnv_features$max_event_length = NA_real_
  }
  cnv_features$max_event_length = as.numeric(cnv_features$max_event_length)
  as.data.frame(cnv_features)
}

# Lower calls that are part of a longer change, for features with a
# max_event_length (e.g. the LymphGen gene features, which count only
# changes of 30 Mb or less). The extent of a call is found by starting from
# the segments overlapping the feature that are at exactly the call's level
# and walking outward along the sample's segments on that chromosome for as
# long as they stay at that level. Uncovered stretches between segments
# (liftOver gaps, centromeres) are stepped over without adding to the
# extent, so a segment split into pieces by liftOver has the same extent as
# the unsplit segment. Segments at a different level stop the walk, as
# LymphGen judges each segment by its own length: a short single-copy loss
# next to a homozygous deletion still counts. If the extent is longer than
# max_event_length, the call drops to the next level down and is checked
# again with the segments at that level; with none overlapping the feature,
# the call drops again (a long homozygous deletion counts as nothing, while
# one inside an arm-level loss keeps level 2).
# With gene_rule = "lymphgen", features with a max_event_length are instead
# called as LymphGen's training data are (see segmented_data_to_cnv_status).
.drop_broad_cnv_calls = function(calls, segs, cnv_features, single_thresh,
                                 high_thresh, event_extent = "run", levels = NULL,
                                 gene_rule = "fraction"){
  limited = cnv_features[!is.na(cnv_features$max_event_length), , drop = FALSE]
  # with LymphGen's rule a feature can be called from any overlap, so every
  # assessed feature is checked, not just those already called
  check = if(gene_rule == "lymphgen"){
    which(calls$name %in% limited$name & !is.na(calls$.event_value))
  } else which(calls$name %in% limited$name & calls$level > 0)
  if(length(check) == 0){
    if(gene_rule == "lymphgen") calls$.rule_status = NA_character_
    return(calls)
  }
  if(!".centre" %in% colnames(segs)) segs$.centre = NA_real_
  if(!".rule_len" %in% colnames(segs)) segs$.rule_len = segs$end - segs$start
  segs_by_sample = split(segs[, c("chrom", "start", "end", ".value", ".centre", ".rule_len")], segs$ID)
  feature = limited[match(calls$name[check], limited$name), , drop = FALSE]
  extent = function(sg, f_start, f_end, in_state){
    hit = which(in_state & sg$end > f_start & sg$start < f_end)
    if(length(hit) == 0){
      return(Inf)
    }
    first = min(hit)
    while(first > 1 && in_state[first - 1]) first = first - 1
    last = max(hit)
    while(last < nrow(sg) && in_state[last + 1]) last = last + 1
    sg$end[last] - sg$start[first]
  }
  sample_segs = function(i, k){
    sg = segs_by_sample[[calls$ID[i]]]
    sg = sg[sg$chrom == feature$chrom[k], , drop = FALSE]
    sg[order(sg$start), , drop = FALSE]
  }
  if(gene_rule == "lymphgen"){
    # one state per gene, decided by one segment no longer than
    # max_event_length: level 2 (HOMDEL, AMP) beats level 1 in either
    # direction, then the largest overlap with the feature; the feature gets
    # that segment's level if it is in the feature's direction. The value is
    # the deciding segment's change, and .rule_status says why the feature
    # is or isn't called
    decided = lapply(seq_along(check), function(k){
      i = check[k]
      sg = sample_segs(i, k)
      overlap = pmin(sg$end, feature$end[k]) - pmax(sg$start, feature$start[k])
      gain_level = .seg_levels(sg, "gain", high_thresh, single_thresh, levels)
      loss_level = .seg_levels(sg, "loss", high_thresh, single_thresh, levels)
      seg_level = pmax(gain_level, loss_level)
      short = sg$.rule_len <= feature$max_event_length[k]
      cand = which(overlap > 0 & short & seg_level > 0)
      if(length(cand) == 0){
        own = if(calls$direction[i] == "gain") gain_level else loss_level
        too_long = any(overlap > 0 & !short & own > 0)
        return(list(level = 0, value = calls$.event_value[i],
                    status = if(too_long) "event too long" else "below threshold"))
      }
      pick = cand[order(-seg_level[cand], -overlap[cand])[1]]
      pick_direction = if(gain_level[pick] > 0) "gain" else "loss"
      if(pick_direction == calls$direction[i]){
        list(level = seg_level[pick], value = sg$.value[pick], status = "called")
      }else{
        list(level = 0, value = sg$.value[pick], status = "outweighed by opposite change")
      }
    })
    calls$.rule_status = NA_character_
    calls$level[check] = vapply(decided, `[[`, 0, "level")
    calls$.event_value[check] = vapply(decided, `[[`, 0, "value")
    calls$.rule_status[check] = vapply(decided, `[[`, "", "status")
    return(calls)
  }
  calls$level[check] = vapply(seq_along(check), function(k){
    i = check[k]
    sg = sample_segs(i, k)
    seg_level = .seg_levels(sg, calls$direction[i], high_thresh, single_thresh, levels)
    level = calls$level[i]
    if(event_extent == "segment"){
      # each overlapping segment on its own raw length (LGenIC, LymphGen)
      hit = sg$end > feature$start[k] & sg$start < feature$end[k] &
        sg$.rule_len <= feature$max_event_length[k]
      return(min(level, max(c(0, seg_level[hit]))))
    }
    while(level > 0){
      if(extent(sg, feature$start[k], feature$end[k], seg_level == level) <=
         feature$max_event_length[k]){
        break
      }
      level = level - 1
    }
    level
  }, numeric(1))
  calls
}

# Each segment's level (0, 1, 2) in one direction
.seg_levels = function(sg, direction, high_thresh, single_thresh, levels){
  if(is.null(levels)){
    signed = if(direction == "gain") sg$.value else -sg$.value
    return(ifelse(signed > high_thresh, 2, ifelse(signed > single_thresh, 1, 0)))
  }
  .cnv_levels(sg$.value, sg$.centre, rep(direction, nrow(sg)), levels)
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
# Baseline ("centre") of each segment's sample, on the scale of its values:
# the sample's length-weighted median ("median"), its dominant CN level
# ("mode", CN only) or the diploid level ("none"). segs needs ID and chrom;
# values and lengths are aligned to its rows. Shared by
# segmented_data_to_cnv_status and calculate_pga so both centre the same way.
.cnv_sample_centre = function(segs, values, lengths, value_scale, centering,
                              autosomes_only){
  if(!centering %in% c("median", "mode", "none")){
    stop("centering must be \"median\", \"mode\" or \"none\"")
  }
  if(centering == "mode" && value_scale != "CN"){
    stop("centering = \"mode\" requires value_scale = \"CN\"")
  }
  diploid = ifelse(value_scale == "log_ratio", 0, 2)
  if(nrow(segs) == 0){
    return(numeric(0))
  }
  # placeholder rows (dummy_segment == 1) never count towards the baseline
  real = .real_segments(segs)
  if(centering == "mode"){
    baseline = segmented_data_to_ploidy_estimate(
      data.frame(ID = segs$ID, chrom = segs$chrom,
                 start = segs$start, end = segs$end, CN = values)[real, , drop = FALSE],
      method = "mode", autosomes_only = autosomes_only
    )
    centre = baseline$ploidy[match(segs$ID, baseline$sample_id)]
  }else if(centering == "median"){
    in_baseline = real
    if(autosomes_only){
      in_baseline = in_baseline & .is_autosome(segs$chrom)
    }
    rows = seq_len(nrow(segs))[in_baseline]
    medians = vapply(split(rows, segs$ID[in_baseline]),
                     function(ix) .weighted_region_value(values[ix], lengths[ix],
                                                         0.5, FALSE),
                     numeric(1))
    centre = unname(medians[as.character(segs$ID)])
  }else{
    centre = rep(diploid, nrow(segs))
  }
  centre[is.na(centre)] = diploid
  centre
}

.prepare_cnv_segments = function(seg_data, cnv_features, genome_build,
                                 these_samples_metadata, value_scale,
                                 centering, centering_autosomes_only,
                                 blacklist, blacklist_min_feature_overlap,
                                 blacklist_drop_fraction, background,
                                 hygiene = NULL){
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

  centre = .cnv_sample_centre(segs, raw, seg_len, value_scale, centering,
                              centering_autosomes_only)

  if(value_scale == "log_ratio"){
    # copy number change relative to diploid, as in DLBClass
    segs$.value = 2^(raw - centre + 1) - 2
  }else{
    segs$.value = raw - centre
  }
  segs$.centre = centre
  segs$.rule_len = segs$end - segs$start
  # optional segment hygiene: isolated short gain spikes are ignored, and
  # small gaps between segments at the same level are bridged
  if(!is.null(hygiene)){
    if(hygiene$min_gain_spike > 0){
      segs = .drop_gain_spikes(segs, hygiene$min_gain_spike, hygiene$single_thresh)
    }
    if(hygiene$bridge_gaps > 0){
      segs = .bridge_small_gaps(segs, hygiene$bridge_gaps, hygiene)
    }
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
#'   \item For features with a `max_event_length` (e.g. the LymphGen gene
#'     features, which count only changes of 30 Mb or less), a call is kept
#'     only if the change is no longer than that. Its extent is measured
#'     from the segments overlapping the feature that are at exactly the
#'     call's level, outward along adjacent segments for as long as they stay
#'     at that level; uncovered stretches between segments (e.g. liftOver
#'     gaps or the centromere) are stepped over without adding to the extent,
#'     so segments split by liftOver give the same result as the unsplit
#'     segment, while a segment at another level ends the change (as in
#'     LymphGen, which judges each segment by its own length). With
#'     `event_extent = "segment"`, each overlapping segment is instead judged
#'     on its own raw length, as LGenIC and LymphGen's training data do: the
#'     call is the highest level among the overlapping segments no longer than
#'     `max_event_length`, so a broad change split into shorter segments
#'     counts in each of them. A level 2 call
#'     that is too long is checked again as level 1 using the level 1
#'     segments overlapping the feature, if any; a deeper focal change inside
#'     a broad one (e.g. a homozygous deletion in an arm-level loss) keeps its
#'     level. With
#'     `encoding = "value"`, the measured change is still reported.
#'   \item With `gene_rule = "lymphgen"`, features with a `max_event_length`
#'     are called as LymphGen's training data are (recovered exactly from the
#'     segment labels behind them): one segment decides each feature, and
#'     only segments no longer than `max_event_length` on their own raw
#'     length that overlap the feature at all count (`min_fraction` and
#'     `event_extent` don't apply). A level 2 segment (HOMDEL or AMP) beats
#'     a level 1 segment in either direction; among segments of the same
#'     level, the one overlapping the feature most decides. The feature gets
#'     that segment's level if it is in the feature's direction and 0
#'     otherwise, so a gene never has both a gain and a loss: a short loss
#'     inside a larger gain over the gene is not a loss. `min_coverage`
#'     still applies. With `encoding = "value"` the value is the deciding
#'     segment's change (positive for a gain, even for a loss feature it
#'     outweighs), and the "long" output gets a `rule_status` column:
#'     "called", "outweighed by opposite change", "event too long" (only
#'     longer changes in the feature's direction overlap it) or "below
#'     threshold" (NA for features without a `max_event_length`).
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
#' @param centering_autosomes_only Use autosomes (1-22) only when computing
#' the sample baseline, as every other GAMBLR baseline does. Default TRUE.
#' `output_format = "gsm"` always uses FALSE (all chromosomes), as DLBClass
#' does, so its output matches DLBClass's.
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
#' `name` as columns), "long" (one row per sample and feature, with columns
#' sample_id, name, chrom, start and end, then the call in a column named
#' after `encoding`: level, binary or value; useful e.g. for accounting for
#' features that share a chromosome when modelling) or "gsm" (DLBClass
#' layout: a classifier_name column using `dlbclass_name`, then one column
#' per sample).
#' @param gene_labels Set to TRUE to add each feature's representative_gene
#' to its name, e.g. "2p16.1_gain" becomes "2p16.1~REL_gain". Features
#' without a representative gene keep their name. Ignored for
#' `output_format = "gsm"`, which always uses the DLBClass names. With
#' `output_format = "long"`, name keeps the feature name and the labelled
#' name goes in an added label column. Default FALSE.
#' @param gene_separator Separator between the region and the gene when
#' `gene_labels = TRUE`. Default "~".
#' @param homdel_max_cn On the CN scale, the copy number at or below which a
#' loss is level 2 (a homozygous deletion), whatever the sample's baseline;
#' level 2 gains then need a change of `high_thresh` scaled by baseline / 2
#' (e.g. 3 copies on a tetraploid baseline). Default 0.5. NULL (or NA) judges
#' level 2 as a change of `high_thresh` from the baseline, as on the log ratio
#' scale (where it always is), and as LymphGen's training data do: its
#' homozygous deletions and amplifications are changes of at least 1.5 copies
#' from the tumour's own baseline (e.g. 4 to 2 copies in a tetraploid tumour).
#' @param bridge_gaps Gaps between neighbouring segments up to this length
#' (bp) are bridged when both segments are at the same level, at the less
#' extreme of their values, so a feature in the gap (e.g. a gene without
#' capture targets) can be assessed. Default 1e6; 0 to leave gaps.
#' @param min_gain_spike Isolated gain segments shorter than this (bp) are
#' ignored: those more than `single_thresh` above both neighbours, whose
#' values agree, e.g. single-bin spikes. Losses are always kept. Default 1e5;
#' 0 to keep every segment.
#' @param event_extent How the length of a change is measured for features
#' with a `max_event_length` (see Details): "run" (default; the change's run
#' of neighbouring segments at the call's level) or "segment" (each segment
#' on its own raw length, as LGenIC and LymphGen's training data).
#' @param gene_rule How features with a `max_event_length` are called:
#' "fraction" (default; the level reached by `min_fraction` of the feature,
#' then `event_extent`) or "lymphgen" (one segment per feature, as
#' LymphGen's training data; see Details).
#'
#' @return A data frame in the requested format. The "long" format is sorted
#' by sample, then by feature in the order of `cnv_features`.
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
#'
#' # one row per sample and feature, with the feature coordinates
#' cnv_long = segmented_data_to_cnv_status(segs, genome_build = "grch37",
#'                                         output_format = "long")
#' }
segmented_data_to_cnv_status = function(seg_data,
                                        cnv_features = NULL,
                                        these_samples_metadata,
                                        genome_build,
                                        value_scale = "log_ratio",
                                        centering = "median",
                                        centering_autosomes_only = TRUE,
                                        single_thresh = NULL,
                                        high_thresh = NULL,
                                        background = "auto",
                                        blacklist = NULL,
                                        blacklist_min_feature_overlap = 0.01,
                                        blacklist_drop_fraction = 0.8,
                                        encoding = "level",
                                        output_format = "matrix",
                                        gene_labels = FALSE,
                                        gene_separator = "~",
                                        event_extent = "run",
                                        gene_rule = "fraction",
                                        homdel_max_cn = 0.5,
                                        bridge_gaps = 1e6,
                                        min_gain_spike = 1e5){
  event_extent = match.arg(event_extent, c("run", "segment"))
  gene_rule = match.arg(gene_rule, c("fraction", "lymphgen"))
  if(!encoding %in% c("level", "binary", "value")){
    stop("encoding must be \"level\", \"binary\" or \"value\"")
  }
  if(!output_format %in% c("matrix", "long", "gsm")){
    stop("output_format must be \"matrix\", \"long\" or \"gsm\"")
  }
  # DLBClass centres on all chromosomes; the GSM follows it
  if(output_format == "gsm" && centering_autosomes_only){
    if(!missing(centering_autosomes_only)){
      message("output_format = \"gsm\" centres on all chromosomes, as DLBClass does ",
              "(centering_autosomes_only = TRUE is ignored)")
    }
    centering_autosomes_only = FALSE
  }
  if(missing(these_samples_metadata)){
    these_samples_metadata = NULL
  }
  if(missing(genome_build)){
    genome_build = NULL
  }
  thresholds = .resolve_cnv_thresholds(value_scale, single_thresh, high_thresh)
  single_thresh = thresholds[["single"]]
  high_thresh = thresholds[["high"]]
  # level 2 on the CN scale is judged against the sample's ploidy (see
  # .cnv_levels); NULL keeps it relative to the baseline
  levels = list(value_scale = value_scale, single = single_thresh, high = high_thresh,
                homdel_max_cn = if(value_scale == "CN" && !is.null(homdel_max_cn) && !is.na(homdel_max_cn)) homdel_max_cn else NULL)
  prepared = .prepare_cnv_segments(
    seg_data, cnv_features, genome_build, these_samples_metadata,
    value_scale, centering, centering_autosomes_only, blacklist,
    blacklist_min_feature_overlap, blacklist_drop_fraction,
    background,
    hygiene = c(levels, list(single_thresh = single_thresh, bridge_gaps = bridge_gaps,
                             min_gain_spike = min_gain_spike))
  )
  cnv_features = prepared$cnv_features
  all_samples = prepared$all_samples

  centres = prepared$segs$.centre[match(unique(prepared$segs$ID), prepared$segs$ID)]
  names(centres) = unique(prepared$segs$ID)
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
      .signed = ifelse(direction == "gain", .event_value, -.event_value)
    )
  calls$level = .cnv_levels(calls$.event_value, centres[calls$ID], calls$direction, levels)
  calls = .drop_broad_cnv_calls(calls, prepared$segs, cnv_features,
                                single_thresh, high_thresh, event_extent, levels, gene_rule)

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
  feature_labels = .label_cnv_feature_names(cnv_features, gene_labels,
                                            gene_separator)
  if(output_format == "long"){
    n_samples = length(all_samples)
    n_features = nrow(cnv_features)
    feature_ix = rep(seq_len(n_features), times = n_samples)
    long = data.frame(sample_id = rep(all_samples, each = n_features),
                      name = cnv_features$name[feature_ix],
                      stringsAsFactors = FALSE)
    if(gene_labels){
      long$label = feature_labels[feature_ix]
    }
    long$chrom = cnv_features$chrom[feature_ix]
    long$start = cnv_features$start[feature_ix]
    long$end = cnv_features$end[feature_ix]
    # rows are sample-major, so read the matrix row by row
    long[[encoding]] = as.vector(t(result))
    if(gene_rule == "lymphgen"){
      # why each feature called by LymphGen's gene rule is or isn't called
      status = matrix(NA_character_, nrow = length(all_samples), ncol = n_features)
      has = !is.na(calls$.rule_status)
      status[cell[has, , drop = FALSE]] = calls$.rule_status[has]
      long$rule_status = as.vector(t(status))
    }
    return(long)
  }
  colnames(result) = feature_labels
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
                        centering_autosomes_only = TRUE,
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

# Event level (0, 1, 2) of copy number changes in a feature's direction.
# Level 1 needs a change beyond `single`. Level 2 needs a change beyond
# `high` on the log ratio scale (relative to a diploid baseline); on the CN
# scale it is judged against the sample's ploidy: a loss is level 2 only when
# the copy number itself is at most homdel_max_cn (a homozygous deletion,
# whatever the baseline), and a gain only when the change reaches `high`
# scaled by baseline / 2 (so 3 copies on a tetraploid baseline). With
# homdel_max_cn NULL, level 2 is relative to the baseline as on the log ratio
# scale.
.cnv_levels = function(change, centre, direction, levels){
  signed = ifelse(direction == "gain", change, -change)
  level = ifelse(is.na(signed), 0, ifelse(signed > levels$single, 1, 0))
  if(is.null(levels$homdel_max_cn)){
    level[!is.na(signed) & signed > levels$high] = 2
    return(level)
  }
  centre = ifelse(is.na(centre), 2, centre)
  absolute = centre + change
  deep_loss = direction == "loss" & !is.na(absolute) & absolute <= levels$homdel_max_cn & level > 0
  big_gain = direction == "gain" & !is.na(signed) & signed > levels$high * pmax(centre, 2) / 2
  level[deep_loss | big_gain] = 2
  level
}

# Segments ignored as artefacts: isolated short gains (shorter than
# max_len, more than `single` above both neighbours, whose own values agree
# within `single`), e.g. single-bin spikes. Losses are kept: a short, deep,
# isolated loss is what a focal homozygous deletion looks like.
.drop_gain_spikes = function(segs, max_len, single){
  segs = segs[order(segs$ID, segs$chrom, segs$start), , drop = FALSE]
  n = nrow(segs)
  if(n < 3) return(segs)
  same_prev = c(FALSE, segs$ID[-1] == segs$ID[-n] & segs$chrom[-1] == segs$chrom[-n])
  same_next = c(same_prev[-1], FALSE)
  v = segs$.value
  prev_v = c(NA, v[-n]); next_v = c(v[-1], NA)
  spike = same_prev & same_next & (segs$end - segs$start) < max_len &
    v - pmax(prev_v, next_v) > single & abs(prev_v - next_v) < single
  spike[is.na(spike)] = FALSE
  segs[!spike, , drop = FALSE]
}

# Small gaps (up to max_gap) between neighbouring segments at the same level
# are bridged with a segment at the less extreme of the two values, so a
# feature in the gap (e.g. a gene without capture targets) can be assessed.
# Gaps between segments at different levels are left: the breakpoint could
# be anywhere in them. A bridge carries the longer flank's length for the
# max_event_length rule.
.bridge_small_gaps = function(segs, max_gap, levels){
  segs = segs[order(segs$ID, segs$chrom, segs$start), , drop = FALSE]
  n = nrow(segs)
  if(n < 2) return(segs)
  i = which(segs$ID[-1] == segs$ID[-n] & segs$chrom[-1] == segs$chrom[-n])
  gap = segs$start[i + 1] - segs$end[i]
  i = i[gap > 0 & gap <= max_gap]
  if(length(i) == 0) return(segs)
  lvl = function(rows, direction) .cnv_levels(segs$.value[rows], segs$.centre[rows], rep(direction, length(rows)), levels)
  same = lvl(i, "gain") == lvl(i + 1, "gain") & lvl(i, "loss") == lvl(i + 1, "loss")
  i = i[same]
  if(length(i) == 0) return(segs)
  milder = ifelse(abs(segs$.value[i]) <= abs(segs$.value[i + 1]), i, i + 1)
  bridge = segs[milder, , drop = FALSE]
  bridge$start = segs$end[i]
  bridge$end = segs$start[i + 1]
  bridge$.rule_len = pmax(segs$end[i] - segs$start[i], segs$end[i + 1] - segs$start[i + 1])
  out = rbind(segs, bridge)
  out[order(out$ID, out$chrom, out$start), , drop = FALSE]
}
