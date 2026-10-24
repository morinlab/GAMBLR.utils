#' Curated CNV features (grch37)
#'
#' Curated arm-level and focal copy number features, in the format used by
#' [segmented_data_to_cnv_status] and [annotate_cnv]. Built by
#' `data-raw/cnv_features.R` from the master table
#' `inst/extdata/cnv_features.grch37.tsv`, which is meant to be extended over
#' time. Use [get_cnv_features] rather than this object directly.
#'
#' The rows tagged "DLBClass" were seeded from the DLBClass arm-level and
#' focal (GISTIC) definitions and are kept fixed, since the DLBClass model
#' was trained on them. They are in DLBClass GSM order. The rows tagged
#' "LymphGen" are the LymphGen copy number features (Wright et al. 2020).
#' Rows with an empty `feature_set` are hand-curated.
#'
#' A feature in more than one set has one row per set (same `name`), so that
#' each named set can be returned complete with its own settings. Where rows
#' are combined across sets ([get_cnv_features] for the core set or "all"),
#' one row is kept per name, by priority: hand-curated first, then LymphGen,
#' then DLBClass. A curated row for a feature therefore overrides a tool's
#' settings for it.
#'
#' grch37/hg19 only. An hg38 version still needs to be created properly
#' (lifted-over focal features and hg38 GISTIC arm boundaries).
#'
#' @format A data frame with one row per feature:
#' \describe{
#'   \item{name}{Event name, e.g. "17p_loss" or "2p16.1_gain"; unique within
#'     each feature set and in the default set.}
#'   \item{type}{Feature type, e.g. "arm" or "focal". [annotate_cnv] returns
#'     one `<type>_events` column per type.}
#'   \item{region_name}{Arm (e.g. "17p") or cytoband (e.g. "2p16.1").}
#'   \item{chrom, start, end}{Coordinates (grch37, no chr prefix). Arm
#'     boundaries (in every set) are at the hg19 cytoband centromere, as in
#'     the DLBClass GISTIC arms, with q arms ending at the chromosome end.}
#'   \item{direction}{"gain" or "loss".}
#'   \item{min_fraction}{Fraction of the covered feature length that must pass
#'     the threshold (0.5 is the length-weighted median).}
#'   \item{min_coverage}{Minimum fraction of the feature that segments must
#'     cover for an event to be called.}
#'   \item{background_arm}{Usually empty. Which arm a feature is measured
#'     relative to is normally decided by the rule described under
#'     `background` in [segmented_data_to_cnv_status] (an arm feature in the
#'     same table that contains it, with the same direction). A value here
#'     records an exception for that feature: an arm-type region_name to
#'     always measure it against, or "none" to always use the genome
#'     baseline. Must be empty for DLBClass rows, which the build script
#'     enforces.}
#'   \item{representative_gene}{A single gene representing the region, used
#'     to label features (e.g. "2p16.1~REL_gain") when `gene_labels = TRUE`.
#'     Left empty where no gene is a convincing target. A few lie just
#'     outside their peak, which `summary` notes.}
#'   \item{genes}{Optional comma-separated additional genes.}
#'   \item{feature_set}{Comma-separated named sets the feature belongs to,
#'     e.g. "DLBClass" or "LymphGen"; empty for hand-curated features.}
#'   \item{dlbclass_name}{Feature name in the DLBClass GSM, e.g. "17P.DEL".}
#'   \item{include}{FALSE retires a row (e.g. a feature redundant with
#'     another) without deleting it: it is left out of the core set and of
#'     `feature_set = "all"`, but still returned with its own named set, so
#'     published sets stay complete.}
#'   \item{core}{TRUE for the rows in the default (core) set returned by
#'     [get_cnv_features] without a `feature_set`: a curated starting point
#'     rather than every feature. Judged per row: marking a feature core on
#'     one set's row uses that row's settings. Initially the hand-curated
#'     rows and the (non-retired) DLBClass rows.}
#'   \item{citation, summary}{Source and description.}
#'   \item{lymphgen_level}{For LymphGen rows, the magnitude of copy number
#'     change (in copies) an event must reach to count for LymphGen: 1 or 2.
#'     Recorded for future use; not yet used when calling features.}
#' }
"cnv_features_grch37"

#' @title Get curated CNV feature definitions
#'
#' @description Return CNV feature definitions for a genome build: the core
#' set (default), every non-redundant feature ("all"), or a named set such
#' as "DLBClass", "LymphGen" or "curated". The
#' result can be used as a template for a custom `cnv_features` table in
#' [segmented_data_to_cnv_status] and [annotate_cnv].
#'
#' @param genome_build Genome build. Only "grch37" (or "hg19") is currently
#' supported.
#' @param feature_set Optional:
#' \itemize{
#'   \item missing (default): the core set, the rows with `core = TRUE`, one
#'     per feature name (see `priority`): a curated starting point;
#'   \item "all": every feature that isn't retired (`include = TRUE`), one
#'     row per feature name: the core set plus everything else, without
#'     redundant rows;
#'   \item a named set ("DLBClass", "LymphGen", or "curated" for the
#'     hand-curated rows): that set in full, regardless of `include` and
#'     `core`, so published sets stay complete.
#' }
#' @param priority For the core set and "all", the order in which a
#' feature's rows from different sets are preferred: "curated" (rows with
#' an empty feature_set) and set names. Default c("curated", "LymphGen",
#' "DLBClass"): a hand-curated row overrides a tool's settings, and LymphGen
#' settings are used over DLBClass ones. Sets not listed come last.
#'
#' @return A data frame of CNV features (see [cnv_features_grch37]).
#'
#' @export
#'
#' @examples
#' cnv_features = get_cnv_features("grch37")
#' dlbclass_features = get_cnv_features("grch37", feature_set = "DLBClass")
#' lymphgen_features = get_cnv_features("grch37", feature_set = "LymphGen")
#' every_feature = get_cnv_features("grch37", feature_set = "all")
get_cnv_features = function(genome_build, feature_set,
                            priority = c("curated", "LymphGen", "DLBClass")){
  if(missing(genome_build)){
    stop("genome_build must be specified (only \"grch37\" is currently supported)")
  }
  if(!genome_build %in% c("grch37", "hg19")){
    # TODO(hg38): needs a properly built hg38 table, see data-raw/cnv_features.R
    stop(paste0("CNV features are currently only available for grch37/hg19 ",
                "(got \"", genome_build, "\")"))
  }
  features = GAMBLR.utils::cnv_features_grch37
  if(missing(feature_set)){
    # the core set: core (and not retired) rows, one per name
    return(.collapse_cnv_features(features[features$include & features$core, , drop = FALSE], priority))
  }
  if(identical(feature_set, "all")){
    # every feature that isn't retired, one row per name
    return(.collapse_cnv_features(features[features$include, , drop = FALSE], priority))
  }
  if(identical(feature_set, "curated")){
    return(features[is.na(features$feature_set), , drop = FALSE])
  }
  in_set = vapply(strsplit(ifelse(is.na(features$feature_set), "",
                                  features$feature_set), ","),
                  function(sets) feature_set %in% trimws(sets),
                  logical(1))
  if(!any(in_set)){
    stop(paste0("no features in feature_set \"", feature_set, "\""))
  }
  features[in_set, , drop = FALSE]
}

# One row per feature name: where a feature has rows in several sets, the
# row whose set comes first in priority ("curated" = empty feature_set; a
# row in several sets ranks by its best one; unlisted sets last), keeping
# the table's row order
.collapse_cnv_features = function(features, priority){
  sets = strsplit(ifelse(is.na(features$feature_set), "curated", features$feature_set), ",")
  rank = vapply(sets, function(s){
    r = match(trimws(s), priority)
    if(all(is.na(r))) length(priority) + 1 else min(r, na.rm = TRUE)
  }, numeric(1))
  best = stats::ave(rank, features$name, FUN = min)
  keep = rank == best & !duplicated(data.frame(features$name, rank))
  features[keep, , drop = FALSE]
}
