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
#' was trained on them. They are in DLBClass GSM order.
#'
#' grch37/hg19 only. An hg38 version still needs to be created properly
#' (lifted-over focal features and hg38 GISTIC arm boundaries).
#'
#' @format A data frame with one row per feature:
#' \describe{
#'   \item{name}{Unique event name, e.g. "17p_loss" or "2p16.1_gain".}
#'   \item{type}{Feature type, e.g. "arm" or "focal". [annotate_cnv] returns
#'     one `<type>_events` column per type.}
#'   \item{region_name}{Arm (e.g. "17p") or cytoband (e.g. "2p16.1").}
#'   \item{chrom, start, end}{Coordinates (grch37, no chr prefix). The
#'     DLBClass arm boundaries are the GISTIC ones and include the centromere.}
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
#'     e.g. "DLBClass".}
#'   \item{dlbclass_name}{Feature name in the DLBClass GSM, e.g. "17P.DEL".}
#'   \item{include}{FALSE retires a feature from the default set without
#'     deleting it.}
#'   \item{citation, summary}{Source and description.}
#' }
"cnv_features_grch37"

#' @title Get curated CNV feature definitions
#'
#' @description Return the curated CNV feature table for a genome build,
#' either the current default set or a named set such as "DLBClass". The
#' result can be used as a template for a custom `cnv_features` table in
#' [segmented_data_to_cnv_status] and [annotate_cnv].
#'
#' @param genome_build Genome build. Only "grch37" (or "hg19") is currently
#' supported.
#' @param feature_set Optional name of a feature set (e.g. "DLBClass"). Named
#' sets are returned in full regardless of `include`, so they stay stable as
#' the curated table changes. When missing, all features with
#' `include = TRUE` are returned.
#'
#' @return A data frame of CNV features (see [cnv_features_grch37]).
#'
#' @export
#'
#' @examples
#' cnv_features = get_cnv_features("grch37")
#' dlbclass_features = get_cnv_features("grch37", feature_set = "DLBClass")
get_cnv_features = function(genome_build, feature_set){
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
    return(features[features$include, , drop = FALSE])
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
