#!/usr/bin/env Rscript
# Build cnv_features_grch37 from the curated master table
# inst/extdata/cnv_features.grch37.tsv. Edit the TSV (add rows, or set
# include = FALSE to retire one) and rerun this script.
#
# core = TRUE marks the rows in the default (core) set; include = FALSE
# retires a row (left out of the core set and of "all", kept in its named
# set).
#
# A feature in more than one set gets one row per set, with the same name
# and only feature_set (and set-specific settings) differing. Hand-curated
# rows (empty feature_set) take priority over tool-specific ones when
# get_cnv_features() is called without a set. Arm boundaries in every set
# are at the hg19 cytoband centromere, as in the DLBClass GISTIC arms (not
# the GAMBLR.data arm table, which leaves out the centromeric gaps).
#
# The rows tagged "DLBClass" in feature_set were seeded from the DLBClass
# arm and focal definitions (DLBCL_broad_significance.18Aug2024 and
# DLBCL_focal_peaks.18Aug2024 in getzlab/DLBCL-Classifier). Keep those rows
# unchanged: the DLBClass model was trained on exactly these features.
#
# TODO(hg38): grch37 only. An hg38 table needs the focal features lifted
# over and hg38 GISTIC arm boundaries (not a naive liftover of the arm
# coordinates).

library(usethis)
library(readr)
library(dplyr)

cnv_features_grch37 = read_tsv(
  "inst/extdata/cnv_features.grch37.tsv",
  col_types = cols(chrom = "c", start = "i", end = "i",
                   include = "l", core = "l", lymphgen_level = "i", max_event_length = "d",
                   .default = "c"),
  na = c("", "NA")
) %>%
  mutate(min_fraction = as.numeric(min_fraction),
         min_coverage = as.numeric(min_coverage)) %>%
  as.data.frame()

in_dlbclass = grepl("(^|,)\\s*DLBClass\\s*(,|$)",
                    ifelse(is.na(cnv_features_grch37$feature_set), "",
                           cnv_features_grch37$feature_set))

# feature names: a feature in several sets has one row per set (same name),
# so names need only be unique within each set (and among the hand-curated
# rows, which have an empty feature_set); get_cnv_features() collapses them
# to one row per name when no set is asked for
set_members = strsplit(ifelse(is.na(cnv_features_grch37$feature_set), "curated",
                              cnv_features_grch37$feature_set), ",")
set_names = unique(trimws(unlist(set_members)))
duplicated_within = unlist(lapply(set_names, function(s){
  in_s = vapply(set_members, function(m) s %in% trimws(m), logical(1))
  n = cnv_features_grch37$name[in_s]
  if(any(duplicated(n))) paste0(s, ": ", paste(unique(n[duplicated(n)]), collapse = ", "))
}))
if(length(duplicated_within) > 0){
  stop("duplicated feature names within a set:\n  ", paste(duplicated_within, collapse = "\n  "))
}
in_lymphgen = grepl("(^|,)\\s*LymphGen\\s*(,|$)",
                    ifelse(is.na(cnv_features_grch37$feature_set), "",
                           cnv_features_grch37$feature_set))

stopifnot(
  all(cnv_features_grch37$direction %in% c("gain", "loss")),
  "chrom must be 1-22, X or Y" = all(cnv_features_grch37$chrom %in% c(1:22, "X", "Y")),
  all(cnv_features_grch37$start < cnv_features_grch37$end),
  all(stats::na.omit(cnv_features_grch37$background_arm) %in%
        c("none", cnv_features_grch37$region_name[cnv_features_grch37$type == "arm"])),
  # DLBClass backgrounds come only from the rule applied within the DLBClass
  # set (background = "auto"); an override here would change the DLBClass GSM
  "DLBClass rows must have an empty background_arm" =
    all(is.na(cnv_features_grch37$background_arm[in_dlbclass])),
  all(cnv_features_grch37$type %in% c("arm", "focal", "gene", "chromosome")), # Make sure only allowed 
  # lymphgen_level: the copy number change (in copies) an event must reach
  # to count for LymphGen, 1 or 2; only on LymphGen rows
  all(cnv_features_grch37$lymphgen_level %in% c(1L, 2L, NA)),
  "lymphgen_level only on LymphGen rows" =
    all(is.na(cnv_features_grch37$lymphgen_level[!in_lymphgen])),
  # LymphGen features are never measured against an arm
  "LymphGen rows must have background_arm none" =
    all(cnv_features_grch37$background_arm[in_lymphgen] %in% "none"),
  # core: the default set; retired rows (include FALSE) can't be core
  all(cnv_features_grch37$core %in% c(TRUE, FALSE)),
  "retired rows (include FALSE) can't be core" =
    !any(cnv_features_grch37$core & !cnv_features_grch37$include),
  # focality limit: a positive length
  all(is.na(cnv_features_grch37$max_event_length) | cnv_features_grch37$max_event_length > 0)
)

usethis::use_data(cnv_features_grch37, overwrite = TRUE)
