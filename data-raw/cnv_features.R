#!/usr/bin/env Rscript
# Build cnv_features_grch37 from the curated master table
# inst/extdata/cnv_features.grch37.tsv. Edit the TSV (add rows, or set
# include = FALSE to retire one) and rerun this script.
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
                   include = "l", .default = "c"),
  na = c("", "NA")
) %>%
  mutate(min_fraction = as.numeric(min_fraction),
         min_coverage = as.numeric(min_coverage)) %>%
  as.data.frame()

in_dlbclass = grepl("(^|,)\\s*DLBClass\\s*(,|$)",
                    ifelse(is.na(cnv_features_grch37$feature_set), "",
                           cnv_features_grch37$feature_set))

stopifnot(
  !any(duplicated(cnv_features_grch37$name)),
  all(cnv_features_grch37$direction %in% c("gain", "loss")),
  all(cnv_features_grch37$start < cnv_features_grch37$end),
  all(stats::na.omit(cnv_features_grch37$background_arm) %in%
        c("none", cnv_features_grch37$region_name[cnv_features_grch37$type == "arm"])),
  # DLBClass backgrounds come only from the rule applied within the DLBClass
  # set (background = "auto"); an override here would change the DLBClass GSM
  "DLBClass rows must have an empty background_arm" =
    all(is.na(cnv_features_grch37$background_arm[in_dlbclass]))
)

usethis::use_data(cnv_features_grch37, overwrite = TRUE)
