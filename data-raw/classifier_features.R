#!/usr/bin/env Rscript
# Build genetic_features_grch37 and classifier_features from the master
# tables inst/extdata/genetic_features.grch37.tsv (one row per shared
# feature definition) and inst/extdata/classifier_features.tsv (which
# features each classifier uses, and how it encodes them). Edit the TSVs and
# rerun this script.
#
# The two tables join on `feature`: a feature is defined once and any number
# of classifiers can use it, so no definitions are duplicated across sets.
# Copy number features point at rows of cnv_features_grch37 (cnv_feature,
# cnv_feature_set) rather than repeating their coordinates.
#
# Seeded from:
#   DLBClass: maf2gsm.py, sv2gsm.py and gsm/feature_order.19Aug2024.txt in
#     DLBclass-tool (163 features, GSM order kept in feature_order), with the
#     MYD88.OTHER encoding taken as intended (a silent MYD88 mutation doesn't
#     replace a non-silent one) rather than as maf2gsm.py implements it.
#   LymphGen: the features its R code (Model9) selects in each of its nine
#     submodels (full, A53, and those for missing copy number or BCL2/BCL6
#     translocation status; the submodels column records which use each),
#     with the position rules and Synon windows LymphGen 2 applies (NOTCH1/
#     NOTCH2 PEST, CD79B ITAM, and Synon in the first 4 kb of the gene from
#     its TSS) and the EZH2 SET rule from the LymphGen instructions, each
#     tagged in position_rule so it can be switched off
#     (assemble_classifier_features(position_rules = ...)).
# See data-raw/classifier_features_seed.R for how the tables were seeded.
#
# TODO(hg38): grch37 coordinates only.

library(usethis)
library(readr)
# for .mutation_classes
devtools::load_all(quiet = TRUE)

genetic_features_grch37 = as.data.frame(read_tsv(
  "inst/extdata/genetic_features.grch37.tsv",
  col_types = cols(start = "i", end = "i", .default = "c"),
  na = c("", "NA")
))

classifier_features = as.data.frame(read_tsv(
  "inst/extdata/classifier_features.tsv",
  col_types = cols(feature_order = "i", min_level = "i", .default = "c"),
  na = c("", "NA")
))

load("data/cnv_features_grch37.rda")

f = genetic_features_grch37
m = classifier_features
classes = trimws(unlist(strsplit(stats::na.omit(f$variant_class), ",")))
cnv_rows = f$data_type == "cnv"
cnv_key = paste(f$cnv_feature, f$cnv_feature_set)[cnv_rows]
cnv_sets = ifelse(is.na(cnv_features_grch37$feature_set), "curated",
                  cnv_features_grch37$feature_set)
cnv_known = unlist(lapply(seq_len(nrow(cnv_features_grch37)), function(i){
  paste(cnv_features_grch37$name[i], trimws(strsplit(cnv_sets[i], ",")[[1]]))
}))
positioned = !is.na(f$chrom) | !is.na(f$start) | !is.na(f$end)

stopifnot(
  "feature names must be unique" = !anyDuplicated(f$feature),
  all(f$data_type %in% c("ssm", "sv", "cnv")),
  "unknown variant_class" = all(classes %in% .mutation_classes),
  "ssm features need a gene and variant_class" =
    all(!is.na(f$gene[f$data_type == "ssm"]) & !is.na(f$variant_class[f$data_type == "ssm"])),
  "sv features need a gene" = all(!is.na(f$gene[f$data_type == "sv"])),
  "cnv features must name a cnv_features row and its set" =
    all(cnv_key %in% cnv_known),
  "positions need chrom, start and end" =
    all(!is.na(f$chrom[positioned]) & !is.na(f$start[positioned]) &
          !is.na(f$end[positioned]) & f$start[positioned] <= f$end[positioned]),
  "map rows must name a known feature" = all(m$feature %in% f$feature),
  "value must be a number or \"level\"" =
    all(m$value == "level" | !is.na(suppressWarnings(as.numeric(m$value)))),
  "each classifier feature has one position in its order" =
    all(tapply(m$feature_order, paste(m$classifier, m$classifier_feature),
               function(x) length(unique(x))) == 1),
  "min_level is 1 or 2" = all(m$min_level %in% c(1L, 2L, NA)),
  "position_rule needs a position" = all(positioned[!is.na(f$position_rule)]),
  "unknown LymphGen submodel" =
    all(trimws(unlist(strsplit(stats::na.omit(m$submodels), ","))) %in%
          c("Full", "NoCGH", "NoBCL2Fus", "NoBCL2CGH", "NoBCL6Fus", "NoBCL6CGH",
            "NoFus", "NoFusCGH", "A53")),
  "each classifier feature has one set of submodels" =
    all(tapply(m$submodels, paste(m$classifier, m$classifier_feature),
               function(x) length(unique(x))) == 1)
)

usethis::use_data(genetic_features_grch37, overwrite = TRUE)
usethis::use_data(classifier_features, overwrite = TRUE)
