#!/usr/bin/env Rscript
# Seed for the classifier feature tables, kept for provenance. Run from the
# GAMBLR.utils root. It writes
#   inst/extdata/genetic_features.grch37.tsv  (shared feature definitions)
#   inst/extdata/classifier_features.tsv      (classifier feature map)
# and replaces the "LymphGen_model" rows of inst/extdata/cnv_features.grch37.tsv.
# After seeding, the TSVs are the master copies: edit them directly and run
# data-raw/cnv_features.R and data-raw/classifier_features.R.
#
# Sources (local copies; not distributed with the package):
#   DLBClass: DLBclass-tool maf2gsm.py, sv2gsm.py and
#     gsm/feature_order.19Aug2024.txt.
#   LymphGen: LymphGen 2, the version run in GAMBL ("Source LymphGen.R",
#     Full.Uber.2020.txt, Fullmat.2020.txt, Featord.2020.txt,
#     Study.samp.txt), with Version = 2 (its default). Model9() refits the
#     model from the training data each run; the features here are those its
#     nine submodels select with complete gene lists: the full model (copy
#     number and both fusions available), the A53 model, and the submodels
#     used when copy number and/or the BCL2 or BCL6 fusion status is missing.
#     classifier_features.tsv records which submodels use each. The position
#     rules and Synon windows are those LymphGen 2 applies.
suppressMessages({library(dplyr); library(readr)})

dlbclass_order_file = "~/git/DLBClass-tool/gsm/feature_order.19Aug2024.txt"
lymphgen_dir = "~/git/Lymphgen_2"

feature_order = readLines(dlbclass_order_file)
cnv = read_tsv("inst/extdata/cnv_features.grch37.tsv", col_types = cols(.default = "c"),
               na = c("", "NA"))

nonsilent = "trunc,missense,inframe,nonstop,splice_region,other_nonsilent"
rows = list(); map = list()
add = function(...) rows[[length(rows) + 1]] <<- tibble(...)
link = function(...) map[[length(map) + 1]] <<- tibble(...)

# DLBClass mutations ------------------------------------------------------
dlbclass_genes = feature_order[!grepl("^X|^SV[.]|^MYD88[.]", feature_order)]
for(g in dlbclass_genes){
  gene = sub("^HLA[.]", "HLA-", g)
  add(feature = paste0(gene, "_nonsilent"), data_type = "ssm", gene = gene,
      variant_class = nonsilent, requires = "maf_data")
  add(feature = paste0(gene, "_silent"), data_type = "ssm", gene = gene,
      variant_class = "silent", requires = "maf_data")
  link(classifier = "DLBClass", classifier_feature = g,
       feature = paste0(gene, c("_nonsilent", "_silent")), value = c("2", "1"))
}
add(feature = "MYD88_L265P", data_type = "ssm", gene = "MYD88",
    variant_class = "missense", protein_change = "p.L265P", requires = "maf_data")
add(feature = "MYD88_nonsilent_other", data_type = "ssm", gene = "MYD88",
    variant_class = nonsilent, protein_change_exclude = "p.L265P", requires = "maf_data")
link(classifier = "DLBClass", classifier_feature = "MYD88.L265P",
     feature = "MYD88_L265P", value = "2")
link(classifier = "DLBClass", classifier_feature = "MYD88.OTHER",
     feature = c("MYD88_nonsilent_other", "MYD88_silent"), value = c("2", "1"))

# SVs (shared) ------------------------------------------------------------
for(g in c("BCL2", "BCL6", "MYC")){
  add(feature = paste0(g, "_SV"), data_type = "sv", gene = g,
      requires = "sv_data or sv_from_metadata")
  link(classifier = "DLBClass", classifier_feature = paste0("SV.", g),
       feature = paste0(g, "_SV"), value = "3")
}

# DLBClass CNV: reference the DLBClass rows of cnv_features ---------------
dcnv = filter(cnv, grepl("DLBClass", feature_set))
for(i in seq_len(nrow(dcnv))){
  f = paste0(dcnv$name[i], "_DLBClass")
  add(feature = f, data_type = "cnv", cnv_feature = dcnv$name[i],
      cnv_feature_set = "DLBClass", requires = "seg_data")
  link(classifier = "DLBClass", classifier_feature = paste0("X", dcnv$dlbclass_name[i]),
       feature = f, value = "level")
}

# LymphGen: features selected by Model9 ------------------------------------
# (the source file also defines Predict9 etc.; only Model9 is used)
suppressWarnings(source(file.path(lymphgen_dir, "Source LymphGen.R")))
rd = function(f) read.table(file.path(lymphgen_dir, f), header = TRUE, sep = "\t")
uber = rd("Full.Uber.2020.txt")
study = rd("Study.samp.txt")
# as in Run_Lymphgen.R
lab = study[, 18]; lab[lab == "SDT"] = "ST2"
featord = rd("Featord.2020.txt"); featord[, 11] = featord[, 11] == 1
# all nine submodels: Model9 picks one per sample by which of copy number and
# the BCL2 and BCL6 fusions are available, and they select different features
# (e.g. BCL2 Synon when BCL2 fusion status is missing)
submodel_names = c("Full", "NoCGH", "NoBCL2Fus", "NoBCL2CGH", "NoBCL6Fus",
                   "NoBCL6CGH", "NoFus", "NoFusCGH", "A53")
model = Model9(lab, uber, as.matrix(rd("Fullmat.2020.txt") == 1), CGH.class = 0,
               mut.genelist = unique(uber[, 3]), CGH.genelist = unique(uber[, 3]),
               Has.CGH = study[, 7] == 1, modelcheck = 1:9,
               Featord.2020 = featord, Version = 2)
selected = NULL
# Full and A53 first, so their features come first in feature_order
for(k in c(1, 9, 2:8)){
  mt = model[[k]]$Mod.type
  mt = mt[!is.na(mt$out) & mt$out > 0, ]   # drop padding rows
  # a selected feature combines up to four training features (out, ot1, ...)
  ix = as.matrix(mt[, grep("^(out|ot1)", names(mt))])
  subtype = if(k == 9) rep("A53", nrow(mt)) else model[[k]]$lst[mt$i]
  for(r in seq_len(nrow(mt))){
    parts = uber[match(stats::na.omit(ix[r, ]), uber$Feature.Index), ]
    selected = rbind(selected, data.frame(submodel = submodel_names[k],
                                          subtype = subtype[r], rank = r,
                                          gene = parts$Gene.Symbol,
                                          entrez = parts$ENTREZ.ID,
                                          type = parts$Feature.type,
                                          # LymphGen codes X and Y as 23 and 24
                                          chrom = sub("^23$", "X", sub("^24$", "Y", parts$Chrom.2018)),
                                          start = parts$Gene.start,
                                          end = parts$Gene.end))
  }
}

# grch37 symbols, where LymphGen uses newer ones
symbol = c(MYD88L265P = "MYD88", BCL6FusionCatcherORInHouse = "BCL6",
           BCL2TransFISHMUTFusion = "BCL2")
selected$gene37 = ifelse(selected$gene %in% names(symbol), symbol[selected$gene], selected$gene)
arm = grepl("^[0-9XY]+[pq]$", selected$gene)

# Synon windows: those LymphGen 2 applies (Feat.start to Feat.end of its
# Synon features, grch37, ends inclusive): the first 4 kb of the gene from
# its TSS, downstream only, and no further than the gene end. LymphGen 2 uses
# its own gene coordinates, so its TSS can differ from get_gene_tss() (e.g.
# ETS1, by 65 kb, and BCL2 and BCL2L1, by about 1 kb).
synon2 = distinct(filter(uber, grepl("Syn", Feature.type)),
                  ENTREZ.ID, Chrom.2018, Feat.start, Feat.end)
stopifnot(!anyDuplicated(synon2$ENTREZ.ID))
tss_window = function(gene, entrez){
  w = synon2[synon2$ENTREZ.ID == entrez, ]
  if(nrow(w) == 0) stop("no LymphGen 2 Synon window for ", gene)
  list(chrom = sub("^23$", "X", sub("^24$", "Y", w$Chrom.2018)),
       start = w$Feat.start, end = w$Feat.end)
}

# position rules as LymphGen 2 applies them (hg19, mutation start; a mutation
# is dropped past the given position): NOTCH1 and NOTCH2 truncations up to
# 139391455 and 120459150, CD79B mutations up to 62007172 and truncations up
# to 62006800. Its EZH2 (SET) rule uses Entrez 2416 rather than 2146, so it
# never takes effect; the SET row is kept but isn't applied by default.
rule_range = list(PEST_NOTCH1 = c("9", 1, 139391455), PEST_NOTCH2 = c("1", 1, 120459150),
                  ITAM_trunc = c("17", 1, 62006800), ITAM_missense = c("17", 1, 62007172),
                  SET = c("7", 148506238, 148508764))
ssm_feature = function(gene, entrez, kind){
  # kind: "trunc", "mutation", "silent_tss4kb" or "utr5_tss4kb"; returns the
  # shared feature name(s), adding definitions as needed
  out = character(0)
  restricted = list(NOTCH1 = "PEST_NOTCH1", NOTCH2 = "PEST_NOTCH2")
  if(kind %in% c("trunc", "mutation") && gene %in% c("NOTCH1", "NOTCH2")){
    # NOTCH1/NOTCH2: truncations only (non-truncating mutations are dropped)
    r = rule_range[[restricted[[gene]]]]
    f = paste0(gene, "_trunc_PEST")
    add(feature = f, data_type = "ssm", gene = gene, entrez = as.character(entrez),
        variant_class = "trunc", chrom = r[1], start = as.integer(r[2]),
        end = as.integer(r[3]), position_rule = "PEST", requires = "maf_data",
        summary = "LymphGen: truncations only, up to the given position (PEST domain)")
    return(f)
  }
  if(kind %in% c("trunc", "mutation") && gene == "CD79B"){
    for(t in if(kind == "trunc") "trunc" else c("trunc", "missense")){
      r = rule_range[[paste0("ITAM_", t)]]
      f = paste0("CD79B_", if(t == "trunc") "trunc" else "missense_inframe", "_ITAM")
      add(feature = f, data_type = "ssm", gene = gene, entrez = as.character(entrez),
          variant_class = if(t == "trunc") "trunc" else "missense,inframe,nonstop",
          chrom = r[1], start = as.integer(r[2]), end = as.integer(r[3]),
          position_rule = "ITAM", requires = "maf_data",
          summary = "LymphGen: up to the given position (ITAM region)")
      out = c(out, f)
    }
    return(out)
  }
  if(kind %in% c("trunc", "mutation") && gene == "EZH2"){
    r = rule_range$SET
    f = paste0("EZH2_", kind, "_SET")
    add(feature = f, data_type = "ssm", gene = gene, entrez = as.character(entrez),
        variant_class = if(kind == "trunc") "trunc" else "trunc,missense,inframe,nonstop",
        chrom = r[1], start = as.integer(r[2]), end = as.integer(r[3]),
        position_rule = "SET", requires = "maf_data",
        summary = "LymphGen: catalytic (SET) domain (not applied by the LymphGen code, which uses the wrong Entrez ID)")
    return(f)
  }
  if(kind == "trunc"){
    f = paste0(gene, "_trunc")
    add(feature = f, data_type = "ssm", gene = gene, entrez = as.character(entrez),
        variant_class = "trunc", requires = "maf_data")
    return(f)
  }
  if(kind == "mutation"){
    f = paste0(gene, "_mutation")
    add(feature = f, data_type = "ssm", gene = gene, entrez = as.character(entrez),
        variant_class = "trunc,missense,inframe,nonstop", requires = "maf_data",
        summary = "LymphGen MUTATION: any truncating or non-truncating coding mutation")
    return(f)
  }
  w = tss_window(gene, entrez)
  cls = if(kind == "silent_tss4kb") "silent" else "utr5"
  f = paste0(gene, "_", kind)
  add(feature = f, data_type = "ssm", gene = gene, entrez = as.character(entrez),
      variant_class = cls, chrom = w$chrom, start = as.integer(w$start),
      end = as.integer(w$end), position_rule = "tss4kb",
      requires = if(cls == "silent") "maf_data" else "noncoding_maf",
      summary = "LymphGen Synon: first 4 kb of the gene from the TSS, within the gene (the LymphGen 2 window)")
  f
}

# copy number rows for cnv_features (set "LymphGen_model"), one per gene or
# arm and direction; HETLOSS/HOMDEL and GAIN/AMP share a row and differ by
# min_level in the map
arm_rows = filter(cnv, type == "arm", is.na(feature_set) | grepl("DLBClass|LymphGen", feature_set)) %>%
  distinct(region_name, .keep_all = TRUE)
cytobands = GAMBLR.data::cytobands_grch37
cnv_new = list()
cn_feature = function(gene, entrez, direction, chrom, start, end, is_arm){
  name = paste0(gene, "_", direction)
  if(is_arm){
    a = filter(arm_rows, region_name == gene)
    if(nrow(a) == 0) stop("no arm coordinates for ", gene)
    chrom = a$chrom[1]; start = as.numeric(a$start[1]); end = as.numeric(a$end[1])
    region = gene
  }else{
    cb = filter(cytobands, cb.chromosome == chrom, cb.start <= start, cb.end >= start)
    region = paste0(chrom, cb$cb.name[1])
  }
  cnv_new[[name]] <<- tibble(
    name = name, type = if(is_arm) "arm" else "gene", region_name = region,
    chrom = as.character(chrom), start = format(start, scientific = FALSE, trim = TRUE),
    end = format(end, scientific = FALSE, trim = TRUE),
    direction = direction,
    # arms: more than 80% of the arm (or a larger change); genes: any overlap
    min_fraction = if(is_arm) "0.8" else "0.01", min_coverage = "0.5",
    background_arm = "none", representative_gene = if(is_arm) NA_character_ else gene,
    feature_set = "LymphGen_model", include = "TRUE", core = "FALSE",
    citation = "LymphGen 2 R code (Source LymphGen.R, Full.Uber.2020.txt)",
    # genes: only focal changes (LymphGen counts changes of 30 Mb or less)
    max_event_length = if(is_arm) NA_character_ else "30000000",
    summary = if(is_arm) "LymphGen arm feature" else
      "LymphGen gene feature (LymphGen counts only focal changes of 30 Mb or less)")
  f = paste0(name, "_LymphGen")
  add(feature = f, data_type = "cnv", gene = if(is_arm) NA_character_ else gene,
      entrez = if(is_arm) NA_character_ else as.character(entrez),
      cnv_feature = name, cnv_feature_set = "LymphGen_model", requires = "seg_data")
  f
}

# one classifier feature per subtype, gene and types, wherever it's selected;
# a feature selected for two subtypes (HIST1H1B:SubSyn, ST2 and MCD) has a
# row for each but one feature_order
selected$row = paste(selected$submodel, selected$subtype, selected$rank)
types = tapply(selected$type, selected$row, paste, collapse = "|")
selected$key = paste(selected$subtype, selected$gene[match(selected$row, selected$row)],
                     types[selected$row])
lymphgen_orders = character(0)
for(key in unique(selected$key)){
  first_row = selected$row[selected$key == key][1]
  s = selected[selected$row == first_row, ]
  in_models = unique(selected$submodel[selected$key == key])
  submodels = paste(submodel_names[submodel_names %in% in_models], collapse = ",")
  gene = s$gene37[1]
  name = paste0(s$gene[1], ":", paste(s$type, collapse = "|"))
  lymphgen_orders = union(lymphgen_orders, name)
  lymphgen_order = match(name, lymphgen_orders)
  features = character(0); min_level = integer(0)
  for(j in seq_len(nrow(s))){
    for(t in strsplit(s$type[j], "/")[[1]]){
      t = sub("^Sub", "", t)
      if(s$gene[j] == "MYD88L265P"){
        new = "MYD88_L265P"; lvl = NA
      }else if(t %in% c("FUSION", "CompUp")){
        new = paste0(gene, "_SV"); lvl = NA
      }else if(t %in% c("TRUNC", "MUTATION")){
        new = ssm_feature(gene, s$entrez[j], tolower(t)); lvl = NA
      }else if(t %in% c("Syn", "Synon")){
        new = c(ssm_feature(gene, s$entrez[j], "mutation"),
                ssm_feature(gene, s$entrez[j], "silent_tss4kb"),
                ssm_feature(gene, s$entrez[j], "utr5_tss4kb")); lvl = NA
      }else if(t %in% c("HETLOSS", "HOMDEL", "GAIN", "AMP")){
        direction = if(t %in% c("HETLOSS", "HOMDEL")) "loss" else "gain"
        new = cn_feature(gene, s$entrez[j], direction, s$chrom[j], s$start[j], s$end[j], arm[match(s$gene[j], selected$gene)])
        lvl = if(t %in% c("HOMDEL", "AMP")) 2L else 1L
      }else{
        stop("unknown LymphGen feature type ", t)
      }
      features = c(features, new); min_level = c(min_level, rep(lvl, length(new)))
    }
  }
  # a feature listed twice (e.g. HETLOSS and HOMDEL) keeps the lower level
  keep = !duplicated(data.frame(features, min_level))
  features = features[keep]; min_level = min_level[keep]
  lowest = tapply(min_level, features, function(x) if(all(is.na(x))) NA else min(x, na.rm = TRUE))
  features = unique(features)
  link(classifier = "LymphGen", classifier_feature = name, group = s$subtype[1],
       feature_order = lymphgen_order, feature = features, value = "1",
       min_level = as.integer(lowest[features]), submodels = submodels)
}

# outputs -----------------------------------------------------------------
features = bind_rows(rows) %>% distinct(feature, .keep_all = TRUE) %>%
  select(feature, data_type, gene, entrez, variant_class, protein_change,
         protein_change_exclude, chrom, start, end, position_rule, cnv_feature,
         cnv_feature_set, requires, summary) %>%
  mutate(start = as.integer(start), end = as.integer(end))
classifier_map = bind_rows(map)
d = classifier_map$classifier == "DLBClass"
classifier_map$feature_order[d] = match(classifier_map$classifier_feature[d], feature_order)
classifier_map = classifier_map %>%
  select(classifier, classifier_feature, group, feature_order, feature, value, min_level,
         submodels) %>%
  arrange(classifier, feature_order)
stopifnot(setequal(classifier_map$classifier_feature[d], feature_order),
          all(classifier_map$feature %in% features$feature),
          !anyDuplicated(unique(classifier_map[, c("classifier", "classifier_feature", "feature_order")])[, 1:2]))

write_tsv(features, "inst/extdata/genetic_features.grch37.tsv", na = "")
write_tsv(classifier_map, "inst/extdata/classifier_features.tsv", na = "")

cnv_out = bind_rows(filter(cnv, is.na(feature_set) | feature_set != "LymphGen_model"),
                    bind_rows(cnv_new))
write_tsv(cnv_out, "inst/extdata/cnv_features.grch37.tsv", na = "")
cat(nrow(features), "features;", nrow(classifier_map), "map rows;",
    length(cnv_new), "LymphGen_model cnv rows\n")
print(classifier_map %>% group_by(classifier, group) %>%
        summarise(n = n_distinct(classifier_feature), .groups = "drop"))
