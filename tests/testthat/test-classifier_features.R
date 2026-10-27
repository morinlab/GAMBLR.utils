maf_row = function(sample, gene, vc, change, chrom = "1", pos = 1){
  data.frame(Tumor_Sample_Barcode = sample, Hugo_Symbol = gene,
             Variant_Classification = vc, HGVSp_Short = change,
             Chromosome = chrom, Start_Position = pos)
}

test_that("registry tables join and DLBClass is complete", {
  dlbclass = get_classifier_features("DLBClass")
  expect_equal(length(unique(dlbclass$classifier_feature)), 163)
  expect_true(all(dlbclass$feature %in% genetic_features_grch37$feature))
  expect_false(anyDuplicated(genetic_features_grch37$feature) > 0)
  expect_error(get_classifier_features("nope"), "unknown classifier")
})

test_that("classify_mutations assigns classes", {
  maf = rbind(maf_row("a", "X", "Nonsense_Mutation", "p.Q1*"),
              maf_row("a", "X", "Splice_Region", ""),
              maf_row("a", "X", "Splice_Region", "p.X5_splice"),
              maf_row("a", "X", "Nonstop_Mutation", "p.*10Lext*5"),
              maf_row("a", "X", "Intron", ""))
  expect_equal(classify_mutations(maf)$mutation_class,
               c("trunc", NA, "splice_region", "other_nonsilent", NA))
})

test_that("DLBClass mutation and SV encoding", {
  meta = data.frame(sample_id = c("s1", "s2", "s3", "s4"),
                    bcl2_ba = c("POS", "NEG", NA, NA),
                    bcl6_ba = NA, myc_ba = NA)
  maf = rbind(
    maf_row("s1", "MYD88", "Missense_Mutation", "p.L265P", "3", 38182641),
    maf_row("s1", "MYD88", "Silent", "p.L10L", "3", 38180000),
    maf_row("s2", "MYD88", "Missense_Mutation", "p.S219C", "3", 38182000),
    maf_row("s2", "MYD88", "Silent", "p.L10L", "3", 38180000),
    maf_row("s3", "CREBBP", "Silent", "p.L10L", "16", 3800000)
  )
  sv = data.frame(tumour_sample_id = "s3", gene = "BCL2", partner = "IGH")
  ft = assemble_classifier_features(meta, classifier = "DLBClass",
                                    maf_data = maf, sv_data = sv)
  gsm = classifier_feature_matrix(ft, "DLBClass")
  expect_equal(colnames(gsm)[1:3], c("SV.BCL6", "DTX1", "CD70"))
  expect_equal(gsm[, "MYD88.L265P"], c(2, 0, 0, 0))
  # a silent mutation doesn't replace a non-silent one (unlike maf2gsm.py)
  expect_equal(gsm[, "MYD88.OTHER"], c(1, 2, 0, 0))
  expect_equal(gsm[, "MYD88"], c(2, 2, 0, 0))
  expect_equal(gsm[, "CREBBP"], c(0, 0, 1, 0))
  # SV from metadata (s1) or sv_data (s3)
  expect_equal(gsm[, "SV.BCL2"], c(3, 0, 3, 0))
  long = classifier_feature_matrix(ft, "DLBClass", output_format = "long",
                                   unassessed_value = NA)
  bcl2 = long[long$classifier_feature == "SV.BCL2", ]
  expect_equal(bcl2$assessed, c(TRUE, TRUE, TRUE, TRUE))
  expect_true(all(is.na(long$value[long$classifier_feature == "X17P.DEL"])))
  expect_equal(ft$source[ft$feature == "BCL2_SV"], c("metadata", NA, "sv_data", NA))
})

test_that("LymphGen features and position rules", {
  lymphgen = get_classifier_features("LymphGen")
  expect_equal(length(unique(lymphgen$classifier_feature)), 115)
  expect_setequal(unique(lymphgen$group), c("BN2", "EZB", "MCD", "N1", "ST2", "A53"))

  meta = data.frame(sample_id = c("s1", "s2", "s3"),
                    ssm_available = c(TRUE, TRUE, FALSE))
  maf = rbind(
    # NOTCH1: truncations in the PEST domain only
    maf_row("s1", "NOTCH1", "Nonsense_Mutation", "p.Q2440*", "9", 139390000),
    maf_row("s2", "NOTCH1", "Nonsense_Mutation", "p.Q100*", "9", 139420000),
    maf_row("s2", "NOTCH1", "Missense_Mutation", "p.R2400W", "9", 139390500),
    maf_row("s1", "MYD88", "Missense_Mutation", "p.L265P", "3", 38182641),
    # EZH2 outside the SET domain
    maf_row("s2", "EZH2", "Missense_Mutation", "p.A100T", "7", 148530000),
    # PIM1 silent mutations near (s1) and far from (s2) the TSS
    maf_row("s1", "PIM1", "Silent", "p.L5L", "6", 37138500),
    maf_row("s2", "PIM1", "Silent", "p.L300L", "6", 37142500)
  )
  ft = assemble_classifier_features(meta, classifier = "LymphGen",
                                    maf_data = maf, sv_from_metadata = NULL)
  lg = classifier_feature_matrix(ft, "LymphGen", unassessed_value = NA)
  expect_equal(lg[, "NOTCH1:MUTATION"], c(1, 0, NA))
  expect_equal(lg[, "MYD88L265P:MUTATION"], c(1, 0, NA))
  expect_equal(lg[, "EZH2:SubMUTATION"], c(0, 0, NA))
  expect_equal(lg[, "PIM1:Synon"], c(1, 0, NA))
  # 5' UTR features need noncoding_maf; Synon is still assessed from the MAF
  expect_false(any(ft$assessed[grepl("_utr5_", ft$feature)]))

  # without the SET and tss4kb rules
  ft = assemble_classifier_features(meta, classifier = "LymphGen",
                                    maf_data = maf, sv_from_metadata = NULL,
                                    position_rules = c("PEST", "ITAM"))
  lg = classifier_feature_matrix(ft, "LymphGen", unassessed_value = NA)
  expect_equal(lg[, "EZH2:SubMUTATION"], c(0, 1, NA))
  expect_equal(lg[, "PIM1:Synon"], c(1, 1, NA))
  expect_equal(lg[, "NOTCH1:MUTATION"], c(1, 0, NA))
})
