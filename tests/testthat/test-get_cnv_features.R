test_that("named sets come back complete, with unique names", {
  for(s in c("DLBClass", "LymphGen")){
    f = get_cnv_features("grch37", feature_set = s)
    expect_false(any(duplicated(f$name)))
    expect_true(all(grepl(s, f$feature_set)))
  }
  expect_equal(nrow(get_cnv_features("grch37", feature_set = "DLBClass")), 61)
  expect_false(any(grepl("_LG$", cnv_features_grch37$name)))
})

test_that("the core set, all, and curated", {
  core = get_cnv_features("grch37")
  expect_false(any(duplicated(core$name)))
  expect_true(all(cnv_features_grch37$core[match(paste(core$name, core$feature_set),
                                                   paste(cnv_features_grch37$name, cnv_features_grch37$feature_set))]))
  # the starting core: the curated rows and the non-retired DLBClass rows
  expect_true(all(is.na(core$feature_set) | core$feature_set == "DLBClass"))
  expect_equal(nrow(core), 68)
  all = get_cnv_features("grch37", feature_set = "all")
  expect_false(any(duplicated(all$name)))
  expect_true(all(core$name %in% all$name))
  expect_gt(nrow(all), nrow(core))
  # retired rows are in neither, but stay in their named set
  retired = cnv_features_grch37$name[!cnv_features_grch37$include]
  expect_false(any(retired %in% all$name))
  expect_true(all(retired %in% get_cnv_features("grch37", feature_set = "DLBClass")$name))
  expect_true(all(is.na(get_cnv_features("grch37", feature_set = "curated")$feature_set)))
})

test_that("rows are chosen by priority curated > LymphGen > DLBClass", {
  all = get_cnv_features("grch37", feature_set = "all")
  # in the curated set and LymphGen: the curated row
  expect_true(is.na(all$feature_set[all$name == "SPIB_gain"]))
  # in LymphGen and DLBClass: the LymphGen row
  expect_equal(all$feature_set[all$name == "6q_loss"], "LymphGen")
  g = get_cnv_features("grch37", feature_set = "all", priority = c("DLBClass", "curated", "LymphGen"))
  expect_equal(g$feature_set[g$name == "6q_loss"], "DLBClass")
  expect_false(any(duplicated(g$name)))
})

test_that("arm boundaries are the same in every set", {
  f = cnv_features_grch37[cnv_features_grch37$type == "arm", ]
  key = paste(f$chrom, f$region_name)
  expect_true(all(tapply(paste(f$start, f$end), key, function(v) length(unique(v))) == 1))
})
