gene_feature = function(name, chrom, start, end, max_event_length = 30e6){
  data.frame(name = name, type = "gene", region_name = name, chrom = chrom,
             start = start, end = end, direction = "loss", min_fraction = 0.01,
             min_coverage = 0.5, background_arm = "none",
             max_event_length = max_event_length)
}
# grch37 coordinates
cdkn2a = gene_feature("CDKN2A_loss", "9", 21967750, 21995300)
tp53 = gene_feature("TP53_loss", "17", 7571719, 7590868)

# one sample's chromosome: CN 2 except for the given pieces (start, end, CN)
chrom_segs = function(id, chrom, chrom_end, pieces){
  breaks = sort(unique(c(0, chrom_end, unlist(lapply(pieces, `[`, 1:2)))))
  segs = data.frame(ID = id, chrom = chrom, start = utils::head(breaks, -1),
                    end = breaks[-1], CN = 2)
  for(p in pieces){
    segs$CN[segs$start >= p[1] & segs$end <= p[2]] = p[3]
  }
  segs
}
chrom9 = function(id, pieces) chrom_segs(id, "9", 141213431, pieces)

call_loss = function(segs, feature){
  out = segmented_data_to_cnv_status(segs, cnv_features = feature,
                                     genome_build = "grch37",
                                     value_scale = "CN", centering = "none")
  stats::setNames(out[[feature$name]], rownames(out))
}

test_that("max_event_length keeps only changes no longer than the limit", {
  segs = rbind(
    chrom9("arm_loss", list(c(0, 49e6, 1))),
    chrom9("focal_loss", list(c(20e6, 25e6, 1))),
    chrom9("loss_25Mb", list(c(9.5e6, 34.5e6, 1))),
    chrom9("loss_35Mb", list(c(5e6, 40e6, 1))),
    chrom9("homdel_in_arm_loss", list(c(0, 21.9e6, 1), c(21.9e6, 22.1e6, 0),
                                      c(22.1e6, 49e6, 1))),
    chrom9("arm_homdel", list(c(0, 49e6, 0)))
  )
  calls = call_loss(segs, cdkn2a)
  expect_equal(calls[["arm_loss"]], 0)
  expect_equal(calls[["focal_loss"]], 1)
  expect_equal(calls[["loss_25Mb"]], 1)
  expect_equal(calls[["loss_35Mb"]], 0)
  expect_equal(calls[["homdel_in_arm_loss"]], 2)
  # too long at level 2 and at level 1
  expect_equal(calls[["arm_homdel"]], 0)

  # without the limit every loss is called
  no_limit = cdkn2a
  no_limit$max_event_length = NA
  expect_true(all(call_loss(segs, no_limit) > 0))
})

test_that("segments split by liftOver give the same extent", {
  # a 9p loss in 5 Mb pieces with uncovered 100 kb gaps between them
  pieces = lapply(seq(0, 45e6, by = 5e6), function(s){
    data.frame(ID = "split_arm_loss", chrom = "9", start = s,
               end = min(s + 4.9e6, 49e6), CN = 1)
  })
  segs = rbind(do.call(rbind, pieces),
               data.frame(ID = "split_arm_loss", chrom = "9", start = 49e6,
                          end = 141213431, CN = 2))
  expect_equal(call_loss(segs, cdkn2a)[["split_arm_loss"]], 0)
})

test_that("short arms count, and the centromere doesn't stop the extent", {
  # 17p is about 22 Mb, so a whole-17p loss counts for TP53 (as in LymphGen's
  # training data); a whole-chromosome-17 loss does not
  centromere = c(22.2e6, 25.8e6)
  segs = rbind(
    data.frame(ID = "loss_17p", chrom = "17", start = c(0, centromere[2]),
               end = c(centromere[1], 81195210), CN = c(1, 2)),
    data.frame(ID = "loss_chr17", chrom = "17", start = c(0, centromere[2]),
               end = c(centromere[1], 81195210), CN = c(1, 1))
  )
  calls = call_loss(segs, tp53)
  expect_equal(calls[["loss_17p"]], 1)
  expect_equal(calls[["loss_chr17"]], 0)
})

test_that("LymphGen_model gene rows have the 30 Mb limit", {
  model = get_cnv_features(feature_set = "LymphGen_model")
  expect_true(all(model$max_event_length[model$type == "gene"] == 30e6))
  expect_true(all(is.na(model$max_event_length[model$type == "arm"])))
  expect_true(all(is.na(get_cnv_features(feature_set = "DLBClass")$max_event_length)))
})

test_that("a change at another level ends the extent", {
  # a 10 Mb single-copy loss over CDKN2A next to a 25 Mb homozygous deletion:
  # each counts on its own length, as in LymphGen
  segs = chrom9("hetloss_next_to_homdel", list(c(15e6, 25e6, 1), c(25e6, 50e6, 0)))
  expect_equal(call_loss(segs, cdkn2a)[["hetloss_next_to_homdel"]], 1)
  # a gain next to an amplification
  gain = cdkn2a
  gain$name = "CDKN2A_gain"
  gain$direction = "gain"
  segs = chrom9("gain_next_to_amp", list(c(15e6, 25e6, 3), c(25e6, 50e6, 5)))
  out = segmented_data_to_cnv_status(segs, cnv_features = gain, genome_build = "grch37",
                                     value_scale = "CN", centering = "none")
  expect_equal(out$CDKN2A_gain, 1)
})
