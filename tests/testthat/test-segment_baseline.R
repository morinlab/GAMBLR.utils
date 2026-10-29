genome_seg = function(id, cn_first = 2, cn_rest = 2){
  arms = as.data.frame(GAMBLR.data::chromosome_arms_grch37)
  arms$chromosome = gsub("^chr", "", arms$chromosome)
  arms = arms[arms$chromosome %in% as.character(1:22), ]
  data.frame(ID = id, chrom = arms$chromosome, start = arms$start, end = arms$end,
             CN = c(rep(cn_first, 10), rep(cn_rest, nrow(arms) - 10)))
}

test_that("whole-genome segments pass the coverage check and get a baseline", {
  seg = rbind(genome_seg("diploid"), genome_seg("tetraploid", cn_first = 3, cn_rest = 4))
  cov = check_segment_coverage(seg)
  expect_true(all(cov$ok))
  expect_true(all(cov$coverage > 0.9))
  b = add_segment_baseline(seg)
  expect_equal(unique(b$baseline_cn[b$ID == "diploid"]), 2)
  expect_equal(unique(b$baseline_cn[b$ID == "tetraploid"]), 4)
  expect_equal(unique(b$baseline_method), "mode")
  m = add_segment_baseline(seg, method = "mean")
  expect_true(unique(m$baseline_cn[m$ID == "tetraploid"]) > 3 & unique(m$baseline_cn[m$ID == "tetraploid"]) < 4)
})

test_that("region-filtered segments fail with a clear error", {
  seg = genome_seg("s1")
  region = seg[seg$chrom == "3", ]
  expect_false(check_segment_coverage(region)$ok)
  expect_error(add_segment_baseline(region), "covering too little of the genome")
})

test_that("an existing baseline of the same method is kept, another is recomputed", {
  seg = genome_seg("s1")
  region = seg[seg$chrom == "3", ]
  region$baseline_cn = 4; region$baseline_method = "mode"
  expect_equal(add_segment_baseline(region)$baseline_cn, rep(4, nrow(region)))
  whole = seg; whole$baseline_cn = 99; whole$baseline_method = "mean"
  expect_message(r <- add_segment_baseline(whole, method = "mode"), "recomputing")
  expect_equal(unique(r$baseline_cn), 2)
})

test_that("overlapping segments are not double counted in coverage", {
  seg = genome_seg("s1")
  doubled = rbind(seg, seg)
  expect_equal(check_segment_coverage(doubled)$coverage, check_segment_coverage(seg)$coverage)
})

test_that("check_segments reports coverage, build and overlap problems per sample", {
  good = genome_seg("good")
  sparse = genome_seg("sparse"); sparse = sparse[sparse$chrom == "3", ]
  hg38_like = genome_seg("hg38_like")
  q17 = which(hg38_like$chrom == "17")
  hg38_like$end[q17[length(q17)]] = 83197981  # the 17q segment past the grch37 end
  dup = genome_seg("dup"); dup = rbind(dup, dup[1, ])
  qc = check_segments(rbind(good, sparse, hg38_like, dup))
  rownames(qc) = qc$ID
  expect_true(qc["good", "ok"])
  expect_false(qc["sparse", "coverage_ok"])
  expect_false(qc["hg38_like", "build_ok"])
  expect_equal(qc["hg38_like", "n_beyond_chrom_end"], 1L)
  expect_false(qc["dup", "overlap_ok"])
  expect_equal(sort(qc$ID[!qc$ok]), c("dup", "hg38_like", "sparse"))
})

test_that("the baseline error lists the samples and how to remove them", {
  seg = genome_seg("s1"); seg = seg[seg$chrom == "3", ]
  expect_error(add_segment_baseline(seg), "check_segments")
})

test_that("dummy (placeholder) rows never count towards baseline, ploidy or coverage", {
  seg = rbind(data.frame(ID = "s1", chrom = as.character(1:22), start = 1, end = 4e7, CN = 3,
                         dummy_segment = 0),
              data.frame(ID = "s1", chrom = as.character(1:22), start = 4e7 + 1, end = 1.4e8,
                         CN = 2, dummy_segment = 1))
  no_flag = seg[, setdiff(colnames(seg), "dummy_segment")]
  expect_lt(check_segment_coverage(seg)$coverage, check_segment_coverage(no_flag)$coverage)
  expect_equal(segmented_data_to_ploidy_estimate(seg)$ploidy, 3)
  expect_equal(segmented_data_to_ploidy_estimate(no_flag)$ploidy, 2)
  expect_equal(unique(add_segment_baseline(seg, min_genome_coverage = 0.1)$baseline_cn), 3)
})

test_that("the three baseline methods differ as documented and are shared", {
  seg = genome_seg("s1", cn_first = 3.3, cn_rest = 1.9)   # neutral at 1.9
  expect_equal(unique(add_segment_baseline(seg)$baseline_cn), 1.9)
  expect_equal(unique(add_segment_baseline(seg, "mode_rounded")$baseline_cn), 2)
  expect_gt(unique(add_segment_baseline(seg, "mean")$baseline_cn), 1.9)
  for(m in baseline_methods()){
    expect_equal(segmented_data_to_ploidy_estimate(seg, method = m)$ploidy,
                 unique(add_segment_baseline(seg, m)$baseline_cn))
  }
  # a "mode" baseline converts to "mode_rounded" without recomputing (works on a region)
  region = add_segment_baseline(seg)[seg$chrom == "3", ]
  expect_equal(unique(add_segment_baseline(region, "mode_rounded")$baseline_cn), 2)
  expect_error(add_segment_baseline(region, "mean"), "covering too little")
})

test_that("ploidy estimates and cn matrices use the chosen method and existing baselines", {
  seg = genome_seg("s1", cn_first = 3.3, cn_rest = 1.9)
  # an existing genome-wide baseline is used on region segments, not re-estimated
  region = add_segment_baseline(seg)[seg$chrom == "3", ]
  region$CN = 3; region$log.ratio = log2(3 / 2)
  expect_equal(segmented_data_to_ploidy_estimate(region)$ploidy, 1.9)
  expect_equal(segmented_data_to_ploidy_estimate(region, method = "mode_rounded")$ploidy, 2)
  expect_message(segmented_data_to_ploidy_estimate(region, method = "mean"), "estimated from the segments")
  bins = data.frame(chrom = "3", start = 1e6, end = 2e6, name = "bin1")
  cn = function(m) segmented_data_to_cn_matrix(region, strategy = "custom_regions", regions = bins,
                                               genome_build = "grch37", adjust_for_ploidy = TRUE,
                                               baseline_method = m, rounded = FALSE)
  expect_equal(unname(cn("mode")[1, 1]), 3.1)
  expect_equal(unname(cn("mode_rounded")[1, 1]), 3)
})

test_that("only autosomes 1-22 count, and copy number is capped at 20, everywhere", {
  seg = genome_seg("s1")
  extra = data.frame(ID = "s1", chrom = c("X", "4_GL000008v2_random"), start = 1, end = 5e9, CN = 1)
  seg = rbind(seg, extra)
  expect_equal(segmented_data_to_ploidy_estimate(seg)$ploidy, 2)
  expect_equal(unique(add_segment_baseline(seg)$baseline_cn), 2)
  amp = genome_seg("s2", cn_first = 50)
  expect_equal(segmented_data_to_ploidy_estimate(amp, method = "mean")$ploidy,
               unique(add_segment_baseline(amp, "mean")$baseline_cn))
})

test_that("samples without a baseline stop the functions that need one, naming them", {
  seg = rbind(genome_seg("s1"), genome_seg("s2"))
  seg = add_segment_baseline(seg)
  seg$baseline_cn[seg$ID == "s2"] = NA
  expect_error(add_segment_baseline(seg), "s2")
  expect_error(center_segments(seg), "no baseline_cn")
  expect_error(segmented_data_to_ploidy_estimate(seg), "no baseline_cn")
  expect_equal(unique(center_segments(seg[seg$ID == "s1", ])$centre_cn), 2)
  # a sample whose segments have no CN can't get one either
  no_cn = genome_seg("s3"); no_cn$CN = NA
  expect_error(add_segment_baseline(rbind(genome_seg("s1"), no_cn)), "s3")
})

test_that("controlfreec samples get a baseline of 2 without the coverage check", {
  cf = genome_seg("cf", cn_first = 3, cn_rest = 3)[1:6, ]   # sparse, gains only
  cf$module = "controlfreec"
  bb = genome_seg("bb", cn_first = 3, cn_rest = 4); bb$module = "battenberg"
  seg = rbind(cf, bb)
  for(m in baseline_methods()){
    b = add_segment_baseline(seg, m)
    expect_equal(unique(b$baseline_cn[b$ID == "cf"]), 2)
    p = segmented_data_to_ploidy_estimate(seg, method = m)
    expect_equal(p$ploidy[p$sample_id == "cf"], 2)
  }
  expect_equal(unique(add_segment_baseline(seg)$baseline_cn[seg$ID == "bb"]), 4)
  expect_equal(unique(center_segments(cf)$CN), 3)
  expect_true(check_segments(seg)$coverage_ok[check_segments(seg)$ID == "cf"])
})
