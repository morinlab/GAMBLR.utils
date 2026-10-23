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
