qc_segs = function(id, cn, chrom = "1", len = 1e7){
  start = 1 + c(0, cumsum(rep(len, length(cn) - 1)))
  data.frame(ID = id, chrom = chrom, start = start, end = start + len - 1, CN = cn,
             log.ratio = log2(cn / 2))
}

test_that("density, steps and altered fraction", {
  clean = qc_segs("clean", c(2, 2, 1, 2, 2))
  noisy = qc_segs("noisy", c(2, 2.1, 1.9, 2, 2.1), len = 5e5)
  qc = segmented_data_qc(rbind(clean, noisy))
  expect_equal(qc$sample_id, c("clean", "noisy"))
  expect_equal(qc$n_segments, c(5, 5))
  expect_equal(qc$short_segment_fraction, c(0, 1))
  # clean: steps 0, 1, 1, 0; noisy: 0.1, 0.2, 0.1, 0.1
  expect_equal(qc$median_step, c(0.5, 0.1))
  expect_equal(qc$baseline_mode, c(2, 2))
  expect_equal(qc$altered_fraction, c(0.2, 0))
  expect_true(is.na(qc$supplied_baseline[1]))
})

test_that("supplied baseline, CN/log ratio consistency, sex chromosomes", {
  ok = qc_segs("ok", c(2, 1, 3, 2, 4))
  # CN running against the log ratios, as from a failed absolute fit
  bad = qc_segs("bad", c(2, 1, 3, 2, 4)); bad$CN = 6 - bad$CN
  x = rbind(ok, bad); x$baseline_cn = 2
  x = rbind(x, transform(qc_segs("ok", 1), chrom = "X", baseline_cn = 2))
  qc = segmented_data_qc(x)
  expect_equal(qc$deep_loss_fraction, c(0, 0))
  # a failed fit: most of the genome at CN 0 while the log ratios are neutral
  failed = qc_segs("failed", c(0, 0, 0, 2, 3))
  failed$log.ratio = c(-0.1, 0, 0.05, 0, 0.6)
  f = segmented_data_qc(failed)
  expect_equal(f$deep_loss_fraction, 0.6)
  expect_equal(f$deep_loss_logratio, 0)
  # bad is mostly at CN 4 against a supplied baseline of 2
  expect_equal(qc$baseline_delta, c(0, -2))
  # chromosome X is left out by default
  expect_equal(qc$n_segments[1], 5)
  expect_equal(segmented_data_qc(x, autosomes_only = FALSE)$n_segments[1], 6)
})

test_that("log ratio scale and comparison with a second segment set", {
  a = qc_segs("s1", c(2, 2, 1, 2, 2))
  # same structure at a doubled baseline agrees fully; a shifted loss doesn't
  b = transform(a, CN = CN * 2)
  c = qc_segs("s1", c(2, 1, 2, 2, 2))
  expect_equal(segmented_data_qc(a, compare_seg_data = b)$compare_agreement, 1)
  qc = segmented_data_qc(a, compare_seg_data = c)
  expect_equal(qc$compare_bins, 50)
  # segments 2 and 3 (20 of 50 bins) disagree
  expect_equal(qc$compare_agreement, 0.6)
  expect_lt(qc$compare_correlation, 0)
  # a doubled profile, and one with changes compressed by low purity,
  # correlate perfectly
  expect_equal(segmented_data_qc(a, compare_seg_data = b)$compare_correlation, 1)
  impure = transform(a, CN = 2 + (CN - 2) * 0.3)
  expect_equal(segmented_data_qc(a, compare_seg_data = impure)$compare_correlation, 1)
  expect_equal(segmented_data_qc(a[, -5], value_scale = "log_ratio")$baseline_mode, 2)
  expect_error(segmented_data_qc(a, value_scale = "x"), "value_scale")
})
