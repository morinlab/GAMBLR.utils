seg_tiled = function(){
  data.frame(ID = rep(c("s1", "s2"), each = 3), chrom = "1",
             start = c(1, 101, 201, 1, 151, 301),
             end = c(100, 200, 300, 150, 300, 400), CN = c(2, 3, 2, 2, 1, 2))
}

test_that("non-overlapping 1-based segments pass", {
  expect_equal(nrow(check_segment_overlaps(seg_tiled())), 0)
})

test_that("a data frame bound to itself is reported and fixed", {
  seg = rbind(seg_tiled(), seg_tiled())
  ov = check_segment_overlaps(seg)
  expect_equal(nrow(ov), 6)
  expect_true(all(ov$duplicate_row))
  fixed = suppressMessages(resolve_segment_overlaps(seg))
  expect_equal(nrow(fixed), 6)
  expect_equal(nrow(check_segment_overlaps(fixed)), 0)
})

test_that("partial overlaps are trimmed and nested ones dropped", {
  seg = data.frame(ID = "s1", chrom = "1", start = c(1, 90, 120, 500),
                   end = c(100, 200, 150, 600), CN = c(2, 3, 4, 2))
  ov = check_segment_overlaps(seg)
  expect_equal(nrow(ov), 2)
  expect_equal(ov$overlap_bp, c(10, 80))
  fixed = suppressMessages(resolve_segment_overlaps(seg))
  expect_equal(fixed$start, c(1, 101, 500))
  expect_equal(nrow(check_segment_overlaps(fixed)), 0)
})

test_that("same coordinates with different values keep the first", {
  seg = data.frame(ID = "s1", chrom = "1", start = c(1, 1), end = c(100, 100), CN = c(2, 3))
  ov = check_segment_overlaps(seg)
  expect_true(ov$same_coordinates)
  expect_false(ov$duplicate_row)
  fixed = suppressMessages(resolve_segment_overlaps(seg))
  expect_equal(fixed$CN, 2)
})

test_that("profiles of different seq_types under one ID are not overlaps", {
  seg = rbind(cbind(seg_tiled(), seg_seq_type = "genome"),
              cbind(seg_tiled(), seg_seq_type = "capture"))
  expect_equal(nrow(check_segment_overlaps(seg)), 0)
  expect_equal(nrow(check_segment_overlaps(seg, by = "ID")), 6)
})

test_that("tolerance, drop_samples and error work", {
  seg = data.frame(ID = c("s1", "s1", "s2", "s2"), chrom = "1",
                   start = c(1, 99, 1, 50), end = c(100, 200, 100, 200), CN = 2)
  expect_equal(nrow(check_segment_overlaps(seg)), 2)
  expect_equal(nrow(check_segment_overlaps(seg, tolerance = 1)), 1)
  kept = suppressMessages(resolve_segment_overlaps(seg, method = "drop_samples", tolerance = 1))
  expect_equal(unique(kept$ID), "s1")
  expect_error(resolve_segment_overlaps(seg, method = "error", verbose = FALSE), "overlapping")
})

test_that("seg_data attributes are kept", {
  seg = create_seg_data(rbind(seg_tiled(), seg_tiled()), "grch37")
  fixed = suppressMessages(resolve_segment_overlaps(seg))
  expect_s3_class(fixed, "seg_data")
  expect_equal(attr(fixed, "genome_build"), "grch37")
})

test_that("handle_segment_overlaps warns, errors, resolves or skips", {
  seg = rbind(seg_tiled(), seg_tiled())
  expect_warning(handle_segment_overlaps(seg), "overlapping segment")
  expect_error(handle_segment_overlaps(seg, "error"), "overlapping segment")
  expect_equal(nrow(suppressMessages(handle_segment_overlaps(seg, "resolve"))), 6)
  expect_equal(nrow(handle_segment_overlaps(seg, "none")), 12)
  expect_silent(handle_segment_overlaps(seg_tiled()))
})

test_that("a nested segment's overlap is its length, and trimming drops it", {
  seg = data.frame(ID = "s1", chrom = "1", start = c(766292, 146348227, 149038009),
                   end = c(121345136, 249222482, 149220070), CN = c(2, 3.25, 2))
  ov = check_segment_overlaps(seg)
  expect_equal(ov$overlap_bp, 149220070 - 149038009)
  fixed = suppressMessages(resolve_segment_overlaps(seg))
  expect_equal(fixed$CN, c(2, 3.25))
})
