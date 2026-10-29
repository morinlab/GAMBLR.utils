new_dir = function(){
  d = tempfile("bb_")
  dir.create(d)
  d
}
cnloh = function(segments, ...){
  d = new_dir()
  write_battenberg_sample(d, "S1", segments,
                          cellularity = data.frame(cellularity = 1, ploidy = 2, psi = 2))
  summarize_one(d, ...)
}

test_that("events are classified by extent", {
  seg = arm_genome()
  seg = set_region(seg, "3", 0, 3e8, 2, 0)          # whole chromosome
  seg = set_region(seg, "7", 5e7, 3e8, 2, 0)        # all of 7q
  seg = set_region(seg, "1", 0, 3e7, 2, 0)          # terminal 1p
  seg = set_region(seg, "5", 1e8, 1.2e8, 2, 0)      # interstitial 5q
  bb = cnloh(seg)
  e = bb$complete_loh_events
  extent = setNames(e$extent, e$chrom)
  expect_equal(unname(extent[c("3", "7", "1", "5")]),
               c("whole_chromosome", "arm", "terminal", "interstitial"))
  s = bb$summary
  expect_equal(s$cnloh_events, 4)
  expect_equal(s$cnloh_whole_chromosome_events, 1)
  expect_equal(s$cnloh_arm_events, 1)
  expect_equal(s$cnloh_terminal_events, 1)
  expect_equal(s$cnloh_interstitial_events, 1)
  expect_equal(s$cnloh_interstitial_mb, 20)
  expect_gt(s$cnloh_fraction, s$cnloh_segmental_fraction)
})

test_that("a whole-chromosome event is merged across the centromere", {
  seg = set_region(arm_genome(), "1", 0, 3e8, 2, 0)
  e = cnloh(seg)$complete_loh_events
  expect_equal(nrow(e), 1)
  expect_equal(e$extent, "whole_chromosome")
})

test_that("cnloh_reference chooses which LOH counts as copy-neutral in a doubled genome", {
  seg = arm_genome(2, 2)
  seg = set_region(seg, "3", 0, 3e8, 4, 0)          # 4+0: neutral vs the sample
  seg = set_region(seg, "5", 1e8, 1.2e8, 2, 0)      # 2+0: neutral vs diploid
  by_sample = cnloh(seg)
  by_diploid = cnloh(seg, cnloh_reference = "diploid")
  e = by_sample$complete_loh_events
  expect_equal(nrow(e), 2)
  expect_equal(e$loh_type[order(e$chrom)], c("neutral", "loss"))
  expect_equal(e$neutral_vs_sample[e$chrom == "3"], TRUE)
  expect_equal(e$neutral_vs_diploid[e$chrom == "5"], TRUE)
  expect_equal(e$counted[order(e$chrom)], c(TRUE, FALSE))
  expect_equal(by_sample$summary$cnloh_whole_chromosome_events, 1)
  expect_equal(by_sample$summary$cnloh_interstitial_events, 0)
  expect_equal(by_diploid$summary$cnloh_whole_chromosome_events, 0)
  expect_equal(by_diploid$summary$cnloh_interstitial_events, 1)
  expect_equal(by_diploid$summary$cnloh_reference, "diploid")
  expect_equal(by_diploid$complete_loh_events$counted[order(e$chrom)], c(FALSE, TRUE))
})

test_that("small events are listed but not counted", {
  seg = set_region(arm_genome(), "5", 1e8, 1.0005e8, 2, 0)
  bb = cnloh(seg)
  expect_equal(nrow(bb$complete_loh_events), 1)
  expect_false(bb$complete_loh_events$counted)
  expect_equal(bb$summary$cnloh_events, 0)
  expect_gt(bb$summary$cnloh_fraction, 0)
})

test_that("unresolved subclonal LOH is left out", {
  seg = set_region(arm_genome(), "5", 1e8, 1.2e8, 2, 0,
                   frac1_A = 0.6, nMaj2_A = 1, nMin2_A = 1, frac2_A = 0.4)
  bb = cnloh(seg)
  expect_null(bb$complete_loh_events)
  expect_equal(bb$summary$cnloh_events, 0)
  expect_equal(bb$summary$cnloh_fraction, 0)
})

test_that("copy-neutral LOH with no events gives zero counts", {
  bb = cnloh(set_region(arm_genome(), "5", 1e8, 1.2e8, 1, 0))
  expect_equal(bb$complete_loh_events$loh_type, "loss")
  expect_equal(bb$summary$cnloh_events, 0)
  expect_equal(bb$summary$cnloh_whole_chromosome_mb, 0)
  expect_equal(bb$summary$cnloh_largest_segmental_mb, 0)
})

test_that("the segmental copy-neutral shares leave out whole-chromosome LOH", {
  seg = arm_genome()
  seg = set_region(seg, "3", 0, 3e8, 2, 0)          # whole-chromosome copy-neutral
  seg = set_region(seg, "4", 0, 3e8, 1, 0)          # monosomy
  seg = set_region(seg, "1", 0, 3e7, 2, 0)          # 30 Mb copy-neutral, terminal
  seg = set_region(seg, "5", 1e8, 1.1e8, 1, 0)      # 10 Mb deletion
  seg = set_region(seg, "6", 1e8, 1.1e8, 1, 0)      # 10 Mb deletion
  s = cnloh(seg)$summary
  expect_equal(s$loh_copy_neutral_share_segmental, 29.99 / 49.99)
  expect_equal(s$loh_copy_neutral_event_share, 1 / 3)
  # the original share includes the whole chromosomes
  expect_false(isTRUE(all.equal(s$loh_copy_neutral_share, 30 / 50)))
})

test_that("the shares are NA without segmental LOH", {
  s = cnloh(set_region(arm_genome(), "3", 0, 3e8, 2, 0))$summary
  expect_true(is.na(s$loh_copy_neutral_share_segmental))
  expect_true(is.na(s$loh_copy_neutral_event_share))
})

test_that("LOH is called at features by type", {
  seg = arm_genome()
  seg = set_region(seg, "6", 0, 3e8, 2, 0)          # whole-chromosome copy-neutral
  seg = set_region(seg, "1", 0, 3e7, 2, 0)          # terminal 1p copy-neutral
  seg = set_region(seg, "17", 0, 2.4e7, 1, 0)       # 17p deletion LOH
  bb = cnloh(seg)
  cn = battenberg_to_loh_status(bb, encoding = "binary")
  expect_equal(unname(unlist(cn["S1", c("6p_cnloh", "6q_cnloh", "17p_cnloh")])), c(1, 1, 0))
  expect_true(all(grepl("_cnloh$", colnames(cn))))
  # 1p is mostly not LOH; a focal feature in the region is
  expect_equal(unname(unlist(cn["S1", "1p_cnloh"])), 0)
  del = battenberg_to_loh_status(bb, loh_type = "deletion", encoding = "binary")
  expect_equal(unname(unlist(del["S1", c("17p_deletion_loh", "6q_deletion_loh")])), c(1, 0))
  seg_only = battenberg_to_loh_status(bb, exclude_whole_chromosome = TRUE,
                                      encoding = "binary")
  expect_equal(unname(unlist(seg_only["S1", c("6p_cnloh", "6q_cnloh")])), c(0, 0))
  any_loh = battenberg_to_loh_status(bb, loh_type = "any", encoding = "binary")
  expect_equal(unname(unlist(any_loh["S1", c("6q_loh", "17p_loh")])), c(1, 1))
})

test_that("LOH status gives long output and rejects fixed arguments", {
  bb = cnloh(set_region(arm_genome(), "6", 0, 3e8, 2, 0))
  long = battenberg_to_loh_status(bb, output_format = "long", encoding = "binary")
  expect_true(all(grepl("_cnloh$", long$name)))
  expect_error(battenberg_to_loh_status(bb, centering = "median"), "set by")
})
