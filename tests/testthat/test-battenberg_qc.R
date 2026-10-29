new_dir = function(){
  d = tempfile("bb_")
  dir.create(d)
  d
}

test_that("the three fitting stages and the final ploidy are extracted", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments(),
                          rho_and_psi = fit_table(ascat_ploidy = 2.27, clonal_ploidy = 2.76,
                                                  clonal_distance = 0.6, ref_seg_ploidy = 2.91),
                          cellularity = data.frame(cellularity = 0.95, ploidy = 2.1, psi = 1.85))
  s = summarize_one(d)$summary
  expect_equal(s$ascat_ploidy, 2.27)
  expect_equal(s$ascat_purity, 1)
  expect_equal(s$clonal_fit_ploidy, 2.76)
  expect_equal(s$clonal_fit_purity, 0.95)
  expect_equal(s$clonal_fit_score, 0.6)
  expect_equal(s$ref_seg_ploidy, 2.91)
  expect_equal(s$final_profile_ploidy, 2.1)
  expect_equal(s$ploidy, 2.1)
  expect_equal(s$purity, 0.95)
})

test_that("ploidy differences are signed (later stage minus earlier)", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments(),
                          rho_and_psi = fit_table(ascat_ploidy = 2.27, clonal_ploidy = 2.76),
                          cellularity = data.frame(cellularity = 0.95, ploidy = 1.89, psi = 1.85))
  s = summarize_one(d)$summary
  expect_equal(s$ploidy_delta_ascat_clonal, 2.76 - 2.27)
  expect_equal(s$ploidy_delta_clonal_final, 1.89 - 2.76)
  expect_equal(s$ploidy_delta_ascat_final, 1.89 - 2.27)
})

test_that("the clonal fit score is FRAC_GENOME's even when is.best marks another row", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments(),
                          rho_and_psi = fit_table(clonal_distance = 0.6, ref_seg_distance = 0.4,
                                                  best = "REF_SEG"),
                          cellularity = data.frame(cellularity = 0.95, ploidy = 2.1, psi = 1.85))
  s = summarize_one(d)$summary
  expect_equal(s$clonal_fit_score, 0.6)
  # fit_distance keeps its meaning: the is.best row's distance
  expect_equal(s$fit_distance, 0.4)
})

test_that("missing fitting files give NA fields and unknown QC", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments())
  s = summarize_one(d)$summary
  expect_true(is.na(s$clonal_fit_score))
  expect_true(is.na(s$final_profile_ploidy))
  expect_true(s$qc_missing_fit_info)
  expect_equal(s$allelic_qc_status, "unknown")
  expect_match(s$allelic_qc_reasons, "missing fitting information")
  # the copy number summaries are still computed
  expect_equal(s$loh_fraction, 0)
})

test_that("a malformed rho_and_psi file is treated as missing, not guessed", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments(), malformed_rho = TRUE,
                          cellularity = data.frame(cellularity = 0.95, ploidy = 2.1, psi = 1.85))
  s = summarize_one(d)$summary
  expect_true(is.na(s$ascat_ploidy))
  expect_true(is.na(s$clonal_fit_score))
  expect_equal(s$allelic_qc_status, "unknown")
})

test_that("refit suggestions are counted", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments(), rho_and_psi = fit_table(),
                          cellularity = data.frame(cellularity = 0.95, ploidy = 2.1, psi = 1.85),
                          refit_rows = data.frame(chrom = c(17, 20), maj = 2, min = 1))
  expect_equal(summarize_one(d)$summary$refit_suggestions, 2)
})

test_that("several failing criteria are all reported, with severe taking precedence", {
  d = new_dir()
  seg = genome_segments()
  # every segment subclonal (above 0.9 subclonal), most with disagreeing
  # states below min_clonal_fraction (low assessable fraction)
  seg$frac1_A = 0.6; seg$nMaj2_A = 2; seg$nMin2_A = 0; seg$frac2_A = 0.4
  seg$frac1_A[1:3] = 0.95; seg$frac2_A[1:3] = 0.05
  write_battenberg_sample(d, "S1", seg,
                          rho_and_psi = fit_table(clonal_ploidy = 2.76, clonal_distance = 0.02),
                          cellularity = data.frame(cellularity = 0.95, ploidy = 1.9, psi = 1.85))
  s = summarize_one(d)$summary
  expect_true(s$qc_fit_severe)
  expect_true(s$qc_high_subclonal)
  expect_true(s$qc_low_assessable)
  expect_true(s$qc_ploidy_discordant)
  expect_equal(s$allelic_qc_status, "severe_review")
  for(reason in c("clonal fit score < 0.1", "subclonal fraction > 0.9",
                  "assessable fraction < 0.7", "|final - clonal fit ploidy| > 0.5")){
    expect_true(grepl(reason, s$allelic_qc_reasons, fixed = TRUE))
  }
})

test_that("statuses: pass, review and unknown take precedence correctly", {
  d = new_dir()
  write_battenberg_sample(d, "PASS", genome_segments(),
                          rho_and_psi = fit_table(clonal_ploidy = 2.0, clonal_distance = 0.8),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  write_battenberg_sample(d, "REVIEW", genome_segments(),
                          rho_and_psi = fit_table(clonal_ploidy = 2.0, clonal_distance = 0.3),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  s = summarize_one(d)$summary
  expect_equal(s$allelic_qc_status[s$sample == "PASS"], "pass")
  expect_equal(s$allelic_qc_reasons[s$sample == "PASS"], "")
  expect_equal(s$allelic_qc_status[s$sample == "REVIEW"], "review")
  expect_true(s$qc_fit_review[s$sample == "REVIEW"])
})

test_that("an unknown distance metric switches off the fit-score criteria", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments(),
                          rho_and_psi = fit_table(clonal_ploidy = 2.0, clonal_distance = 0.02),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  bb = summarize_one(d)
  expect_equal(bb$summary$allelic_qc_status, "severe_review")
  s = battenberg_qc(bb, clonality_dist_metric = NA)$summary
  expect_false(s$qc_fit_severe)
  expect_equal(s$allelic_qc_status, "pass")
  expect_match(s$allelic_qc_reasons, "fit score not used")
})

test_that("thresholds are configurable and unknown ones are rejected", {
  d = new_dir()
  write_battenberg_sample(d, "S1", genome_segments(),
                          rho_and_psi = fit_table(clonal_ploidy = 2.0, clonal_distance = 0.3),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  bb = summarize_one(d)
  expect_equal(bb$summary$allelic_qc_status, "review")
  relaxed = battenberg_qc(bb, qc_thresholds = list(review_fit_score = 0.2))
  expect_equal(relaxed$summary$allelic_qc_status, "pass")
  expect_error(battenberg_qc(bb, qc_thresholds = list(foo = 1)), "unknown qc_thresholds")
})

test_that("no assessable sequence gives NA, not zero, and excluded LOH is reported", {
  d = new_dir()
  seg = genome_segments()
  # majority state 2+0 (LOH) but the second state disagrees: nothing assessable
  seg$nMaj1_A = 2; seg$nMin1_A = 0; seg$frac1_A = 0.6
  seg$nMaj2_A = 1; seg$nMin2_A = 1; seg$frac2_A = 0.4
  write_battenberg_sample(d, "S1", seg, rho_and_psi = fit_table(),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  s = summarize_one(d)$summary
  expect_equal(s$assessable_fraction, 0)
  expect_true(is.na(s$complete_loh_fraction))
  expect_true(is.na(s$relative_loss_fraction))
  expect_equal(s$excluded_loh_fraction, 1)
  expect_equal(s$loh_fraction, 1)
  expect_true(s$qc_low_assessable)
})

test_that("original LOH summaries are unchanged by the allele-specific additions", {
  d = new_dir()
  seg = genome_segments()
  seg$nMaj1_A[1] = 2; seg$nMin1_A[1] = 0    # copy-neutral LOH on chr1
  seg$nMaj1_A[2] = 1; seg$nMin1_A[2] = 0    # LOH with loss on chr2
  write_battenberg_sample(d, "S1", seg, rho_and_psi = fit_table(),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  bb = summarize_one(d)
  s = bb$summary
  expect_equal(s$modal_cn, 2)
  expect_equal(s$loh_neutral_fraction, 1 / 22)
  expect_equal(s$loh_loss_fraction, 1 / 22)
  expect_equal(s$complete_loh_fraction, 2 / 22)
  expect_equal(s$excluded_loh_fraction, 0)
  expect_setequal(bb$loh_segments$class, c("loh_neutral", "loh_loss"))
  expect_true(all(c("loh_fraction", "loh_neutral_events", "loh_loss_terminal_fraction",
                    "whole_chromosome_loh", "segments_per_gb") %in% colnames(s)))
  states = bb$allelic_states
  expect_equal(states$assessable_fraction[states$allelic_state == "2+0"], 1 / 22)
})

test_that("QC filtering and the report don't modify the original", {
  d = new_dir()
  write_battenberg_sample(d, "PASS", genome_segments(),
                          rho_and_psi = fit_table(clonal_ploidy = 2.0, clonal_distance = 0.8),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  write_battenberg_sample(d, "SEVERE", genome_segments(),
                          rho_and_psi = fit_table(clonal_ploidy = 2.0, clonal_distance = 0.02),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  bb = summarize_one(d)
  kept = battenberg_qc_filter(bb)
  expect_equal(kept$summary$sample, "PASS")
  expect_true(all(kept$segments$sample == "PASS"))
  expect_equal(nrow(bb$summary), 2)
  report = battenberg_qc_report(bb, groups = c(PASS = "C6", SEVERE = "C8"))
  expect_equal(sum(report$status$n), 2)
  expect_equal(report$status$percent_pass[report$status$group == "C6"], 100)
})

test_that("94-15772_tumorB (real data) is flagged for severe review", {
  root = "/Users/rmorin/Dropbox/projects/dlbcl_cn_clustering/battenberg_check"
  f = list.files(root, "^94-15772_tumorB_subclones\\.txt$", recursive = TRUE, full.names = TRUE)
  skip_if(length(f) == 0, "local Battenberg check files not available")
  s = summarize_battenberg(subclones_files = f, loh_genes = NULL)$summary
  expect_equal(s$allelic_qc_status, "severe_review")
  expect_equal(s$clonal_fit_score, 0.0241948525051919, tolerance = 1e-9)
  expect_equal(s$clonal_fit_ploidy, 2.76265406143488, tolerance = 1e-9)
  expect_equal(s$final_profile_ploidy, 1.88779851115333, tolerance = 1e-9)
  expect_equal(s$ascat_ploidy, 2.27313480668601, tolerance = 1e-9)
  expect_true(s$qc_fit_severe && s$qc_high_subclonal && s$qc_ploidy_discordant)
  expect_true(s$excluded_loh_fraction > 0.08)
})

test_that("large homozygous deletions are measured and flagged", {
  d = new_dir()
  fit = list(rho_and_psi = fit_table(clonal_ploidy = 2.0, clonal_distance = 0.8),
             cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  small = genome_segments(); small$nMaj1_A[5] = 0; small$nMin1_A[5] = 0
  small$endpos[5] = small$startpos[5] + 2e6                      # 2 Mb homdel
  write_battenberg_sample(d, "SMALL", small, fit$rho_and_psi, fit$cellularity)
  mid = genome_segments(); mid$nMaj1_A[5] = 0; mid$nMin1_A[5] = 0  # 45 Mb homdel
  write_battenberg_sample(d, "MID", mid, fit$rho_and_psi, fit$cellularity)
  big = genome_segments(); big$nMaj1_A[5:6] = 0; big$nMin1_A[5:6] = 0
  big$endpos[5] = 1.2e8                                            # 115 Mb homdel
  write_battenberg_sample(d, "BIG", big, fit$rho_and_psi, fit$cellularity)
  s = summarize_one(d)$summary
  get = function(id, col) s[[col]][s$sample == id]
  expect_equal(get("SMALL", "homdel_mb"), 2)
  expect_equal(get("SMALL", "largest_homdel_mb"), 2)
  expect_false(get("SMALL", "qc_homdel_review"))
  expect_equal(get("SMALL", "allelic_qc_status"), "pass")
  expect_equal(get("MID", "largest_homdel_mb"), 45)
  expect_true(get("MID", "qc_homdel_review"))
  expect_equal(get("MID", "allelic_qc_status"), "review")
  expect_equal(get("BIG", "largest_homdel_mb"), 115)
  expect_true(get("BIG", "qc_homdel_severe"))
  expect_equal(get("BIG", "allelic_qc_status"), "severe_review")
  expect_match(get("BIG", "allelic_qc_reasons"), "homozygous deletion > 50 Mb")
})

test_that("mean_total_cn ignores segments without an exact total", {
  d = new_dir()
  seg = genome_segments()
  seg$frac1_A = 1; seg$nMaj2_A = NA; seg$nMin2_A = NA; seg$frac2_A = NA
  # chr1: subclonal 1+1 vs 1+0 at 0.6/0.4 -> exact total unknown
  seg$frac1_A[1] = 0.6; seg$nMaj2_A[1] = 1; seg$nMin2_A[1] = 0; seg$frac2_A[1] = 0.4
  seg$nMaj1_A[2] = 2; seg$nMin1_A[2] = 1                       # chr2: 2+1, total 3
  write_battenberg_sample(d, "S1", seg, rho_and_psi = fit_table(),
                          cellularity = data.frame(cellularity = 0.9, ploidy = 2.0, psi = 2))
  s = summarize_one(d)$summary
  expect_false(is.na(s$mean_total_cn))
  expect_equal(s$mean_total_cn, (20 * 2 + 3) / 21)
})
