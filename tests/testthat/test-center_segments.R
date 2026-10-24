tetraploid = function(){
  data.frame(ID = "s1", chrom = as.character(1:22), start = 1, end = 1e8,
             CN = c(3, 6, rep(4, 20)))
}

test_that("centring puts the modal state at CN 2 and leaves log.ratio alone by default", {
  seg = tetraploid()
  seg$log.ratio = c(-0.3, 0.4, rep(-0.05, 20))   # caller's own scale: neutral near 0
  centred = center_segments(seg)
  expect_equal(centred$CN, c(1, 4, rep(2, 20)))
  expect_equal(centred$log.ratio, seg$log.ratio)
  expect_equal(centred$CN_uncentered, seg$CN)
  expect_equal(unique(centred$centre_cn), 4)
  # optionally centred on its own neutral segments, not by log2(4 / 2)
  lr = center_segments(seg, center_log_ratio = TRUE)
  expect_equal(lr$log.ratio[3:22], rep(0, 20))
  expect_equal(lr$log.ratio[1], -0.25)
  expect_error(center_segments(seg[, c("ID", "chrom", "start", "end", "log.ratio")]), "needs columns")
})

test_that("fractional offsets are removed by mode but not by mode_rounded", {
  seg = tetraploid()
  seg$CN = seg$CN - 1.7  # neutral at 2.3
  expect_equal(center_segments(seg)$CN[3], 2)
  expect_equal(center_segments(seg, method = "mode_rounded")$CN[3], 2.3)
})

test_that("an existing baseline_cn column is used, and centring twice stops", {
  seg = tetraploid()
  seg$baseline_cn = 3; seg$baseline_method = "mode"
  expect_equal(center_segments(seg)$CN[3], 3)
  expect_error(center_segments(center_segments(seg)), "already centred")
})

test_that("seg_data_to_igv writes an IGV seg table", {
  f = tempfile(fileext = ".seg")
  out = seg_data_to_igv(center_segments(tetraploid()), f, chr_prefix = TRUE)
  back = utils::read.delim(f)
  expect_equal(colnames(back)[1:5], c("ID", "chrom", "loc.start", "loc.end", "num.mark"))
  expect_equal(nrow(back), 22)
  expect_true(all(startsWith(back$chrom, "chr")))
  expect_equal(colnames(back)[6], "log2_CN")         # the default with a CN column
  expect_equal(back$log2_CN[back$chrom == "chr3"], 0)  # centred CN 2
  expect_equal(back$log2_CN[back$chrom == "chr1"], -1) # centred CN 1
  expect_equal(colnames(seg_data_to_igv(tetraploid(), value = "CN"))[6], "CN")
  lr = transform(tetraploid(), log.ratio = 0.1)
  expect_equal(seg_data_to_igv(lr, value = "log.ratio")$seg.mean, rep(0.1, 22))
  expect_error(seg_data_to_igv(tetraploid(), value = "log.ratio"), "no log.ratio")
})

test_that("seg_data_to_igv writes an IGV sample information file", {
  seg = rbind(tetraploid(), transform(tetraploid(), ID = "s2"))
  meta = data.frame(sample_id = c("s1", "s2", "s3"), seq_type = "genome", cohort = "c1",
                    pathology = c("BL", "DLBCL", "BL"), lymphgen = c("EZB", NA, "MCD"))
  f = tempfile(fileext = ".seg")
  res = suppressMessages(seg_data_to_igv(seg, f, these_samples_metadata = meta,
                                         additional_columns = "lymphgen"))
  info_file = sub("\\.seg$", ".sample_info.txt", f)
  expect_true(file.exists(info_file))
  lines = readLines(info_file)
  expect_equal(lines[1], "#sampleTable")
  expect_equal(lines[2], "ID\tseq_type\tcohort\tpathology\tlymphgen")
  expect_equal(nrow(res$sample_info), 2)          # s3 has no segments
  # segments of samples not in the metadata are dropped
  res2 = suppressMessages(seg_data_to_igv(seg, these_samples_metadata = meta[meta$sample_id == "s1", ]))
  expect_equal(unique(res2$seg$ID), "s1")
  expect_equal(res2$sample_info$ID, "s1")
  expect_true("#colors" %in% lines)
  expect_true(any(startsWith(lines, "pathology\tBL\t")))
  expect_error(seg_data_to_igv(seg, these_samples_metadata = meta, additional_columns = "nope"),
               "has column")
  # per-sample columns from the segments
  seg$module = ifelse(seg$ID == "s1", "battenberg", "controlfreec")
  seg$module[seg$ID == "s2"][1] = "sequenza"
  res3 = suppressMessages(seg_data_to_igv(seg, these_samples_metadata = meta,
                                          additional_columns = c("lymphgen", "module")))
  expect_equal(colnames(res3$sample_info), c("ID", "seq_type", "cohort", "pathology", "lymphgen", "module"))
  expect_equal(res3$sample_info$module, c("battenberg", "controlfreec+sequenza"))
})
