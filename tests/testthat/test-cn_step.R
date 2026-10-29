step_seg = function(id, cn18 = c(2, 2), cn14 = c(2, 2), module = "battenberg"){
  data.frame(ID = id, chrom = c("18", "18", "14", "14"),
             start = c(5.5e7, 60.8e6, 1e8, 106.33e6), end = c(60.8e6, 7e7, 106.33e6, 107349540),
             CN = c(cn18, cn14), module = module)
}

test_that("cn_step measures each side and the step", {
  seg = rbind(step_seg("up", cn18 = c(3, 2)), step_seg("flat"))
  st = cn_step(seg, chrom = "chr18", pos = 60.8e6)
  expect_equal(st$cn_left, c(3, 2)); expect_equal(st$cn_right, c(2, 2))
  expect_equal(st$step, c(-1, 0))
  # sparse sides are NA; metadata samples without segments are included
  st2 = cn_step(seg, chrom = "18", pos = 60.8e6, flank = 2e7, these_samples_metadata =
                  data.frame(sample_id = c("up", "none")))
  expect_true(is.na(st2$cn_left[1]))
  expect_true(all(is.na(st2[st2$sample_id == "none", c("cn_left", "cn_right")])))
})

test_that("controlfreec gaps count as neutral", {
  cf = step_seg("cf", cn18 = c(3, 3), module = "controlfreec")[1, ]  # only the left side reported
  st = cn_step(cf, chrom = "18", pos = 60.8e6)
  expect_equal(c(st$cn_left, st$cn_right), c(3, 2))
})

test_that("cn_step_concordance counts opposite steps by group", {
  seg = rbind(step_seg("d18", cn18 = c(3, 2), cn14 = c(2, 3)),   # der(18) gain: -, +
              step_seg("d14", cn18 = c(2, 3), cn14 = c(3, 2)),   # der(14) gain: +, -
              step_seg("same", cn18 = c(3, 2), cn14 = c(3, 2)),
              step_seg("flat"))
  meta = data.frame(sample_id = c("d18", "d14", "same", "flat"), ba = c("POS", "POS", "NEG", "NEG"))
  conc = cn_step_concordance(seg, a = list(chrom = "18", pos = 60.8e6),
                             b = list(chrom = "14", pos = 106.33e6),
                             these_samples_metadata = meta, group_by = "ba")
  expect_equal(conc$steps$pattern, c("expected", "expected", "other", "none"))
  expect_equal(as.character(conc$steps$step_type), c("a down, b up", "a up, b down", "a down, b down", "none"))
  named = cn_step_concordance(seg, a = list(name = "BCL2", chrom = "18", pos = 60.8e6),
                              b = list(name = "IGH", chrom = "14", pos = 106.33e6),
                              these_samples_metadata = meta)
  expect_equal(levels(named$steps$step_type)[1:2], c("BCL2 down, IGH up", "BCL2 up, IGH down"))
  pos = conc$summary[conc$summary$group == "POS", ]
  expect_equal(c(pos$n_expected, pos$n_other), c(2, 0))
  expect_equal(dim(conc$fisher$table), c(2, 2))
})

test_that("cn_transition finds each sample's largest transition and its position", {
  seg = data.frame(ID = c("s1", "s1", "s1", "s2"), chrom = "18",
                   start = c(5e7, 60.9e6, 61.3e6, 5e7), end = c(60.9e6, 61.3e6, 7e7, 7e7),
                   CN = c(3, 2, 1, 2))
  tr = cn_transition(seg, chrom = "18", from = 60e6, to = 62e6, flank = 3e5)
  expect_equal(tr$pos, c(60.9e6, NA))   # s2 has no boundary in the range
  expect_equal(tr$step[2], 0)            # but is covered: measured, no transition
  expect_equal(tr$step[1], -1)
  # in cn_step_concordance, a range gives positions
  seg14 = data.frame(ID = "s1", chrom = "14", start = c(1e8, 106.5e6), end = c(106.5e6, 107349540),
                     CN = c(2, 3))
  conc = cn_step_concordance(rbind(seg, seg14), a = list(chrom = "18", from = 60e6, to = 62e6, flank = 3e5),
                             b = list(chrom = "14", from = 106e6, to = 107e6, flank = 3e5),
                             these_samples_metadata = data.frame(sample_id = c("s1", "s2")))
  expect_equal(conc$steps$a_pos[1], 60.9e6); expect_equal(conc$steps$b_pos[1], 106.5e6)
  expect_equal(conc$steps$pattern[1], "expected")
})

test_that("several partners: each sample gets its best partner, p-arm partners expect the same sign", {
  # MYC (8q) with IGH (14q) or IGK (2p)
  myc = function(id, cn) data.frame(ID = id, chrom = "8", start = c(1.2e8, 128.75e6),
                                    end = c(128.75e6, 1.4e8), CN = cn)
  igh = function(id, cn) data.frame(ID = id, chrom = "14", start = c(1e8, 106.33e6),
                                    end = c(106.33e6, 107349540), CN = cn)
  igk = function(id, cn) data.frame(ID = id, chrom = "2", start = c(8e7, 89.2e6),
                                    end = c(89.2e6, 9e7), CN = cn)
  seg = rbind(myc("t814", c(3, 2)), igh("t814", c(2, 3)), igk("t814", c(2, 2)),   # der(8) gain
              myc("t28", c(3, 2)), igh("t28", c(2, 2)), igk("t28", c(2, 3)))      # IGK: same sign
  # t(2;8) der(8) gain: chr8 centromeric side up (step down), and 2p
  # telomeric side (lower coordinates) up, also a step down
  seg$CN[seg$ID == "t28" & seg$chrom == "2"] = c(3, 2)
  conc = cn_step_concordance(seg, a = list(name = "MYC", chrom = "8", from = 128e6, to = 129e6, flank = 3e5),
                             b = list(list(name = "IGH", chrom = "14", from = 106e6, to = 107e6, flank = 3e5),
                                      list(name = "IGK", chrom = "2", from = 88.8e6, to = 89.6e6, flank = 3e5)),
                             these_samples_metadata = data.frame(sample_id = c("t814", "t28")))
  expect_equal(unname(conc$expect), c("opposite", "same"))
  expect_equal(conc$steps$b_partner, c("IGH", "IGK"))
  expect_equal(conc$steps$pattern, c("expected", "expected"))
  expect_equal(as.character(conc$steps$step_type), c("MYC down, IGH up", "MYC down, IGK down"))
  expect_equal(nrow(conc$partner_summary), 2)
})

test_that("samples a partner can't measure keep their row", {
  seg = data.frame(ID = c("s1", "s1", "s2", "s2"), chrom = "8", start = c(1.2e8, 128.75e6, 1.2e8, 128.75e6),
                   end = c(128.75e6, 1.4e8, 128.75e6, 1.4e8), CN = c(3, 2, 2, 2))
  conc = cn_step_concordance(seg, a = list(chrom = "8", from = 128e6, to = 129e6, flank = 3e5),
                             b = list(list(name = "IGH", chrom = "14", from = 106e6, to = 107e6),
                                      list(name = "IGL", chrom = "22", from = 22.3e6, to = 23.4e6)),
                             these_samples_metadata = data.frame(sample_id = c("s1", "s2", "s3")))
  expect_equal(nrow(conc$steps), 3)
  expect_true(all(is.na(conc$steps$pattern)))
})
