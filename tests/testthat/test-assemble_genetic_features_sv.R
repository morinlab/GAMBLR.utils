# public DLBCL genome data through the GAMBLR.open accessors
open_dlbcl_meta = function(){
  meta = GAMBLR.open::get_gambl_metadata()
  meta[meta$pathology == "DLBCL" & meta$seq_type == "genome", ]
}

sv_test_maf = function(n = 4){
  maf = suppressMessages(GAMBLR.open::get_coding_ssm(these_samples_metadata = open_dlbcl_meta()))
  ids = utils::head(unique(maf$Tumor_Sample_Barcode[maf$Hugo_Symbol == "EZH2"]), n)
  list(maf = maf[maf$Tumor_Sample_Barcode %in% ids, ], ids = ids)
}

test_that("sv_data fills in SV status; only rearranged genes with a partner count", {
  skip_if_not_installed("GAMBLR.open")
  d = sv_test_maf()
  ids = d$ids
  meta = data.frame(sample_id = ids,
                    myc_ba = c("POS", "NEG", NA, NA),
                    bcl2_ba = NA, bcl6_ba = NA)
  sv = data.frame(
    tumour_sample_id = c(ids[2], ids[3], ids[3], ids[4], ids[4]),
    gene = c("BCL2", "MYC", "FOXP1", "BCL6", "MYC"),
    partner = c("IGH", "BCL6", "BCL6", NA, "IGH")
  )
  out = suppressWarnings(assemble_genetic_features(d$maf, genes_coding = "EZH2",
                                                   genes_noncoding = NULL,
                                                   these_samples_metadata = meta,
                                                   sv_data = sv))
  out = out[ids, ]
  # MYC: s1 from FISH, s3 (BCL6 as partner) and s4 (IGH) from sv_data
  expect_equal(out$MYC_SV, c(2, 0, 2, 2))
  expect_equal(out$BCL2_SV, c(0, 2, 0, 0))
  # BCL6 as the partner of MYC or FOXP1 doesn't count; nor does an SV with no partner
  expect_equal(out$BCL6_SV, c(0, 0, 0, 0))
})

test_that("without sv_data, SV status comes from metadata as before", {
  skip_if_not_installed("GAMBLR.open")
  d = sv_test_maf()
  meta = data.frame(sample_id = d$ids, myc_ba = c("POS", NA, NA, NA),
                    bcl2_ba = NA, bcl6_ba = NA)
  out = suppressWarnings(assemble_genetic_features(d$maf, genes_coding = "EZH2",
                                                   genes_noncoding = NULL,
                                                   these_samples_metadata = meta))
  expect_equal(out[d$ids, "MYC_SV"], c(2, 0, 0, 0))
  # a missing metadata field still stops when there's no sv_data
  meta$bcl6_ba = NULL
  expect_error(suppressWarnings(assemble_genetic_features(d$maf, genes_coding = "EZH2",
                                                          genes_noncoding = NULL,
                                                          these_samples_metadata = meta)),
               "bcl6_ba")
})

test_that("a missing metadata field is a warning when sv_data is given", {
  skip_if_not_installed("GAMBLR.open")
  d = sv_test_maf()
  meta = data.frame(sample_id = d$ids, myc_ba = NA, bcl2_ba = NA)
  sv = data.frame(tumour_sample_id = d$ids[1], gene = "BCL6", partner = "IGH")
  expect_warning(out <- assemble_genetic_features(d$maf, genes_coding = "EZH2",
                                                  genes_noncoding = NULL,
                                                  these_samples_metadata = meta,
                                                  sv_data = sv),
                 "no bcl6_ba column")
  expect_equal(out[d$ids, "BCL6_SV"], c(2, 0, 0, 0))
})

test_that("assemble_classifier_features ignores genes named only as the partner", {
  meta = data.frame(sample_id = c("s1", "s2", "s3"), bcl2_ba = NA, bcl6_ba = NA, myc_ba = NA)
  sv = data.frame(tumour_sample_id = c("s1", "s2", "s3"),
                  gene = c("MYC", "BCL6", "BCL6"),
                  partner = c("BCL6", NA, "IGH"))
  ft = assemble_classifier_features(meta, classifier = "DLBClass", sv_data = sv)
  bcl6 = ft[ft$feature == "BCL6_SV", ]
  expect_equal(bcl6$value[match(c("s1", "s2", "s3"), bcl6$sample_id)], c(0, 0, 1))
})

test_that("noncoding_window counts only non-coding mutations in each gene's TSS window", {
  skip_if_not_installed("GAMBLR.open")
  skip_if_not("get_ssm_by_genes" %in% getNamespaceExports("GAMBLR.open"),
              "GAMBLR.open has no get_ssm_by_genes yet")
  # coding mutations plus every mutation across BCL2 and PIM1
  meta = utils::head(open_dlbcl_meta(), 40)
  coding = suppressMessages(GAMBLR.open::get_coding_ssm(these_samples_metadata = meta,
                                                        include_silent = TRUE))
  regional = suppressMessages(GAMBLR.open::get_ssm_by_genes(c("BCL2", "PIM1"),
                                                            these_samples_metadata = meta))
  maf = rbind(as.data.frame(coding), as.data.frame(regional)[, colnames(coding)])
  maf = maf[!duplicated(maf[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position",
                                "Tumor_Seq_Allele2")]), ]
  ids = utils::head(unique(maf$Tumor_Sample_Barcode[maf$Hugo_Symbol %in% c("BCL2", "PIM1") &
                                                       maf$Variant_Classification == "Intron"]), 15)
  maf = maf[maf$Tumor_Sample_Barcode %in% ids, ]
  meta = data.frame(sample_id = ids)
  run = function(window){
    out = assemble_genetic_features(maf, genes_coding = c("EZH2", "BCL2", "PIM1"),
                                    genes_noncoding = c("BCL2", "PIM1"),
                                    noncoding_window = window,
                                    these_samples_metadata = meta, sv_from_metadata = NULL,
                                    return_everything = TRUE)
    out$full_status = out$full_status[ids, ]
    out$collapsed_status = out$collapsed_status[ids, ]
    out
  }
  whole = run(NULL)
  windowed = run(c(0, 4000))
  # coding indicators don't change
  coding = setdiff(colnames(whole$full_status), c("BCL2_noncoding", "PIM1_noncoding"))
  expect_equal(windowed$full_status[, coding], whole$full_status[, coding])
  # non-coding columns match subset_ssm_by_tss directly
  nc = maf[!maf$Variant_Classification %in% GAMBLR.helpers::vc_nonSynonymous, ]
  in_window = subset_ssm_by_tss(nc, genes = c("BCL2", "PIM1"), start = 0, end = 4000)
  # genes with no mutations in these samples get no column: count them as 0
  column = function(df, name) if(name %in% colnames(df)) unname(df[[name]]) else rep(0, nrow(df))
  for(g in c("BCL2", "PIM1")){
    expected = ifelse(ids %in% in_window$Tumor_Sample_Barcode[in_window$tss_gene == g], 1, 0)
    expect_equal(column(windowed$full_status, paste0(g, "_noncoding")), expected)
    # the gene's column is 2 with a coding mutation, else 1 with a windowed one
    coding_hit = column(windowed$full_status, paste0(g, "_coding")) > 0
    expect_equal(column(windowed$collapsed_status, g), ifelse(coding_hit, 2, expected))
  }
  # the whole-gene count is at least as large here
  expect_true(all(whole$full_status$BCL2_noncoding >= windowed$full_status$BCL2_noncoding))
  expect_true(any(whole$full_status$BCL2_noncoding > windowed$full_status$BCL2_noncoding))
  expect_error(run(4000), "c\\(start, end\\)")
})

test_that("one get_ssm_by_genes MAF, or a coding MAF plus noncoding_maf, give the same result", {
  skip_if_not_installed("GAMBLR.open")
  skip_if_not("get_ssm_by_genes" %in% getNamespaceExports("GAMBLR.open"),
              "GAMBLR.open has no get_ssm_by_genes yet")
  meta = utils::head(open_dlbcl_meta(), 40)
  genes = c("BCL2", "PIM1", "EZH2")
  gene_maf = suppressMessages(GAMBLR.open::get_ssm_by_genes(genes, these_samples_metadata = meta))
  coding_maf = suppressMessages(GAMBLR.open::get_coding_ssm(these_samples_metadata = meta,
                                                            include_silent = TRUE))
  coding_maf = coding_maf[coding_maf$Hugo_Symbol %in% genes, ]
  run = function(maf, noncoding_maf = NULL, window = c(0, 4000)){
    suppressWarnings(assemble_genetic_features(maf, genes_coding = genes,
                                               genes_noncoding = c("BCL2", "PIM1"),
                                               noncoding_maf = noncoding_maf,
                                               noncoding_window = window,
                                               these_samples_metadata = meta,
                                               sv_from_metadata = NULL))
  }
  one_maf = run(gene_maf)
  two_mafs = run(coding_maf, noncoding_maf = gene_maf)
  expect_equal(two_mafs[rownames(one_maf), colnames(one_maf)], one_maf)
  # the coding MAF alone only has silent mutations to offer as non-coding
  coding_only = run(coding_maf)
  expect_true(sum(coding_only$BCL2 == 1) < sum(one_maf$BCL2 == 1))
  # without a window too
  expect_equal(run(coding_maf, gene_maf, NULL)[rownames(one_maf), ],
               run(gene_maf, window = NULL)[rownames(one_maf), ])
})

test_that("seg_data adds copy number columns; unassessed_value marks missing data", {
  skip_if_not_installed("GAMBLR.open")
  meta = utils::head(open_dlbcl_meta(), 12)
  maf = suppressMessages(GAMBLR.open::get_coding_ssm(these_samples_metadata = meta))
  segs = suppressMessages(GAMBLR.open::get_cn_segments(these_samples_metadata = meta))
  no_cn = meta$sample_id[1]
  segs = segs[segs$ID != no_cn, ]
  meta$bcl2_ba[2] = NA
  run = function(...) suppressWarnings(suppressMessages(assemble_genetic_features(
    maf, genes_coding = c("EZH2", "CREBBP", "KMT2D"), genes_noncoding = NULL,
    these_samples_metadata = meta, seg_data = segs, ...)))
  zero = run()
  mut_only = suppressWarnings(suppressMessages(assemble_genetic_features(
    maf, genes_coding = c("EZH2", "CREBBP", "KMT2D"), genes_noncoding = NULL,
    these_samples_metadata = meta)))
  cnv_cols = setdiff(colnames(zero), colnames(mut_only))
  expect_equal(length(cnv_cols), nrow(get_cnv_features()))
  expect_equal(rownames(zero), rownames(mut_only))
  expect_equal(zero[, colnames(mut_only)], mut_only)
  expect_false(anyNA(zero))
  expect_true(all(zero[no_cn, cnv_cols] == 0))

  na = run(unassessed_value = NA)
  expect_true(all(is.na(na[no_cn, cnv_cols])))
  expect_false(anyNA(na[no_cn, setdiff(colnames(mut_only), c("BCL2_SV", "BCL6_SV", "MYC_SV"))]))
  expect_true(is.na(na[meta$sample_id[2], "BCL2_SV"]))
  # assessed values are unchanged
  expect_equal(na[!is.na(na)], zero[!is.na(na)])
  # mutations are only marked with an ssm_available column
  meta$ssm_available = meta$sample_id != meta$sample_id[3]
  na2 = suppressWarnings(suppressMessages(assemble_genetic_features(
    maf, genes_coding = c("EZH2", "CREBBP", "KMT2D"), genes_noncoding = NULL,
    these_samples_metadata = meta, seg_data = segs, unassessed_value = NA)))
  expect_true(all(is.na(na2[meta$sample_id[3], intersect(c("EZH2", "CREBBP", "KMT2D"), colnames(na2))])))
})
