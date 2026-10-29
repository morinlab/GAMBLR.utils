tss_maf = function(pos, vc = "Intron", chrom = "1"){
  data.frame(Tumor_Sample_Barcode = paste0("s", seq_along(pos)), Hugo_Symbol = "X",
             Chromosome = chrom, Start_Position = pos, Variant_Classification = vc)
}
tss_table = data.frame(gene = c("PLUS", "MINUS"), chrom = c("1", "2"),
                       strand = c("+", "-"), tss = c(10000, 50000))

test_that("windows follow the direction of transcription, ends inclusive", {
  # + strand: TSS at 10000, first 4 kb is 10000-14000
  plus = tss_maf(c(9999, 10000, 14000, 14001))
  out = subset_ssm_by_tss(plus, tss = tss_table)
  expect_equal(out$Start_Position, c(10000, 14000))
  expect_equal(out$tss_distance, c(0, 4000))
  expect_equal(out$tss_gene, c("PLUS", "PLUS"))
  # - strand: TSS at 50000, first 4 kb is 46000-50000
  minus = tss_maf(c(45999, 46000, 50000, 50001), chrom = "2")
  out = subset_ssm_by_tss(minus, tss = tss_table)
  expect_equal(out$Start_Position, c(46000, 50000))
  expect_equal(out$tss_distance, c(4000, 0))
})

test_that("negative start reaches upstream of the TSS", {
  plus = tss_maf(c(7999, 8000, 9000))
  out = subset_ssm_by_tss(plus, start = -2000, end = 0, tss = tss_table)
  expect_equal(out$tss_distance, c(-2000, -1000))
  minus = tss_maf(c(51000, 52000, 52001), chrom = "2")
  out = subset_ssm_by_tss(minus, start = -2000, end = 0, tss = tss_table)
  expect_equal(out$tss_distance, c(-1000, -2000))
})

test_that("coding mutations are dropped unless noncoding_only = FALSE", {
  maf = tss_maf(c(10100, 10200, 10300), vc = c("Missense_Mutation", "Silent", "5'UTR"))
  expect_equal(subset_ssm_by_tss(maf, tss = tss_table)$Variant_Classification,
               c("Silent", "5'UTR"))
  expect_equal(nrow(subset_ssm_by_tss(maf, tss = tss_table, noncoding_only = FALSE)), 3)
})

test_that("mutations are assigned by position, once per overlapping window", {
  # two genes whose windows overlap (bidirectional promoter)
  tss = data.frame(gene = c("A", "B"), chrom = "1", strand = c("+", "-"), tss = c(10000, 9000))
  maf = tss_maf(c(8500, 9500, 12000, 20000), chrom = "chr1")
  # A (+, TSS 10000) covers 9000-14000; B (-, TSS 9000) covers 5000-10000
  out = subset_ssm_by_tss(maf, start = -1000, end = 4000, tss = tss)
  out = out[order(out$Start_Position, out$tss_gene), ]
  expect_equal(out$Start_Position, c(8500, 9500, 9500, 12000))
  expect_equal(out$tss_gene, c("B", "A", "B", "A"))
  expect_equal(out$tss_distance, c(500, -500, -500, 2000))
  # restricting to genes
  expect_equal(unique(subset_ssm_by_tss(maf, genes = "B", start = -1000, end = 4000,
                                        tss = tss)$tss_gene), "B")
})

test_that("bad windows and empty results", {
  expect_error(subset_ssm_by_tss(tss_maf(1), start = 10, end = 0, tss = tss_table), "start <= end")
  out = subset_ssm_by_tss(tss_maf(1), tss = tss_table)
  expect_equal(nrow(out), 0)
  expect_true(all(c("tss_gene", "tss_distance") %in% colnames(out)))
})

test_that("get_gene_tss gives one TSS per gene, within the gene", {
  tss = get_gene_tss(c("BCL2", "MYC"))
  expect_equal(tss$gene, c("BCL2", "MYC"))
  expect_equal(tss$strand, c("-", "+"))
  expect_true(all(tss$tss >= tss$start & tss$tss <= tss$end))
  # old and current symbols find the same gene, by alias or Entrez ID
  expect_equal(get_gene_tss("H1-3")$tss, get_gene_tss("HIST1H1D")$tss)
  expect_equal(get_gene_tss("TSPOAP1")$tss, get_gene_tss("BZRAP1")$tss)
  expect_equal(get_gene_tss("TSPOAP1")$gene, "TSPOAP1")
})

test_that("aSHM genes use their most mutated promoter", {
  # ETS1's longest isoform starts 65 kb upstream of the hypermutated promoter
  ets1 = get_gene_tss("ETS1", which = "all")
  expect_gt(nrow(ets1), 1)
  expect_equal(get_gene_tss("ETS1")$tss, 128392205)
  expect_equal(ets1$tss[ets1$primary], 128392205)
  expect_equal(ets1$tss_source[ets1$primary], "ashm_mutations")
  # PAX5 has two active promoters; "active" always includes the primary
  pax5 = get_gene_tss("PAX5", which = "active")
  expect_equal(nrow(pax5), 2)
  expect_equal(sum(pax5$primary), 1)
  # which doesn't apply to other genes
  expect_equal(nrow(get_gene_tss("MYD88", which = "all")), 1)
  expect_true(is.na(get_gene_tss("MYD88")$active))
})

test_that("a mutation in several windows of one gene counts once, for the nearest TSS", {
  tss = data.frame(gene = "A", chrom = "1", strand = "+", tss = c(10000, 11000))
  out = subset_ssm_by_tss(tss_maf(c(10500, 11500, 14500)), tss = tss)
  expect_equal(out$Start_Position, c(10500, 11500, 14500))
  expect_equal(out$tss_position, c(10000, 11000, 11000))
  expect_equal(out$tss_distance, c(500, 500, 3500))
})

test_that("a maf_data input comes back as maf_data", {
  maf = create_maf_data(tss_maf(c(10000, 14000, 20000)), "grch37")
  out = subset_ssm_by_tss(maf, tss = tss_table)
  expect_s3_class(out, "maf_data")
  expect_equal(get_genome_build(out), "grch37")
  expect_equal(nrow(out), 2)
  # also when nothing is in a window, and from the attribute rather than the argument
  out = subset_ssm_by_tss(create_maf_data(tss_maf(1), "hg38"), tss = tss_table)
  expect_s3_class(out, "maf_data")
  expect_equal(get_genome_build(out), "hg38")
  # plain data frames stay plain
  expect_false(inherits(subset_ssm_by_tss(tss_maf(10000), tss = tss_table), "maf_data"))
})

test_that("get_gene_tss returns windows as regions or bed", {
  df = get_gene_tss(c("BCL2", "MYC"), window = c(0, 4000))
  # BCL2 is on the - strand, MYC on the +
  expect_equal(df$window_start, c(df$tss[1] - 4000, df$tss[2]))
  expect_equal(df$window_end, c(df$tss[1], df$tss[2] + 4000))
  region = get_gene_tss(c("BCL2", "MYC"), window = c(0, 4000), return_as = "region")
  expect_equal(names(region), c("BCL2", "MYC"))
  expect_equal(unname(region), paste0(df$chrom, ":", df$window_start, "-", df$window_end))
  bed = get_gene_tss(c("BCL2", "MYC"), window = c(0, 4000), return_as = "bed")
  expect_equal(colnames(bed), c("chromosome", "start", "end", "hugo_symbol"))
  # hg38 regions have a chr prefix, as from gene_to_region
  expect_match(get_gene_tss("MYC", "hg38", return_as = "region"), "^chr8:")
  # the windows are those subset_ssm_by_tss uses
  maf = tss_maf(c(df$window_start[1], df$window_end[1]), chrom = df$chrom[1])
  expect_equal(nrow(subset_ssm_by_tss(maf, genes = "BCL2")), 2)
})
