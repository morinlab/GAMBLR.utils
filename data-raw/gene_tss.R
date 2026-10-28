# Build gene_tss: the transcription start sites used by get_gene_tss() and
# subset_ssm_by_tss().
#
# Every gene in GAMBLR.data's gene coordinates gets one TSS: that of its
# longest isoform (longest spliced transcript in UCSC knownGene: UCSC Genes
# for grch37, GENCODE V44 for hg38), or the gene span's start where knownGene
# has no transcript for it. Genes in the aSHM catalogue
# (GAMBLR.data::somatic_hypermutation_locations_*_v_latest) instead get every
# distinct transcript TSS (from both builds' annotations, lifted over, so
# the two builds have the same promoters), each with the number of GAMBLR.open genomes
# mutated in its first 4 kb, counting only mutations in the gene's curated
# aSHM regions. The most mutated is the gene's primary TSS and those with
# enough mutations are active. Only aSHM genes are counted because
# GAMBLR.open bundles non-coding mutations only for them.
#
# Needs, at build time only: TxDb.Hsapiens.UCSC.hg19.knownGene,
# TxDb.Hsapiens.UCSC.hg38.knownGene, GenomicFeatures and GAMBLR.open.

library(GenomicFeatures)
library(dplyr) # after GenomicFeatures, so dplyr's select() isn't masked

count_window = c(0, 4000)   # the window mutations are counted in
min_samples = 10            # active: at least this many mutated genomes
min_fraction = 0.1          # and at least this fraction of the gene's best TSS
merge_gap = 2000            # TSSs closer than this are merged (half the window)

genomes = GAMBLR.open::get_gambl_metadata() %>% filter(seq_type == "genome")

# the genes, transcripts and primary TSSs of a build, and the TSSs of its
# aSHM genes' transcripts
prepare_build = function(build){
  txdb = if(build == "grch37"){
    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  }else{
    TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  coords = as.data.frame(if(build == "grch37") GAMBLR.data::grch37_gene_coordinates else GAMBLR.data::hg38_gene_coordinates)
  coords$chromosome = sub("^chr", "", coords$chromosome)
  coords = coords[coords$chromosome %in% c(1:22, "X", "Y"), ]
  coords$symbol = ifelse(is.na(coords$hugo_symbol), coords$gene_name, coords$hugo_symbol)
  coords = coords[!is.na(coords$symbol), ]
  # one span per gene, on its first chromosome
  spans = coords %>%
    group_by(symbol) %>%
    filter(chromosome == chromosome[1]) %>%
    summarise(chrom = chromosome[1], start = min(start), end = max(end),
              ensembl_gene_id = ensembl_gene_id[1], .groups = "drop")

  # Entrez IDs (knownGene's gene IDs): symbol, then Ensembl ID, then alias.
  # Symbol first: an Ensembl ID can map to several Entrez IDs (MEF2B's also
  # maps to the MEF2BNB-MEF2B readthrough)
  org = org.Hs.eg.db::org.Hs.eg.db
  eg = suppressMessages(AnnotationDbi::mapIds(org, spans$symbol, "ENTREZID", "SYMBOL",
                                              multiVals = "first"))
  todo = is.na(eg)
  eg[todo] = suppressMessages(AnnotationDbi::mapIds(org, spans$ensembl_gene_id[todo], "ENTREZID",
                                                    "ENSEMBL", multiVals = "first"))
  todo = is.na(eg)
  eg[todo] = suppressMessages(AnnotationDbi::mapIds(org, spans$symbol[todo], "ENTREZID", "ALIAS",
                                                    multiVals = "first"))
  spans$entrez = unname(eg)

  # every transcript with its spliced length, kept if it's on the gene's
  # chromosome and within 100 kb of its span
  gr = transcripts(txdb, columns = "tx_id")
  tx = data.frame(tx_id = mcols(gr)$tx_id, chrom = sub("^chr", "", as.character(seqnames(gr))),
                  tx_start = start(gr), tx_end = end(gr), strand = as.character(strand(gr)))
  lens = transcriptLengths(txdb)
  tx = tx %>%
    inner_join(lens[, c("tx_id", "gene_id", "tx_len")], by = "tx_id") %>%
    filter(!is.na(gene_id)) %>%
    transmute(entrez = as.character(gene_id), chrom, strand, tx_len,
              tss = ifelse(strand == "+", tx_start, tx_end), width = tx_end - tx_start + 1) %>%
    inner_join(spans %>% select(gene = symbol, entrez, chrom, start, end),
               by = c("entrez", "chrom"), relationship = "many-to-many") %>%
    filter(tss >= start - 1e5, tss <= end + 1e5)

  # longest isoform (ties: wider, then more upstream)
  longest = tx %>%
    group_by(gene) %>%
    arrange(desc(tx_len), desc(width), ifelse(strand == "+", tss, -tss), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(gene, strand, tss, tss_source = "longest_isoform")

  # genes without a knownGene transcript: the start of the gene span, with
  # the strand from org.Hs.eg.db
  rest = spans %>% filter(!symbol %in% longest$gene, !is.na(entrez))
  loc = AnnotationDbi::mget(unique(rest$entrez), org.Hs.eg.db::org.Hs.egCHRLOC, ifnotfound = NA)
  sgn = vapply(loc, function(x) if(all(is.na(x))) NA_real_ else sign(x[!is.na(x)][1]), numeric(1))
  rest = rest %>%
    mutate(strand = ifelse(sgn[entrez] > 0, "+", ifelse(sgn[entrez] < 0, "-", NA))) %>%
    filter(!is.na(strand)) %>%
    transmute(gene = symbol, strand, tss = ifelse(strand == "+", start, end), tss_source = "gene_span")

  primary = bind_rows(longest, rest) %>%
    mutate(primary = TRUE)

  # aSHM genes: every distinct TSS, with the genomes mutated in its window
  catalogue = if(build == "grch37") GAMBLR.data::somatic_hypermutation_locations_GRCh37_v_latest else GAMBLR.data::somatic_hypermutation_locations_GRCh38_v_latest
  # catalogue genes matched to the coordinate table's genes by symbol, then
  # by Entrez ID (the two can use different names, e.g. TSPOAP1 and BZRAP1)
  ashm_symbols = unique(catalogue$gene)
  ashm_entrez = suppressMessages(AnnotationDbi::mapIds(org, ashm_symbols, "ENTREZID", "SYMBOL",
                                                       multiVals = "first"))
  todo = is.na(ashm_entrez)
  ashm_entrez[todo] = suppressMessages(AnnotationDbi::mapIds(org, ashm_symbols[todo], "ENTREZID",
                                                             "ALIAS", multiVals = "first"))
  matched = ifelse(ashm_symbols %in% spans$symbol, ashm_symbols,
                   spans$symbol[match(ashm_entrez, spans$entrez)])
  if(any(is.na(matched))){
    message(build, ": aSHM genes without coordinates: ", paste(ashm_symbols[is.na(matched)], collapse = ", "))
  }
  ashm_genes = unique(stats::na.omit(matched))
  ashm_tss = tx %>%
    filter(gene %in% ashm_genes) %>%
    distinct(gene, entrez, chrom, strand, tss)
  list(build = build, spans = spans, longest = longest, primary = primary, catalogue = catalogue,
       ashm_symbols = ashm_symbols, matched = matched, ashm_genes = ashm_genes, ashm_tss = ashm_tss)
}

# aSHM gene TSSs lifted from one build to the other, kept if they land on the
# same gene (by Entrez ID) and strand, within 100 kb of its span
lift_tss = function(from, to){
  chain = rtracklayer::import.chain(system.file("extdata",
    if(from$build == "grch37") "hg19ToHg38.over.chain" else "hg38ToHg19.over.chain",
    package = "GAMBLR.data"))
  x = from$ashm_tss
  gr = GenomicRanges::GRanges(paste0("chr", x$chrom), IRanges::IRanges(x$tss, width = 1), strand = x$strand)
  lifted = rtracklayer::liftOver(gr, chain)
  ok = lengths(lifted) == 1
  lifted = unlist(lifted[ok])
  data.frame(entrez = x$entrez[ok], chrom = sub("^chr", "", as.character(GenomicRanges::seqnames(lifted))),
             strand = as.character(GenomicRanges::strand(lifted)), tss = GenomicRanges::start(lifted)) %>%
    inner_join(to$spans %>% select(gene = symbol, entrez, gene_chrom = chrom, start, end), by = "entrez") %>%
    inner_join(to$ashm_tss %>% distinct(gene, gene_strand = strand), by = "gene") %>%
    filter(chrom == gene_chrom, strand == gene_strand, tss >= start - 1e5, tss <= end + 1e5) %>%
    distinct(gene, entrez, chrom, strand, tss)
}

# one build's table: the primary TSS of every gene, and for aSHM genes every
# TSS (its own transcripts' and those lifted from the other build) with the
# genomes mutated in its window
build_tss = function(prep, lifted){
  build = prep$build
  spans = prep$spans
  longest = prep$longest
  primary = prep$primary
  catalogue = prep$catalogue
  ashm_symbols = prep$ashm_symbols
  matched = prep$matched
  ashm_genes = prep$ashm_genes
  candidates = bind_rows(prep$ashm_tss, lifted) %>%
    distinct(gene, chrom, strand, tss) %>%
    left_join(longest %>% select(gene, longest_tss = tss), by = "gene") %>%
    mutate(name = paste(gene, tss, sep = "_"),
           w_start = ifelse(strand == "+", tss + count_window[1], tss - count_window[2]),
           w_end = ifelse(strand == "+", tss + count_window[2], tss - count_window[1]))
  # mutations in each gene's curated aSHM regions, then counted in each of
  # its TSS windows: only the parts of a window inside the gene's regions
  # count, so coding and other mutations bundled for the gene (e.g. the
  # NFKBIZ 3' UTR hotspot) don't
  regions = catalogue %>%
    mutate(gene = matched[match(gene, ashm_symbols)]) %>%
    filter(!is.na(gene))
  bed = data.frame(chrom = if(build == "hg38") regions$chr_name else sub("^chr", "", regions$chr_name),
                   start = regions[[2]], end = regions[[3]], name = regions$gene)
  muts = as.data.frame(GAMBLR.open::get_ssm_by_regions(
    these_samples_metadata = genomes,
    regions_bed = GAMBLR.utils::create_bed_data(bed, genome_build = build),
    projection = build, streamlined = FALSE, use_name_column = TRUE)) %>%
    distinct(gene = region_name, sample = Tumor_Sample_Barcode, pos = Start_Position)
  counts = candidates %>%
    select(name, gene, w_start, w_end) %>%
    inner_join(muts, by = "gene", relationship = "many-to-many") %>%
    filter(pos >= w_start, pos <= w_end) %>%
    group_by(name) %>%
    summarise(n_mut = n(), n_samples = n_distinct(sample), .groups = "drop")
  # the most mutated first (ties: more mutations, then the longest isoform's
  # TSS, then the most upstream)
  by_mutations = function(x){
    arrange(x, desc(n_samples), desc(n_mut), tss != longest_tss, ifelse(strand == "+", tss, -tss),
            .by_group = TRUE)
  }
  ashm = candidates %>%
    left_join(counts, by = "name") %>%
    mutate(n_mut = coalesce(n_mut, 0L), n_samples = coalesce(n_samples, 0L)) %>%
    # TSSs within merge_gap of each other share most of their window: keep
    # the most mutated of each such group
    group_by(gene) %>%
    arrange(tss, .by_group = TRUE) %>%
    mutate(group = cumsum(c(TRUE, diff(tss) > merge_gap))) %>%
    group_by(gene, group) %>%
    by_mutations() %>%
    slice(1) %>%
    # primary: the most mutated
    group_by(gene) %>%
    by_mutations() %>%
    mutate(primary = row_number() == 1,
           active = n_samples >= min_samples & n_samples >= min_fraction * max(n_samples),
           tss_source = ifelse(primary & n_samples == 0, "longest_isoform", "ashm_mutations")) %>%
    ungroup() %>%
    select(gene, strand, tss, tss_source, primary, active, n_samples, n_mut)

  bind_rows(primary %>% filter(!gene %in% ashm_genes), ashm) %>%
    left_join(spans %>% select(gene = symbol, entrez, chrom, start, end), by = "gene") %>%
    mutate(genome_build = build, ashm = gene %in% ashm_genes) %>%
    select(genome_build, gene, entrez, chrom, start, end, strand, tss, tss_source, primary, ashm,
           active, n_samples, n_mut) %>%
    arrange(gene, desc(primary), tss)
}

prep = list(grch37 = prepare_build("grch37"), hg38 = prepare_build("hg38"))
gene_tss = bind_rows(build_tss(prep$grch37, lift_tss(prep$hg38, prep$grch37)),
                     build_tss(prep$hg38, lift_tss(prep$grch37, prep$hg38)))
gene_tss = as.data.frame(gene_tss)
attr(gene_tss, "counted_genomes") = nrow(genomes)
attr(gene_tss, "built") = as.character(Sys.Date())

usethis::use_data(gene_tss, overwrite = TRUE, compress = "xz")
