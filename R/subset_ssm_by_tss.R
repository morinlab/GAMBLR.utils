#' @title Get gene transcription start sites
#'
#' @description Return the transcription start site (TSS) and strand of
#' genes, from the gene coordinates in GAMBLR.data
#' (`grch37_gene_coordinates` or `hg38_gene_coordinates`) and the strand in
#' \pkg{org.Hs.eg.db}.
#'
#' @details The gene coordinates are Ensembl gene spans (the union of the
#' gene's transcripts), so the TSS is the most upstream transcription start:
#' the gene start on the + strand and the gene end on the - strand. Only
#' genes on chromosomes 1-22, X and Y are returned. The strand is looked up
#' by Ensembl gene ID, then by symbol or alias; genes without one are left
#' out with a message.
#'
#' @param genes Optional gene symbols. Default NULL returns every gene with
#' coordinates and a strand. Symbols are matched to the coordinate table's
#' `hugo_symbol` or `gene_name`, then through [expand_gene_aliases] (so e.g.
#' "H1-3" finds HIST1H1D in grch37).
#' @param genome_build "grch37" (default) or "hg38".
#' @param window Optional window around each TSS, \code{c(start, end)} in bp
#' as for [subset_ssm_by_tss] (0 is the TSS, positive values are into the
#' gene, negative values upstream), e.g. \code{c(0, 4000)} for the first 4 kb
#' of the gene. Default NULL: the TSS itself.
#' @param return_as "df" (default), "region" or "bed", as for
#' [gene_to_region]: "region" gives a named vector of "chrom:start-end"
#' strings (with a "chr" prefix for hg38), e.g. for
#' `GAMBLR.open::get_ssm_by_regions(regions_list = ...)`, and "bed" a data
#' frame with chromosome, start, end and hugo_symbol. Both cover `window`
#' (or the TSS alone).
#'
#' @return With "df", a data frame with one row per gene: gene (as
#' requested, or the coordinate table's symbol), chrom (without a "chr"
#' prefix), start and end (the gene), strand ("+" or "-") and tss, plus
#' window_start and window_end (genome coordinates) when `window` is given.
#' With "region" or "bed", the window (or TSS) of each gene as described
#' for `return_as`.
#'
#' @export
#'
#' @examples
#' get_gene_tss(c("BCL2", "MYC", "PIM1"))
#'
#' # regions covering the first 4 kb of each gene, to fetch the mutations
#' # there (e.g. with GAMBLR.open::get_ssm_by_regions)
#' get_gene_tss(c("BCL2", "MYC", "PIM1"), window = c(0, 4000), return_as = "region")
get_gene_tss = function(genes = NULL, genome_build = "grch37", window = NULL,
                        return_as = "df"){
  return_as = match.arg(return_as, c("df", "region", "bed"))
  if(!is.null(window) && (!is.numeric(window) || length(window) != 2 || window[1] > window[2])){
    stop("window must be c(start, end) with start <= end, in bp relative to the TSS")
  }
  build = tolower(genome_build)
  if(build %in% c("grch37", "hg19")){
    coords = GAMBLR.data::grch37_gene_coordinates
  }else if(build %in% c("hg38", "grch38")){
    coords = GAMBLR.data::hg38_gene_coordinates
  }else{
    stop("genome_build must be grch37 or hg38")
  }
  coords = as.data.frame(coords)
  coords$chromosome = sub("^chr", "", coords$chromosome)
  coords = coords[coords$chromosome %in% c(1:22, "X", "Y"), , drop = FALSE]
  coords$symbol = ifelse(is.na(coords$hugo_symbol), coords$gene_name, coords$hugo_symbol)

  if(!is.null(genes)){
    genes = unique(genes)
    rows = lapply(genes, function(g){
      hit = which(coords$hugo_symbol %in% g | coords$gene_name %in% g)
      if(length(hit) == 0){
        aliases = setdiff(expand_gene_aliases(g), g)
        hit = which(coords$hugo_symbol %in% aliases | coords$gene_name %in% aliases)
      }
      hit
    })
    missing = genes[lengths(rows) == 0]
    if(length(missing) > 0){
      message(length(missing), " gene(s) not in the ", build, " gene coordinates: ",
              paste(utils::head(missing, 10), collapse = ", "),
              if(length(missing) > 10) ", ..." else "")
    }
    coords = coords[unlist(rows), , drop = FALSE]
    coords$symbol = rep(genes, lengths(rows))
  }

  # strand from org.Hs.eg.db: Ensembl ID, then symbol, then alias
  entrez = suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = coords$ensembl_gene_id, column = "ENTREZID",
    keytype = "ENSEMBL", multiVals = "first"))
  for(keytype in c("SYMBOL", "ALIAS")){
    todo = is.na(entrez)
    if(!any(todo)) break
    entrez[todo] = suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = coords$symbol[todo], column = "ENTREZID",
      keytype = keytype, multiVals = "first"))
  }
  loc = AnnotationDbi::mget(unique(stats::na.omit(entrez)), org.Hs.eg.db::org.Hs.egCHRLOC,
                            ifnotfound = NA)
  sign = vapply(loc, function(x) if(all(is.na(x))) NA_real_ else sign(x[!is.na(x)][1]), numeric(1))
  coords$strand = ifelse(sign[as.character(entrez)] > 0, "+",
                         ifelse(sign[as.character(entrez)] < 0, "-", NA))
  no_strand = unique(coords$symbol[is.na(coords$strand)])
  if(length(no_strand) > 0 && !is.null(genes)){
    message(length(no_strand), " gene(s) without a strand in org.Hs.eg.db, left out: ",
            paste(utils::head(no_strand, 10), collapse = ", "),
            if(length(no_strand) > 10) ", ..." else "")
  }
  coords = coords[!is.na(coords$strand), , drop = FALSE]

  # one row per gene: its span on its (first) chromosome
  out = do.call(rbind, lapply(split(coords, coords$symbol), function(x){
    x = x[x$chromosome == x$chromosome[1], , drop = FALSE]
    data.frame(gene = x$symbol[1], chrom = x$chromosome[1], start = min(x$start),
               end = max(x$end), strand = x$strand[1], stringsAsFactors = FALSE)
  }))
  if(is.null(out)){
    if(return_as == "region") return(character(0))
    return(data.frame(gene = character(), chrom = character(), start = numeric(),
                      end = numeric(), strand = character(), tss = numeric()))
  }
  out$tss = ifelse(out$strand == "+", out$start, out$end)
  if(!is.null(genes)) out = out[order(match(out$gene, genes)), , drop = FALSE]
  rownames(out) = NULL
  if(return_as == "df" && is.null(window)) return(out)
  w = .tss_window(out, if(is.null(window)) c(0, 0) else window)
  if(return_as == "df"){
    out$window_start = w$start
    out$window_end = w$end
    return(out)
  }
  chrom = if(build %in% c("hg38", "grch38")) paste0("chr", out$chrom) else out$chrom
  if(return_as == "bed"){
    return(data.frame(chromosome = chrom, start = w$start, end = w$end,
                      hugo_symbol = out$gene, stringsAsFactors = FALSE))
  }
  stats::setNames(paste0(chrom, ":", w$start, "-", w$end), out$gene)
}

# genome coordinates of a window c(start, end) relative to each TSS (0 is
# the TSS, positive values downstream in the direction of transcription)
.tss_window = function(tss, window){
  plus = tss$strand == "+"
  list(start = ifelse(plus, tss$tss + window[1], tss$tss - window[2]),
       end = ifelse(plus, tss$tss + window[2], tss$tss - window[1]))
}

#' @title Subset mutations to a window around gene transcription start sites
#'
#' @description Keep the mutations (by default the non-coding ones) that lie
#' within a window around the transcription start site (TSS) of each gene,
#' e.g. the first 4 kb of the gene (the LymphGen "Synon" definition, Wright
#' et al. 2020), where aberrant somatic hypermutation concentrates.
#'
#' @details The window is set relative to the TSS in the direction of
#' transcription: 0 is the TSS, positive values are downstream (into the
#' gene) and negative values upstream. `start = 0, end = 4000` keeps the
#' first 4 kb of the gene; `start = -2000, end = 4000` also keeps 2 kb of
#' the promoter. Both ends are inclusive.
#'
#' Mutations are assigned to genes by position, not by their
#' `Hugo_Symbol`, so an upstream mutation annotated to a neighbouring gene
#' (e.g. as 5' flank) still counts for the gene whose window it falls in.
#' A mutation in the windows of two genes (e.g. a bidirectional promoter)
#' gives one row for each. Use a source that includes non-coding mutations,
#' e.g. `get_ssm_by_genes` (GAMBLR.open or GAMBLR.results), or
#' `get_ssm_by_regions` with the windows from [get_gene_tss] as regions; a
#' coding-only MAF only has silent and UTR mutations to offer. Capture data covers little beyond exons and UTRs, so windows are
#' only partly assessed there.
#'
#' @param maf_data MAF data frame (Chromosome, Start_Position,
#' Variant_Classification).
#' @param genes Optional gene symbols whose windows to use. Default NULL
#' uses every gene from [get_gene_tss] (or in `tss`).
#' @param start,end The window, in bp relative to the TSS (see Details).
#' Default 0 and 4000.
#' @param genome_build Genome build of `maf_data`, "grch37" or "hg38".
#' Default: `maf_data`'s `genome_build` attribute if it has one (as GAMBLR's
#' MAF getters set), otherwise "grch37".
#' @param noncoding_only Keep only mutations that don't change the protein:
#' `Variant_Classification` not in `GAMBLR.helpers::vc_nonSynonymous` (i.e.
#' silent, UTR, flank, intron and similar). Default TRUE; FALSE keeps every
#' mutation in the window.
#' @param tss Optional TSS table to use instead of [get_gene_tss], with
#' columns gene, chrom, strand ("+" or "-") and tss.
#'
#' @return The rows of `maf_data` in a window, with two added columns:
#' `tss_gene` (the gene whose window it's in) and `tss_distance` (its
#' position relative to that gene's TSS, as for `start` and `end`). A
#' `maf_data` input gives a `maf_data` result (see [create_maf_data]).
#'
#' @export
#'
#' @examples
#' meta = GAMBLR.open::get_gambl_metadata()
#' meta = head(meta[meta$pathology == "DLBCL" & meta$seq_type == "genome", ], 40)
#' genes = c("BCL2", "PIM1", "MYC")
#'
#' # every mutation in the genes, coding and non-coding
#' maf = GAMBLR.open::get_ssm_by_genes(genes, these_samples_metadata = meta)
#'
#' # the non-coding ones in the first 4 kb of each gene (LymphGen's Synon
#' # window), assigned to genes by position
#' synon = subset_ssm_by_tss(maf, genes = genes)
#' table(synon$tss_gene, synon$Variant_Classification)
#'
#' # 2 kb upstream to 6 kb downstream instead
#' wider = subset_ssm_by_tss(maf, genes = genes, start = -2000, end = 6000)
#' range(wider$tss_distance)
subset_ssm_by_tss = function(maf_data,
                             genes = NULL,
                             start = 0,
                             end = 4000,
                             genome_build = "grch37",
                             noncoding_only = TRUE,
                             tss = NULL){
  if(!is.numeric(start) || !is.numeric(end) || length(start) != 1 || length(end) != 1 || start > end){
    stop("start and end must be single numbers with start <= end")
  }
  if(missing(genome_build) && !is.null(get_genome_build(maf_data))){
    genome_build = get_genome_build(maf_data)
  }
  # a maf_data input is returned as maf_data (as.data.frame and subsetting
  # drop its class and genome_build)
  as_maf_data = inherits(maf_data, "maf_data")
  restore = function(x){
    if(!as_maf_data) return(x)
    create_maf_data(x, .maf_data_build(genome_build))
  }
  maf_data = as.data.frame(maf_data)
  needed = c("Chromosome", "Start_Position", "Variant_Classification")
  if(!all(needed %in% colnames(maf_data))){
    stop("maf_data needs columns ", paste(needed, collapse = ", "))
  }
  if(is.null(tss)){
    tss = get_gene_tss(genes, genome_build)
  }else{
    tss = as.data.frame(tss)
    if(!all(c("gene", "chrom", "strand", "tss") %in% colnames(tss))){
      stop("tss needs columns gene, chrom, strand and tss")
    }
    tss$chrom = sub("^chr", "", as.character(tss$chrom))
    if(!is.null(genes)) tss = tss[tss$gene %in% genes, , drop = FALSE]
  }
  out_cols = c(colnames(maf_data), "tss_gene", "tss_distance")
  empty = cbind(maf_data[0, , drop = FALSE], tss_gene = character(), tss_distance = numeric())
  if(noncoding_only){
    maf_data = maf_data[!maf_data$Variant_Classification %in% GAMBLR.helpers::vc_nonSynonymous, ,
                        drop = FALSE]
  }
  if(nrow(maf_data) == 0 || nrow(tss) == 0) return(restore(empty))

  # windows in genome coordinates
  w = .tss_window(tss, c(start, end))
  windows = data.table::data.table(
    chrom = tss$chrom, gene = tss$gene, gene_tss = tss$tss, strand = tss$strand,
    start = w$start, end = w$end)
  muts = data.table::data.table(
    chrom = sub("^chr", "", as.character(maf_data$Chromosome)),
    start = as.numeric(maf_data$Start_Position), end = as.numeric(maf_data$Start_Position),
    row = seq_len(nrow(maf_data)))
  data.table::setkey(windows, chrom, start, end)
  hits = data.table::foverlaps(muts, windows, type = "within", nomatch = NULL)
  if(nrow(hits) == 0) return(restore(empty))
  hits = hits[order(hits$row), ]
  out = maf_data[hits$row, , drop = FALSE]
  out$tss_gene = hits$gene
  out$tss_distance = ifelse(hits$strand == "+", hits$i.start - hits$gene_tss,
                            hits$gene_tss - hits$i.start)
  rownames(out) = NULL
  restore(out[, out_cols])
}

# genome build names as create_maf_data() accepts them ("grch37" or "hg38")
.maf_data_build = function(genome_build){
  build = tolower(genome_build)
  if(build %in% c("hg19", "grch37", "hs37d5")) return("grch37")
  if(build %in% c("hg38", "grch38")) return("hg38")
  build
}
