#' Gene transcription start sites
#'
#' The transcription start sites (TSS) used by [get_gene_tss] and
#' [subset_ssm_by_tss], for grch37 and hg38, built by `data-raw/gene_tss.R`
#' from UCSC knownGene (UCSC Genes for grch37, GENCODE V44 for hg38).
#'
#' Every gene in GAMBLR.data's gene coordinates (with an Entrez ID and a
#' strand) has one primary TSS: that of its longest isoform (longest spliced
#' transcript), or the start of its gene span if knownGene has no transcript
#' for it. Genes in the aSHM catalogue
#' (`GAMBLR.data::somatic_hypermutation_locations_*_v_latest`) instead have
#' one row per distinct promoter (transcript TSSs within 2 kb of each other
#' merged), each with the number of GAMBLR.open genomes mutated in its first
#' 4 kb, counting only mutations in the gene's curated aSHM regions (so a
#' window reaching e.g. a 3' UTR hotspot doesn't count it). Their primary TSS
#' is the most mutated, and a TSS is active when at
#' least 10 genomes, and at least 10% as many as the primary TSS, are mutated
#' there. Only aSHM genes are counted because GAMBLR.open bundles non-coding
#' mutations only for them.
#'
#' @format A data frame with columns
#' \describe{
#'   \item{genome_build}{"grch37" or "hg38".}
#'   \item{gene, entrez}{Gene symbol (as in that build's gene coordinates)
#'     and Entrez ID.}
#'   \item{chrom, start, end}{The gene span (chromosome without "chr").}
#'   \item{strand, tss}{Strand ("+" or "-") and TSS position.}
#'   \item{tss_source}{"longest_isoform", "gene_span" or "ashm_mutations"
#'     (chosen by mutation counts).}
#'   \item{primary}{TRUE for the gene's primary TSS.}
#'   \item{ashm}{TRUE for genes in the aSHM catalogue.}
#'   \item{active, n_samples, n_mut}{For aSHM genes: whether the TSS is
#'     active, and the genomes and mutations in its first 4 kb; NA
#'     otherwise.}
#' }
"gene_tss"

#' @title Get gene transcription start sites
#'
#' @description Return the transcription start site(s) (TSS) and strand of
#' genes, from the bundled [gene_tss] table.
#'
#' @details Most genes have one TSS: that of their longest isoform. Genes in
#' the aSHM catalogue (`GAMBLR.data::somatic_hypermutation_locations_*`)
#' have one per distinct promoter instead, ranked by how many GAMBLR.open
#' genomes are mutated in their first 4 kb: the most mutated is the gene's
#' primary TSS, and the well-mutated ones are active (see [gene_tss]). This
#' picks the hypermutated promoter where it isn't the longest isoform's (e.g.
#' ETS1, SGK1, BLK), and finds genes with more than one (e.g. PAX5, IRF8).
#' `which` only applies to aSHM genes; other genes always give their one TSS.
#' Only genes on chromosomes 1-22, X and Y are included.
#'
#' @param genes Optional gene symbols. Default NULL returns every gene.
#' Symbols are matched to the table's gene names, then by Entrez ID (so
#' current and former names both work, e.g. "TSPOAP1" finds BZRAP1 in
#' grch37) and through [expand_gene_aliases] (e.g. "H1-3" finds HIST1H1D).
#' @param genome_build "grch37" (default) or "hg38".
#' @param which For aSHM genes: "primary" (default; one TSS per gene, the
#' most mutated), "active" (every well-mutated TSS, always including the
#' primary one) or "all" (every distinct TSS).
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
#' @return With "df", a data frame with one row per TSS: gene (as
#' requested), chrom (without a "chr" prefix), start and end (the gene),
#' strand ("+" or "-"), tss, tss_source ("longest_isoform", "gene_span" or
#' "ashm_mutations"), primary, active and n_samples (the last two NA for
#' genes outside the aSHM catalogue), plus window_start and window_end
#' (genome coordinates) when `window` is given. With "region" or "bed", the
#' window (or TSS) of each row as described for `return_as`.
#'
#' @export
#'
#' @examples
#' get_gene_tss(c("BCL2", "MYC", "PIM1"))
#'
#' # ETS1's hypermutated promoter is not its longest isoform's
#' get_gene_tss("ETS1", which = "all")
#'
#' # PAX5 has two hypermutated promoters
#' get_gene_tss("PAX5", which = "active")
#'
#' # regions covering the first 4 kb of each gene, to fetch the mutations
#' # there (e.g. with GAMBLR.open::get_ssm_by_regions)
#' get_gene_tss(c("BCL2", "MYC", "PIM1"), window = c(0, 4000), return_as = "region")
get_gene_tss = function(genes = NULL, genome_build = "grch37", window = NULL,
                        return_as = "df", which = "primary"){
  return_as = match.arg(return_as, c("df", "region", "bed"))
  which = match.arg(which, c("primary", "active", "all"))
  if(!is.null(window) && (!is.numeric(window) || length(window) != 2 || window[1] > window[2])){
    stop("window must be c(start, end) with start <= end, in bp relative to the TSS")
  }
  build = .maf_data_build(genome_build)
  if(!build %in% c("grch37", "hg38")){
    stop("genome_build must be grch37 or hg38")
  }
  tab = GAMBLR.utils::gene_tss
  tab = tab[tab$genome_build == build, , drop = FALSE]
  keep = switch(which,
                primary = tab$primary,
                active = tab$primary | tab$active %in% TRUE,
                all = rep(TRUE, nrow(tab)))
  tab = tab[keep, , drop = FALSE]

  if(!is.null(genes)){
    genes = unique(genes)
    # by name, then Entrez ID, then alias
    found = genes[genes %in% tab$gene]
    names(found) = found
    todo = setdiff(genes, found)
    if(length(todo) > 0){
      org = org.Hs.eg.db::org.Hs.eg.db
      eg = suppressMessages(AnnotationDbi::mapIds(org, todo, "ENTREZID", "SYMBOL", multiVals = "first"))
      if(any(is.na(eg))){
        eg[is.na(eg)] = suppressMessages(AnnotationDbi::mapIds(org, todo[is.na(eg)], "ENTREZID",
                                                               "ALIAS", multiVals = "first"))
      }
      by_entrez = tab$gene[match(eg, tab$entrez)]
      for(i in seq_along(todo)){
        hit = by_entrez[i]
        if(is.na(hit)) hit = intersect(expand_gene_aliases(todo[i]), tab$gene)[1]
        if(!is.na(hit)) found[todo[i]] = hit
      }
    }
    missing = setdiff(genes, names(found))
    if(length(missing) > 0){
      message(length(missing), " gene(s) not in the ", build, " TSS table: ",
              paste(utils::head(missing, 10), collapse = ", "),
              if(length(missing) > 10) ", ..." else "")
    }
    found = found[intersect(genes, names(found))]
    rows = lapply(found, function(g) which(tab$gene == g))
    tab = tab[unlist(rows), , drop = FALSE]
    tab$gene = rep(names(found), lengths(rows))
  }
  out = data.frame(gene = tab$gene, chrom = tab$chrom, start = tab$start, end = tab$end,
                   strand = tab$strand, tss = tab$tss, tss_source = tab$tss_source,
                   primary = tab$primary, active = tab$active, n_samples = tab$n_samples,
                   stringsAsFactors = FALSE)
  if(nrow(out) == 0 && return_as == "region") return(character(0))
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
#' Genes in the aSHM catalogue can have more than one hypermutated promoter
#' (see [get_gene_tss]); by default (`which = "active"`) a mutation counts if
#' it is in the window of any of them, and is assigned to the nearest.
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
#' @param which Which TSSs of aSHM genes to use, as for [get_gene_tss]:
#' "active" (default; every well-mutated promoter), "primary" or "all".
#' Other genes have one TSS whatever this is.
#' @param tss Optional TSS table to use instead of [get_gene_tss], with
#' columns gene, chrom, strand ("+" or "-") and tss (several rows per gene
#' are allowed).
#'
#' @return The rows of `maf_data` in a window, with three added columns:
#' `tss_gene` (the gene whose window it's in), `tss_position` (the TSS of
#' that window) and `tss_distance` (its position relative to that TSS, as
#' for `start` and `end`). A `maf_data` input gives a `maf_data` result (see
#' [create_maf_data]).
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
                             which = "active",
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
    tss = get_gene_tss(genes, genome_build, which = which)
  }else{
    tss = as.data.frame(tss)
    if(!all(c("gene", "chrom", "strand", "tss") %in% colnames(tss))){
      stop("tss needs columns gene, chrom, strand and tss")
    }
    tss$chrom = sub("^chr", "", as.character(tss$chrom))
    if(!is.null(genes)) tss = tss[tss$gene %in% genes, , drop = FALSE]
  }
  out_cols = c(colnames(maf_data), "tss_gene", "tss_position", "tss_distance")
  empty = cbind(maf_data[0, , drop = FALSE], tss_gene = character(), tss_position = numeric(),
                tss_distance = numeric())
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
  hits$distance = ifelse(hits$strand == "+", hits$i.start - hits$gene_tss,
                         hits$gene_tss - hits$i.start)
  # once per gene: a mutation in the windows of several of a gene's TSSs
  # goes to the nearest
  hits = hits[order(hits$row, hits$gene, abs(hits$distance)), ]
  hits = hits[!duplicated(hits[, c("row", "gene")]), ]
  out = maf_data[hits$row, , drop = FALSE]
  out$tss_gene = hits$gene
  out$tss_position = hits$gene_tss
  out$tss_distance = hits$distance
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
