#' @title Get gene models (strand, exons) from mygene.info
#'
#' @description Retrieve gene coordinates, strand and exon structure for a
#' set of genes, or for every gene in a region, from the mygene.info web
#' service (RefSeq transcripts). The result can be passed to the
#' `gene_models` argument of GAMBLR.viz region plots (e.g.
#' `GAMBLR.viz::plot_ssm_density`) to draw exons, introns and strand.
#'
#' @details Queries `https://mygene.info/v3/query` (no key needed), so an
#' internet connection is required; results are cached for the R session
#' (and on disk when `cache_dir` is given), so repeated plots don't query
#' again. For grch37 the hg19 coordinates are used, for hg38 the current
#' ones. Genes on alternative contigs are left out, and hits without a gene
#' symbol are skipped.
#'
#' `transcripts` chooses which exons are returned:
#' \itemize{
#'   \item "union" (default): the exons of all of a gene's transcripts
#'     merged, i.e. all exonic sequence of the gene (transcript = "union").
#'   \item "longest": the transcript spanning the most sequence.
#'   \item "all": every transcript.
#' }
#' Exons are split into their coding and untranslated parts (`feature`:
#' "coding", "UTR", or "noncoding" for genes without a coding transcript),
#' from each transcript's coding start and end; with "union", a stretch is
#' coding if it is coding in any transcript. A gene's start and end are those
#' of its transcripts' exons (RefSeq), which can differ from the gene record
#' (e.g. Ensembl's, which may include other transcripts). Coordinates are as
#' given by mygene.info (UCSC-style, 0-based starts), so they can differ by a
#' base from gene coordinates elsewhere.
#'
#' @param genes Gene symbols. Give this or `region`.
#' @param region A region as "chrom:start-end" (e.g. "3:187000000-189000000");
#' all genes overlapping it are returned.
#' @param genome_build "grch37" (default; hg19 coordinates) or "hg38".
#' @param transcripts "union" (default), "longest" or "all". See Details.
#' @param cache_dir Optional directory for an on-disk cache of responses
#' (e.g. `tools::R_user_dir("GAMBLR", "cache")`), for reuse across sessions
#' and offline.
#' @param timeout Seconds to wait for the service. Default 30.
#'
#' @return A data frame with one row per exon part (genes without exon data
#' get one row with NA exon coordinates): gene, chrom, start, end (the gene),
#' strand ("+"/"-"), transcript, exon_start, exon_end and feature ("coding",
#' "UTR" or "noncoding").
#'
#' @export
#'
#' @examples
#' \dontrun{
#' models = get_gene_models(c("BCL6", "LPP"), genome_build = "grch37")
#' models = get_gene_models(region = "3:187000000-189000000")
#' GAMBLR.viz::plot_ssm_density(maf, genes = "BCL6", flank = 3e4,
#'                              gene_models = models)
#' }
get_gene_models = function(genes = NULL, region = NULL, genome_build = "grch37",
                           transcripts = "union", cache_dir = NULL, timeout = 30){
  if(!requireNamespace("jsonlite", quietly = TRUE)){
    stop("get_gene_models() needs the jsonlite package")
  }
  transcripts = match.arg(transcripts, c("union", "longest", "all"))
  build = tolower(genome_build)
  if(build %in% c("grch37", "hg19")){
    pos_field = "genomic_pos_hg19"; exon_field = "exons_hg19"; prefix = "hg19."
  }else if(build %in% c("hg38", "grch38")){
    pos_field = "genomic_pos"; exon_field = "exons"; prefix = ""
  }else{
    stop("genome_build must be grch37 or hg38")
  }
  if(is.null(genes) == is.null(region)){
    stop("give genes or region")
  }
  region_chrom = NULL
  if(!is.null(region)){
    clean = gsub(",", "", region)
    parts = regmatches(clean, regexec("^(chr)?([0-9XY]+):([0-9.e+]+)-([0-9.e+]+)$", clean))[[1]]
    if(length(parts) == 0) stop("region must look like \"3:187000000-189000000\"")
    region_chrom = parts[3]
    queries = paste0(prefix, "chr", parts[3], ":", format(as.numeric(parts[4]), scientific = FALSE),
                     "-", format(as.numeric(parts[5]), scientific = FALSE))
  }else{
    genes = unique(as.character(genes))
    chunks = split(genes, ceiling(seq_along(genes) / 50))
    queries = vapply(chunks, function(g) paste0("symbol:", g, collapse = " OR "), character(1))
  }
  fields = paste("symbol", pos_field, exon_field, sep = ",")
  hits = unlist(lapply(queries, function(q){
    url = paste0("https://mygene.info/v3/query?q=", utils::URLencode(q, reserved = TRUE),
                 "&species=human&size=1000&fields=", fields)
    .mygene_get(url, cache_dir, timeout)$hits
  }), recursive = FALSE)

  main = c(as.character(1:22), "X", "Y")
  pick = function(x){
    # a single record or a list of them (e.g. alternative contigs)
    if(is.null(x)) return(list())
    if(!is.null(x$chr)) list(x) else x
  }
  rows = lapply(hits, function(h){
    symbol = h$symbol
    if(is.null(symbol)) return(NULL)
    if(!is.null(genes) && !toupper(symbol) %in% toupper(genes)) return(NULL)
    pos = Filter(function(p) as.character(p$chr) %in% main &&
                   (is.null(region_chrom) || as.character(p$chr) == region_chrom),
                 pick(h[[pos_field]]))
    if(length(pos) == 0) return(NULL)
    pos = pos[[1]]
    strand = if(is.null(pos$strand)) NA_character_ else ifelse(pos$strand > 0, "+", "-")
    gene = data.frame(gene = symbol, chrom = as.character(pos$chr), start = pos$start,
                      end = pos$end, strand = strand, stringsAsFactors = FALSE)
    tx = Filter(function(t) as.character(t$chr) == gene$chrom, pick(h[[exon_field]]))
    if(length(tx) == 0){
      return(cbind(gene, transcript = NA_character_, exon_start = NA_real_, exon_end = NA_real_,
                   feature = NA_character_))
    }
    # each transcript's exons, and its coding part (cdsstart..cdsend; none
    # for non-coding transcripts)
    txs = lapply(tx, function(t){
      p = do.call(rbind, lapply(t$position, function(e) as.numeric(unlist(e))))
      cds = c(as.numeric(t$cdsstart), as.numeric(t$cdsend))
      list(id = as.character(t$transcript), exons = p[, 1:2, drop = FALSE],
           cds = if(length(cds) == 2 && !anyNA(cds) && cds[2] > cds[1]) cds else NULL,
           span = as.numeric(t$txend) - as.numeric(t$txstart))
    })
    if(transcripts == "longest"){
      # longest span; ties go to the transcript with more exons
      span = vapply(txs, `[[`, numeric(1), "span")
      n_exons = vapply(txs, function(t) nrow(t$exons), numeric(1))
      txs = txs[order(-span, -n_exons)[1]]
    }
    coding_gene = any(!vapply(txs, function(t) is.null(t$cds), logical(1)))
    pieces = function(id, exon, coding){
      # exon intervals split into coding and UTR (or non-coding) parts
      cds = .interval_intersect(exon, coding)
      utr = .interval_subtract(exon, cds)
      rbind(if(nrow(cds) > 0) data.frame(transcript = id, exon_start = cds[, 1],
                                         exon_end = cds[, 2], feature = "coding"),
            if(nrow(utr) > 0) data.frame(transcript = id, exon_start = utr[, 1],
                                         exon_end = utr[, 2],
                                         feature = if(coding_gene) "UTR" else "noncoding"))
    }
    coding_of = function(t) if(is.null(t$cds)) matrix(numeric(0), 0, 2) else
      .interval_intersect(t$exons, matrix(t$cds, 1))
    exons = if(transcripts == "union"){
      pieces("union", .interval_union(do.call(rbind, lapply(txs, `[[`, "exons"))),
             .interval_union(do.call(rbind, lapply(txs, coding_of))))
    }else{
      do.call(rbind, lapply(txs, function(t) pieces(t$id, t$exons, coding_of(t))))
    }
    exons = unique(exons[order(exons$transcript, exons$exon_start), , drop = FALSE])
    # the gene's extent from its transcripts, so the drawn model and the
    # coordinates agree (the gene record can span other annotations)
    gene$start = min(exons$exon_start); gene$end = max(exons$exon_end)
    cbind(gene[rep(1, nrow(exons)), , drop = FALSE], exons)
  })
  out = do.call(rbind, rows)
  if(is.null(out)){
    message("no genes found")
    return(data.frame(gene = character(0), chrom = character(0), start = numeric(0),
                      end = numeric(0), strand = character(0), transcript = character(0),
                      exon_start = numeric(0), exon_end = numeric(0), feature = character(0)))
  }
  if(!is.null(genes)){
    missing = setdiff(toupper(genes), toupper(out$gene))
    if(length(missing) > 0){
      message("not found on a main chromosome: ", paste(utils::head(missing, 10), collapse = ", "))
    }
  }
  out = out[order(out$chrom, out$start, out$gene, out$exon_start), , drop = FALSE]
  rownames(out) = NULL
  out
}

# Cached GET of a mygene.info query (session cache, plus optional disk cache)
.gene_model_cache = new.env(parent = emptyenv())
.mygene_get = function(url, cache_dir = NULL, timeout = 30){
  key = url
  if(exists(key, envir = .gene_model_cache, inherits = FALSE)){
    return(get(key, envir = .gene_model_cache))
  }
  file = NULL
  if(!is.null(cache_dir)){
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    file = file.path(cache_dir, paste0("mygene_", .url_hash(url), ".json"))
    if(file.exists(file)){
      value = jsonlite::fromJSON(file, simplifyVector = FALSE)
      assign(key, value, envir = .gene_model_cache)
      return(value)
    }
  }
  old = options(timeout = timeout)
  on.exit(options(old))
  text = tryCatch(paste(readLines(url, warn = FALSE), collapse = "\n"),
                  error = function(e) stop("couldn't reach mygene.info (", conditionMessage(e),
                                           "); check the internet connection", call. = FALSE))
  value = jsonlite::fromJSON(text, simplifyVector = FALSE)
  if(!is.null(file)) writeLines(text, file)
  assign(key, value, envir = .gene_model_cache)
  value
}

# A short stable hash of a string, for cache file names
.url_hash = function(x){
  bytes = as.integer(charToRaw(x))
  h = 5381
  for(b in bytes) h = (h * 33 + b) %% 2147483647
  sprintf("%010.0f", h)
}

# Interval helpers (two-column matrices of start, end)
.interval_union = function(x){
  if(is.null(x) || nrow(x) == 0) return(matrix(numeric(0), 0, 2))
  x = x[order(x[, 1]), , drop = FALSE]
  group = cumsum(c(1, x[-1, 1] > cummax(x[, 2])[-nrow(x)]))
  cbind(as.numeric(tapply(x[, 1], group, min)), as.numeric(tapply(x[, 2], group, max)))
}
.interval_intersect = function(a, b){
  a = .interval_union(a); b = .interval_union(b)
  if(nrow(a) == 0 || nrow(b) == 0) return(matrix(numeric(0), 0, 2))
  out = do.call(rbind, lapply(seq_len(nrow(a)), function(i){
    lo = pmax(a[i, 1], b[, 1]); hi = pmin(a[i, 2], b[, 2])
    keep = hi > lo
    if(any(keep)) cbind(lo[keep], hi[keep]) else NULL
  }))
  if(is.null(out)) matrix(numeric(0), 0, 2) else out
}
.interval_subtract = function(a, b){
  a = .interval_union(a); b = .interval_union(b)
  if(nrow(a) == 0) return(a)
  if(nrow(b) == 0) return(a)
  out = do.call(rbind, lapply(seq_len(nrow(a)), function(i){
    pieces = matrix(a[i, ], 1)
    for(j in seq_len(nrow(b))){
      pieces = do.call(rbind, lapply(seq_len(nrow(pieces)), function(k){
        lo = pieces[k, 1]; hi = pieces[k, 2]
        if(b[j, 2] <= lo || b[j, 1] >= hi) return(matrix(c(lo, hi), 1))
        rbind(if(b[j, 1] > lo) c(lo, b[j, 1]), if(b[j, 2] < hi) c(b[j, 2], hi))
      }))
      if(is.null(pieces)) break
    }
    pieces
  }))
  if(is.null(out)) matrix(numeric(0), 0, 2) else matrix(out, ncol = 2)
}
