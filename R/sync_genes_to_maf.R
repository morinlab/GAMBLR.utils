#' @title Sync Genes To MAF.
#'
#' @description Given a MAF-like data frame and a candidate gene list (which
#' may mix old/new HGNC names, e.g. from
#' \code{\link{expand_gene_aliases}}), return just the subset of names that
#' actually appear at least once as a \code{Hugo_Symbol} value in that MAF.
#' Lets a user lock in a gene list known to match their specific data,
#' rather than assuming a canonical/expanded gene list is fully present.
#'
#' @details Some HGNC renames (e.g. \code{C10orf12}/\code{LCOR}) don't have
#' a clean 1:1 relationship within a single MAF -- because symbol
#' annotation vintage is independent of genome build, a single MAF can
#' contain rows under both the old and new name for what is really one
#' locus. This function makes no attempt to consolidate such cases: if both
#' names are present, both are returned, since silently collapsing them
#' could hide the fact that your data labels the same gene inconsistently.
#' If you want a single, unique name for a gene like this, make that
#' decision explicitly by renaming the rows you consider equivalent, e.g.
#' \code{dplyr::mutate(maf_data, Hugo_Symbol = ifelse(Hugo_Symbol ==
#' "C10orf12", "LCOR", Hugo_Symbol))}, then re-run
#' \code{sync_genes_to_maf()} -- with only one of the two names left in
#' \code{maf_data}, only one will come back.
#'
#' @param maf_data A data frame (or \code{maf_data} object) with a
#' \code{Hugo_Symbol} column (and, if \code{summarize = TRUE}, a
#' \code{Variant_Classification} column).
#' @param genes A character vector of candidate gene symbols to check, e.g.
#' \code{expand_gene_aliases(my_genes)}.
#' @param summarize Default FALSE (return the character vector described
#' below). If TRUE, instead return a data frame with the number of coding
#' variants per gene -- useful for judging, at a glance, whether a gene
#' that's technically "present" actually has meaningful coding variant
#' support in this MAF, rather than just a handful of non-coding calls.
#' @param coding_classes Variant_Classification values counted as coding
#' when \code{summarize = TRUE}. Default \code{GAMBLR.helpers::coding_vc}.
#'
#' @return If \code{summarize = FALSE} (the default): the subset of
#' \code{genes} that appear at least once in \code{maf_data$Hugo_Symbol}, in
#' the same order as supplied. If \code{summarize = TRUE}: a data frame with
#' one row per gene in that subset and a column \code{n_coding} giving its
#' count of coding variants (0 if the gene is present but only via
#' non-coding rows), sorted descending by \code{n_coding}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' genome_meta = get_gambl_metadata() %>%
#'   filter(pathology=="FL",seq_type=="genome")
#' maf_grch37 = get_coding_ssm(these_samples_metadata = genome_meta, projection = "grch37")
#' candidates <- expand_gene_aliases(c("C10orf12", "MYD88", "FAKEGENE"))
#' found_genes <- sync_genes_to_maf(maf_grch37, candidates)
#' cat(found_genes)
#'
#' # same query genes, but as per-gene coding variant counts
#' sync_genes_to_maf(maf_grch37, candidates, summarize = TRUE)
#' }
sync_genes_to_maf = function(maf_data, genes, summarize = FALSE,
                              coding_classes = GAMBLR.helpers::coding_vc){
  if(!"Hugo_Symbol" %in% colnames(maf_data)){
    stop("maf_data must have a Hugo_Symbol column")
  }
  present = unique(maf_data$Hugo_Symbol)
  missing = setdiff(genes, present)
  if(length(missing) > 0){
    message(
      length(missing), " of ", length(genes),
      " requested gene(s) not found in maf_data$Hugo_Symbol: ",
      paste(missing, collapse = ", ")
    )
  }
  found = genes[genes %in% present]

  if(!summarize){
    return(found)
  }

  if(!"Variant_Classification" %in% colnames(maf_data)){
    stop("maf_data must have a Variant_Classification column when summarize = TRUE")
  }

  counts = maf_data %>%
    dplyr::filter(Hugo_Symbol %in% found, Variant_Classification %in% coding_classes) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::count(name = "n_coding") %>%
    dplyr::ungroup()

  # left_join against `found` (not just the counted genes) so a gene that's
  # present but has zero coding variants still shows up, with n_coding = 0,
  # instead of silently disappearing from the summary.
  data.frame(Hugo_Symbol = found) %>%
    dplyr::left_join(counts, by = "Hugo_Symbol") %>%
    dplyr::mutate(n_coding = dplyr::coalesce(n_coding, 0L)) %>%
    dplyr::arrange(dplyr::desc(n_coding))
}
