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
#' \code{Hugo_Symbol} column.
#' @param genes A character vector of candidate gene symbols to check, e.g.
#' \code{expand_gene_aliases(my_genes)}.
#'
#' @return The subset of \code{genes} that appear at least once in
#' \code{maf_data$Hugo_Symbol}, in the same order as supplied.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' candidates <- expand_gene_aliases(c("C10orf12", "MYD88"))
#' locked_in <- sync_genes_to_maf(my_maf, candidates)
#' }
sync_genes_to_maf = function(maf_data, genes){
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
  genes[genes %in% present]
}
