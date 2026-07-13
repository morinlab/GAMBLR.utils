#' @title Expand Gene Aliases.
#'
#' @description Given a vector of gene symbols, return the same symbols
#' plus any known alternate name for genes that have one, so symbol-based
#' filtering (\code{grep}, \code{dplyr::filter(Hugo_Symbol \%in\% ...)},
#' etc.) against a fixed gene list doesn't silently miss rows that are
#' annotated under the gene's other name.
#'
#' @details HGNC has renamed several genes over time (most notably the
#' entire histone gene family, e.g. \code{HIST1H1C} -> \code{H1-2}), and
#' different annotation pipelines/genome builds don't consistently agree on
#' which name they use. \code{aliases} (an internal hg19/hg38 lookup table)
#' records the old/new pairs currently known to affect GAMBL data; this
#' function expands a gene list to include both names for any gene present
#' in that table, checked in both directions, so it's safe to call
#' regardless of which nomenclature the input already uses.
#'
#' @param gene_symbol A vector of gene symbols.
#'
#' @return A character vector: the input symbols plus any known aliases for
#' genes among them, deduplicated. The original input order is preserved;
#' any aliases are appended at the end.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' expand_gene_aliases(c("HIST1H1C", "TP53"))
#'
expand_gene_aliases = function(gene_symbol){
  # Qualified even though this function lives in the same package: lazy-loaded
  # data/*.rda objects are not reliably bound in the package namespace itself
  # (only via the search path once the package is attached, or via an
  # explicit pkg::name lookup), so a bare `aliases` reference here can fail
  # with "object 'aliases' not found" depending on how GAMBLR.utils was
  # installed/loaded, even when GAMBLR.utils::aliases resolves fine.
  al = GAMBLR.utils::aliases
  alias_lookup = c(
    setNames(al$hg38, al$hg19),
    setNames(al$hg19, al$hg38)
  )
  extra = unname(alias_lookup[gene_symbol])
  extra = extra[!is.na(extra)]
  unique(c(gene_symbol, extra))
}
