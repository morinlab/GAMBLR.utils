#' @title Drop mRNA-seq samples from sample metadata
#'
#' @description Remove the samples with seq_type "mrna" from a metadata
#' table, with a warning. Mutations are never associated with mRNA-seq
#' samples in GAMBL, so functions that summarise mutations per sample use
#' this to keep them out of their sample universe (where they would appear
#' unmutated).
#'
#' @param these_samples_metadata Metadata table, or NULL.
#' @param caller Name of the calling function, for the warning.
#'
#' @return `these_samples_metadata` without its mrna rows (unchanged when it
#' has none, has no seq_type column, or is NULL).
#'
#' @export
#'
#' @examples
#' meta = data.frame(sample_id = c("a", "b"), seq_type = c("genome", "mrna"))
#' drop_mrna_samples(meta)
drop_mrna_samples = function(these_samples_metadata, caller = NULL){
  if(is.null(these_samples_metadata) ||
     !"seq_type" %in% colnames(these_samples_metadata)){
    return(these_samples_metadata)
  }
  mrna = as.character(these_samples_metadata$seq_type) %in% "mrna"
  if(any(mrna)){
    warning(if(!is.null(caller)) paste0(caller, ": ") else "",
            "dropping ", sum(mrna), " mrna sample(s) from these_samples_metadata ",
            "(mutations are never associated with mRNA-seq samples)", call. = FALSE)
    these_samples_metadata = these_samples_metadata[!mrna, , drop = FALSE]
  }
  these_samples_metadata
}
