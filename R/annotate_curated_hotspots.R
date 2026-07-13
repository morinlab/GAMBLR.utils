 #' Assemble genetic features
#'
#' Build a sample-by-feature matrix from a MAF by annotating curated drivers,
#' summarizing coding/non-coding status, and collapsing per-gene columns.
#'
#' @param unannotated_maf MAF data frame to annotate and summarize.
#' @param genes_coding Character vector of genes to include for coding mutations.
#' @param genes_noncoding Character vector of genes to include for non-coding mutations.
#' @param genes_hotspot Character vector of genes with hotspot annotations (default: MYD88).
#' @param genes_driver Character vector of curated driver genes (default: CD79B, EZH2, NOTCH1, NOTCH2).
#' @param encoding_policy Named list of values used for coding/driver/noncoding/SV indicators.
#' @param these_samples_metadata Optional metadata for samples. If provided with
#' \code{sv_from_metadata}, SV status columns are appended via
#' \code{encode_sv_status()}.
#' @param sv_from_metadata Named vector or list mapping new SV column names to
#' metadata fields. If names do not end with \code{_SV}, that suffix is added.
#' @param sv_positive_values Values in metadata fields that indicate positive SV
#' status (default: \code{"POS"}).
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param return_everything Logical. If TRUE, return a list with
#' \code{annotated_maf}, \code{full_status}, and \code{collapsed_status}.
#'
#' @return A wide matrix-like data frame of per-sample mutation features.
#' If \code{return_everything = TRUE}, returns a list with intermediate outputs.
#'
#' @export

#'
#' @examples
#' \dontrun{
#' dlbcl_meta = get_gambl_metadata() %>%
#'   filter(seq_type != "mrna", pathology == "DLBCL")
#' # Get coding AND non-coding variants for all genes (including non-Tier 1)
#' # to allow for flexible summarization 
#' 
#' dlbcl_maf = get_all_coding_ssm(
#'   these_samples_metadata = dlbcl_meta,
#'   include_silent = TRUE
#' )
#' # use the Tier 1 DLBCL list (includes almost all LymphGen features and more)
#' dlbcl_genes = lymphoma_genes %>%
#'   filter(DLBCL_Tier == 1) %>%
#'   pull(Gene)
#' 
#' # genes affected by aSHM are not all in the Tier 1 DLBCL list
#' ashm_genes = unique(grch37_ashm_regions$gene)
#' ashm_genes = ashm_genes[ashm_genes %in% dlbcl_genes]
#'
#' assembled = assemble_genetic_features(
#'   dlbcl_maf,
#'   dlbcl_genes,
#'   ashm_genes,
#'   these_samples_metadata = dlbcl_meta
#' )
#' }
assemble_genetic_features <- function(unannotated_maf,
                                     genes_coding,
                                     genes_noncoding,
                                     genes_hotspot = c("MYD88"),
                                     genes_driver = c("CD79B","EZH2","NOTCH1","NOTCH2"),
                                     encoding_policy = list(coding=2,driver=2,noncoding=1,sv=2),
                                     these_samples_metadata = NULL,
                                     sv_from_metadata = c(BCL2="bcl2_ba", BCL6="bcl6_ba", MYC="myc_ba"),
                                     sv_positive_values = "POS",
                                     verbose = FALSE,
                                     return_everything = FALSE
                                     ){
  # to remove redundancy, we will drop all GENE_other for
  # the genes in genes_driver
  # to remove redundancy, we will drop all GENE_driver for
  # the genes in genes_hotspot
  #
  # Defaults are for LymphGen based on the description in this document
  # https://llmpp.ccr.cancer.gov/lymphgen/LymphGenInstructions.pdf?v=1771002307
  # EZH2: SET domain
  # 
  # 
  if(verbose){
    message("Annotating hotspots and drivers using annotate_curated_drivers...")
  }
  annotated_maf= annotate_curated_drivers(maf_data = unannotated_maf,
                                           genes_of_interest = unique(c(genes_hotspot,
                                                                        genes_driver)))
  genes_full = unique(c(genes_coding, genes_hotspot, genes_driver))
  if(verbose){
    message("Summarizing mutation status into a wide matrix with summarise_mutation_status...")
    message("encoding policy: coding=", encoding_policy$coding, ", driver=", encoding_policy$driver, ", noncoding=", encoding_policy$noncoding)
  }
  mut_status = summarise_mutation_status(
    annotated_maf,
    genes_coding = genes_coding,
    genes_noncoding = genes_noncoding,
    encoding_policy = encoding_policy
  )
  collapse_genes = setdiff(genes_full, genes_hotspot)
  collapse_genes = setdiff(collapse_genes, genes_driver)
  
  
  status_collapse  = collapse_columns_by_regex(
    mut_status,
    genes = collapse_genes,
    suffix = ""
  )
  drop_other = paste0(genes_driver,"_other") 
  status_collapse = status_collapse %>% dplyr::select(-any_of(drop_other))
  drop_driver = paste0(genes_hotspot,"_driver") 
  status_collapse = status_collapse %>% dplyr::select(-any_of(drop_driver))

  driver_cols <- grep("_driver$", colnames(status_collapse), value = TRUE)
  if (length(driver_cols) > 0) {
    new_names <- sub("_driver$", "", driver_cols)
    status_collapse <- status_collapse %>%
      dplyr::rename(!!!setNames(driver_cols, new_names))
  }

  if (!is.null(these_samples_metadata) && !is.null(sv_from_metadata)) {
    if (any(grepl("_SV$", names(sv_from_metadata)))) {
      mapping <- as.list(sv_from_metadata)
    } else {
      mapping <- setNames(as.list(sv_from_metadata), paste0(names(sv_from_metadata), "_SV"))
    }
    sv_mat <- encode_sv_status(
      mat = status_collapse,
      metadata = these_samples_metadata,
      mapping = mapping,
      positive_values = sv_positive_values,
      encoding_policy = list(positive = encoding_policy$sv, negative = 0),
      verbose = verbose
    )
    for (col_name in names(mapping)) {
      if (col_name %in% colnames(sv_mat) && all(sv_mat[[col_name]] == 0, na.rm = TRUE)) {
        warning(paste("SV encoding: no POS values found for", col_name))
      }
    }
    status_collapse <- sv_mat
  }
  if(return_everything){
    return(list(
      annotated_maf = annotated_maf,
      full_status = mut_status,
      collapsed_status = status_collapse
    ))
  }
  return(status_collapse)
  
}



#' @title Summarize mutation status
#'
#' @description Build a wide sample-by-feature binary (or weighted) matrix from
#' an annotated MAF. Each column represents one feature — either a gene-level
#' summary or a specific mutation alias — and each cell contains a numeric
#' indicator defined by \code{encoding_policy} (0 = absent, non-zero = present).
#'
#' @details
#' ## Input requirements
#' The input must be the output of \code{\link{annotate_curated_drivers}}, which
#' adds the \code{driver_alias}, \code{mutation_alias}, \code{mutation_alias_all},
#' and \code{hotspot_alias} columns used here.
#'
#' ## Column naming conventions
#' The function produces columns named \code{GENE_suffix}, where the suffix
#' depends on the gene's reporting policy:
#'
#' \describe{
#'   \item{\code{GENE_coding}}{Any mutation present in the gene (the base
#'     indicator). Produced for genes with no driver annotations, genes under
#'     the \code{coding} policy, genes under the \code{detailed} policy that
#'     lack alias annotations (fallback), and genes under the \code{maximal}
#'     policy.}
#'   \item{\code{GENE_driver} / \code{GENE_other}}{Produced for genes that have
#'     driver alias annotations (e.g. \code{MYD88_L265P}, \code{MYD88_other}).
#'     \code{_driver} columns reflect the specific annotated driver alias;
#'     \code{_other} is the catch-all for unannotated mutations.}
#'   \item{Per-alias columns (e.g. \code{EZH2_Y641}, \code{EZH2_other})}{
#'     Produced for genes under the \code{detailed} or \code{maximal} policy
#'     when \code{mutation_alias_all} or \code{hotspot_alias} entries exist.}
#'   \item{\code{GENE_noncoding}}{Produced when the gene appears in
#'     \code{genes_noncoding}. Counts mutations with
#'     \code{Variant_Classification} NOT in \code{vc_nonSynonymous} (i.e.
#'     Silent, UTR, Flank). The input MAF must include these classes.}
#'   \item{\code{GENE_any}}{Produced only for genes under the \code{maximal}
#'     policy. The element-wise maximum across all \code{GENE_*} columns for
#'     that gene, including \code{GENE_noncoding} when the gene is also in
#'     \code{genes_noncoding}.}
#' }
#'
#' ## Default reporting behaviour
#' Genes with driver alias annotations (i.e. \code{driver_alias} contains a
#' \code{_driver} suffix entry) produce \code{GENE_driver}/\code{GENE_other}
#' columns and do \strong{not} get a \code{GENE_coding} column. Genes without
#' driver annotations produce only \code{GENE_coding}. Use
#' \code{gene_reporting_policy} to override this per gene.
#'
#' ## Policy interactions with \code{genes_coding}
#' Genes listed under any \code{gene_reporting_policy} key are automatically
#' added to the internal coding gene universe regardless of whether they appear
#' in \code{genes_coding}. This means you do not need to manually add policy
#' genes to \code{genes_coding}.
#'
#' @param annotated_maf A data frame in MAF format produced by
#'   \code{\link{annotate_curated_drivers}}. Must include \code{Hugo_Symbol},
#'   \code{Tumor_Sample_Barcode}, \code{Variant_Classification},
#'   \code{driver_alias}, \code{mutation_alias}, \code{mutation_alias_all},
#'   and \code{hotspot_alias}.
#' @param genes_coding Optional character vector of genes for which to produce
#'   \code{GENE_coding} or driver-alias columns. Defaults to all Tier 1
#'   lymphoma genes. Genes from \code{gene_reporting_policy} are unioned into
#'   this set automatically.
#' @param genes_noncoding Optional character vector of genes for which to also
#'   count non-coding mutations. Produces a \code{GENE_noncoding} column.
#'   When \code{noncoding_maf} is \code{NULL}, non-coding mutations are taken
#'   from \code{annotated_maf} by excluding rows whose
#'   \code{Variant_Classification} is in \code{vc_nonSynonymous} (Silent, UTR,
#'   Flank), so the input MAF must include those classes. When
#'   \code{noncoding_maf} is provided, all rows in that MAF for genes in
#'   \code{genes_noncoding} are counted regardless of
#'   \code{Variant_Classification}. Default \code{NULL} (no non-coding
#'   columns).
#' @param noncoding_maf Optional MAF data frame supplying non-coding mutations
#'   separately from \code{annotated_maf}. Must contain at least
#'   \code{Hugo_Symbol} and \code{Tumor_Sample_Barcode}. Intended for
#'   region-restricted non-coding mutation sets (e.g. aSHM hotspot windows)
#'   where extracting non-coding variants for all genes from a single MAF
#'   would be inefficient. When provided, every row whose \code{Hugo_Symbol}
#'   is in \code{genes_noncoding} contributes to the \code{GENE_noncoding}
#'   count; no \code{Variant_Classification} filter is applied. Ignored if
#'   \code{genes_noncoding} is \code{NULL}. Default \code{NULL}.
#' @param gene_reporting_policy Optional named list controlling per-gene column
#'   output. Valid keys:
#'   \describe{
#'     \item{\code{detailed}}{Produce per-alias columns (from
#'       \code{mutation_alias_all} and \code{hotspot_alias}) plus a
#'       \code{GENE_other} catch-all for mutations that match no alias.
#'       Genes with no alias data at all fall back to \code{_coding} with a
#'       warning.}
#'     \item{\code{coding}}{Produce only \code{GENE_coding}, suppressing
#'       driver-alias columns. Useful to coarsen a gene that has driver
#'       annotations but where you do not want alias-level resolution.}
#'     \item{\code{maximal}}{Produce the full column set: \code{GENE_coding},
#'       per-alias columns, \code{GENE_driver}/\code{GENE_other}, and
#'       \code{GENE_any}. If the gene is also in \code{genes_noncoding},
#'       \code{GENE_noncoding} is included in \code{GENE_any}.}
#'   }
#'   A gene may appear under only one key. Default \code{NULL} (all genes use
#'   the default driver/coding behaviour).
#' @param genes_drop_unannotated Optional character vector of genes for which
#'   the catch-all \code{GENE_other} column should be removed from the output.
#'   Shorthand for \code{gene_drop_policy = list(unannotated = ...)}.
#'   Mutually exclusive with \code{gene_drop_policy}.
#' @param gene_drop_policy Optional named list with keys \code{all} and/or
#'   \code{unannotated}. \code{all}: remove every \code{GENE_*} column for the
#'   listed genes. \code{unannotated}: remove only the \code{GENE_other}
#'   catch-all column. Mutually exclusive with \code{genes_drop_unannotated}.
#' @param genes_use_counts Optional character vector of gene symbols for which
#'   every \code{GENE_*} column reports the raw number of distinct mutation
#'   occurrences observed in that sample, instead of the fixed
#'   \code{encoding_policy} value. For example, if a sample has 3 separate
#'   \code{KMT2D_trunc} mutations, \code{KMT2D_trunc} is \code{3} for that
#'   sample rather than \code{encoding_policy$driver}. Genes not listed here
#'   are unaffected and continue to use \code{encoding_policy} as usual. Under
#'   the \code{maximal} policy, \code{GENE_any} is still the element-wise
#'   maximum across that gene's columns, so it becomes the largest single
#'   per-feature count rather than a binary indicator. Default \code{NULL}
#'   (all genes use \code{encoding_policy}).
#' @param encoding_policy Named list of numeric values for each indicator type.
#'   Default \code{list(coding = 2, noncoding = 1, driver = 2)}. Setting all
#'   values to 1 produces a standard binary 0/1 matrix.
#' @param use_all_aliases Logical. If \code{TRUE} (default), per-alias columns
#'   are built from \code{mutation_alias_all} (semicolon-separated specific
#'   aliases, no catch-alls). If \code{FALSE}, uses \code{mutation_alias}
#'   (which may contain \code{GENE_other} catch-all values and produce fewer
#'   specific alias columns).
#' @param include_hotspot_alias Logical. If \code{TRUE} (default), also builds
#'   alias columns from \code{hotspot_alias} for genes under the \code{detailed}
#'   or \code{maximal} policy, and for all other annotated genes via a separate
#'   \code{hotspot_count} table. Default \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, prints messages about which policy
#'   and drop decisions are being applied. Default \code{FALSE}.
#' @param these_samples_metadata Optional data frame with a \code{sample_id}
#'   column defining the complete sample universe. When supplied: (1) any
#'   \code{Tumor_Sample_Barcode} in the MAF that is absent from
#'   \code{sample_id} is dropped with a warning; (2) every \code{sample_id}
#'   present in the metadata receives a row in the output, filled with 0 for
#'   all features if no mutations were observed. Default \code{NULL} (sample
#'   universe derived from the MAF).
#'
#' @return A data frame with one row per sample (row names = sample IDs) and
#'   one column per feature. Cell values are 0 (absent) or a positive integer
#'   set by \code{encoding_policy}. Column types that may appear:
#'   \code{GENE_coding}, \code{GENE_noncoding}, \code{GENE_driver},
#'   \code{GENE_other}, per-alias columns (e.g. \code{EZH2_Y641}),
#'   \code{GENE_any} (maximal policy only).
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage: default Tier 1 genes, driver-alias columns where available
#' maf = GAMBLR.open::get_coding_ssm(
#'   these_samples_metadata = GAMBLR.open::get_gambl_metadata(),
#'   this_seq_type = "capture"
#' )
#' maf_anno = annotate_curated_drivers(maf)
#' mat = summarise_mutation_status(maf_anno)
#' # Each row is a sample; columns like MYD88_L265P, MYD88_other, TP53_coding
#' mat[1:5, 1:8]
#' }
#'
#' \dontrun{
#' # Per-gene reporting policies:
#' #   EZH2 -> maximal: gets _coding, per-alias, _driver, _other, _any columns
#' #   MYD88 -> detailed: gets only per-alias columns (MYD88_L265P, etc.)
#' #   NOTCH1 -> coding: collapses everything to a single NOTCH1_coding column
#' # Genes not listed (TP53, KMT2D, ...) use default driver/coding behaviour.
#' # Drop the unannotated catch-all for KMT2D to keep only annotated drivers.
#' maf_anno = annotate_curated_drivers(GAMBLR.open::get_coding_ssm())
#' mat = summarise_mutation_status(
#'   maf_anno,
#'   gene_reporting_policy = list(
#'     maximal  = c("EZH2"),
#'     detailed = c("MYD88"),
#'     coding   = c("NOTCH1")
#'   ),
#'   gene_drop_policy = list(unannotated = c("KMT2D"))
#' )
#' grep("EZH2", colnames(mat), value = TRUE)
#' # [1] "EZH2_coding" "EZH2_Y641"  "EZH2_A677"  "EZH2_driver" "EZH2_other" "EZH2_any"
#' }
#'
#' \dontrun{
#' # Non-coding mutations for aSHM target genes.
#' # The input MAF must include non-coding variant classes (Silent, UTR, Flank).
#' maf_full = GAMBLR.open::get_ssm_by_samples(
#'   these_samples_metadata = GAMBLR.open::get_gambl_metadata(),
#'   coding_only = FALSE
#' )
#' maf_anno = annotate_curated_drivers(maf_full)
#' mat = summarise_mutation_status(
#'   maf_anno,
#'   genes_noncoding = c("PIM1", "SOCS1", "BCL6"),
#'   # use binary encoding for all indicator types
#'   encoding_policy = list(coding = 1, noncoding = 1, driver = 1)
#' )
#' # PIM1_coding: any coding mutation in PIM1
#' # PIM1_noncoding: Silent/UTR/Flank mutations in PIM1
#' grep("PIM1", colnames(mat), value = TRUE)
#' }
#'
#' \dontrun{
#' # Report raw occurrence counts for KMT2D instead of the encoding_policy
#' # value. A sample with 3 separate KMT2D_trunc mutations gets 3, not 2.
#' # All other genes are unaffected and still use encoding_policy.
#' mat = summarise_mutation_status(
#'   maf_anno,
#'   genes_use_counts = c("KMT2D")
#' )
#' }
#'
summarise_mutation_status = function(annotated_maf,
                                     genes_coding = NULL,
                                     genes_noncoding = NULL,
                                     noncoding_maf = NULL,
                                     gene_reporting_policy = NULL,
                                     genes_drop_unannotated = NULL,
                                     gene_drop_policy = NULL,
                                     encoding_policy = list(coding = 2, noncoding = 1, driver = 2),
                                     genes_use_counts = NULL,
                                     use_all_aliases = TRUE,
                                     include_hotspot_alias = TRUE,
                                     verbose = FALSE,
                                     these_samples_metadata = NULL){

  # per-row "mutated" value: the raw occurrence count for genes listed in
  # genes_use_counts, otherwise the fixed encoding_policy value. Gene is
  # recovered from the variant name's prefix (e.g. "KMT2D_trunc" -> "KMT2D"),
  # which holds for both GENE_<suffix> columns and alias strings alike.
  resolve_mutated_value = function(variant, n, policy_value, genes_use_counts){
    if(is.null(genes_use_counts) || length(genes_use_counts) == 0){
      return(rep(policy_value, length(variant)))
    }
    gene = sub("_.*$", "", variant)
    ifelse(gene %in% genes_use_counts, n, policy_value)
  }

  if(is.null(genes_coding)){
    genes_coding = lymphoma_genes %>%
    dplyr::filter(FL_Tier==1 | MCL_Tier== 1 | BL_Tier==1 | DLBCL_Tier==1) %>%
    pull(Gene) %>%
    expand_gene_aliases()
  }

  # TODO similarly fill in counts for aSHM genes (should restrict to aSHM coordinates)

  if(is.null(gene_reporting_policy)){
    gene_reporting_policy = list()
  }
  detailed_genes = gene_reporting_policy$detailed
  if(is.null(detailed_genes)){
    detailed_genes = character(0)
  }
  coding_policy_genes = gene_reporting_policy$coding
  if(is.null(coding_policy_genes)){
    coding_policy_genes = character(0)
  }
  maximal_genes = gene_reporting_policy$maximal
  if(is.null(maximal_genes)){
    maximal_genes = character(0)
  }
  overlap = intersect(detailed_genes, coding_policy_genes)
  if(length(overlap) > 0){
    stop(paste0("Genes cannot be in both gene_reporting_policy$detailed and $coding: ", paste(overlap, collapse=", ")))
  }
  overlap = intersect(maximal_genes, union(detailed_genes, coding_policy_genes))
  if(length(overlap) > 0){
    stop(paste0("Genes cannot be in gene_reporting_policy$maximal and another policy: ", paste(overlap, collapse=", ")))
  }

  if(verbose){
    if(length(detailed_genes) > 0){
      message("Will keep detailed (alias-level) columns for genes: ", paste(detailed_genes, collapse = ", "))
      message("Will drop _driver/_other columns for genes: ", paste(detailed_genes, collapse = ", "))
    }
    if(length(coding_policy_genes) > 0){
      message("Will keep _coding columns (and drop driver/other and detailed) for genes: ", paste(coding_policy_genes, collapse = ", "))
    }
    if(length(maximal_genes) > 0){
      message("Will produce maximal columns (_coding, _driver, _other, aliases, _any) for genes: ", paste(maximal_genes, collapse = ", "))
    }
  }


  if(!is.null(genes_drop_unannotated) && !is.null(gene_drop_policy)){
    stop("Provide only one of genes_drop_unannotated or gene_drop_policy")
  }
  if(is.null(gene_drop_policy) && !is.null(genes_drop_unannotated)){
    gene_drop_policy = list(unannotated = genes_drop_unannotated)
  }
  if(is.null(gene_drop_policy)){
    gene_drop_policy = list()
  }
  if(!is.null(gene_drop_policy$all) && !is.null(gene_drop_policy$unannotated)){
    overlap = intersect(gene_drop_policy$all, gene_drop_policy$unannotated)
    if(length(overlap) > 0){
      stop(paste0("Genes cannot be in both gene_drop_policy$all and $unannotated: ", paste(overlap, collapse=", ")))
    }
  }


  # ── sample universe from metadata ─────────────────────────────────────────
  meta_sample_ids = NULL
  if (!is.null(these_samples_metadata)) {
    if (!"sample_id" %in% colnames(these_samples_metadata)) {
      stop("these_samples_metadata must contain a 'sample_id' column.")
    }
    meta_sample_ids = these_samples_metadata$sample_id
    extra = setdiff(unique(annotated_maf$Tumor_Sample_Barcode), meta_sample_ids)
    if (length(extra) > 0) {
      warning(
        length(extra), " Tumor_Sample_Barcode(s) in MAF absent from metadata; dropped."
      )
      annotated_maf = dplyr::filter(annotated_maf, Tumor_Sample_Barcode %in% meta_sample_ids)
    }
  }

  #expects the output of Annotate curated drivers
  if(!is.null(genes_noncoding) && length(genes_noncoding) > 0){

    noncoding_source = if (!is.null(noncoding_maf)) {
      dplyr::filter(noncoding_maf, Hugo_Symbol %in% genes_noncoding)
    } else {
      dplyr::filter(annotated_maf,
                    !Variant_Classification %in% vc_nonSynonymous,
                    Hugo_Symbol %in% genes_noncoding)
    }

    noncoding_count = noncoding_source %>%
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      dplyr::count(Hugo_Symbol, Tumor_Sample_Barcode, name = "n") %>%
      mutate(Variant=paste0(Hugo_Symbol,"_noncoding")) %>%
      dplyr::select(-Hugo_Symbol) %>%
      mutate(mutated = resolve_mutated_value(Variant, n, encoding_policy$noncoding, genes_use_counts)) %>%
      dplyr::select(-n) # will fill missing with 0 at the join stage

  }

  coding_maf = annotated_maf
  #coding_maf = dplyr::filter(annotated_maf, Variant_Classification %in% vc_nonSynonymous)

  #separately count:
  # - coding variants
  # - drivers
  # - non-drivers
  # - each distinct hotspot

  # genes under any reporting policy must be in coding_maf so that detailed-policy
  # genes with no aliases can fall back to _coding and coding-policy genes are
  # always represented
  genes_coding = union(genes_coding,
                       union(detailed_genes, union(coding_policy_genes, maximal_genes)))

  coding_maf = coding_maf %>%
    dplyr::filter(Hugo_Symbol %in% genes_coding) %>%
    strip_genomic_classes()

  

  # detailed_maf: detailed_genes only (used for skip_genes_coding fallback check)
  detailed_maf = coding_maf %>%
    dplyr::filter(Hugo_Symbol %in% detailed_genes)

  # alias_maf: detailed + maximal genes both get per-alias columns
  alias_genes = union(detailed_genes, maximal_genes)
  alias_maf = coding_maf %>%
    dplyr::filter(Hugo_Symbol %in% alias_genes)

  if(use_all_aliases && "mutation_alias_all" %in% colnames(alias_maf)){
    detailed_alias_long = alias_maf %>%
      dplyr::select(mutation_alias_all, Tumor_Sample_Barcode) %>%
      tidyr::separate_rows(mutation_alias_all, sep = ";") %>%
      dplyr::mutate(mutation_alias_all = trimws(mutation_alias_all)) %>%
      dplyr::filter(!is.na(mutation_alias_all), mutation_alias_all != "") %>%
      dplyr::rename(mutation_alias = mutation_alias_all)
  }else{
    detailed_alias_long = alias_maf %>%
      dplyr::select(mutation_alias, Tumor_Sample_Barcode) %>%
      dplyr::filter(!is.na(mutation_alias), mutation_alias != "")
  }

  if(include_hotspot_alias && "hotspot_alias" %in% colnames(alias_maf)){
    hotspot_alias_long = alias_maf %>%
      dplyr::select(hotspot_alias, Tumor_Sample_Barcode) %>%
      dplyr::filter(!is.na(hotspot_alias), hotspot_alias != "") %>%
      dplyr::rename(mutation_alias = hotspot_alias)
    detailed_alias_long = dplyr::bind_rows(detailed_alias_long, hotspot_alias_long)
  }

  detailed_driver_count = detailed_alias_long %>%
    dplyr::count(mutation_alias, Tumor_Sample_Barcode, name = "n") %>%
    rename(Variant=mutation_alias) %>%
    mutate(mutated = resolve_mutated_value(Variant, n, encoding_policy$driver, genes_use_counts)) %>%
    dplyr::select(-n) # will fill missing with 0 at the join stage

  driver_report_genes = dplyr::filter(coding_maf, grepl("_driver", driver_alias)) %>%
    pull(Hugo_Symbol) %>%
    unique()
  if(length(detailed_genes) > 0){
    driver_report_genes = setdiff(driver_report_genes, detailed_genes)
  }
  if(length(coding_policy_genes) > 0){
    driver_report_genes = setdiff(driver_report_genes, coding_policy_genes)
  }
  if(length(maximal_genes) > 0){
    driver_report_genes = setdiff(driver_report_genes, maximal_genes)
  }
  driver_count = coding_maf %>%
    dplyr::filter(Hugo_Symbol %in% driver_report_genes) %>%
    dplyr::select(driver_alias, Tumor_Sample_Barcode) %>%
    dplyr::count(driver_alias, Tumor_Sample_Barcode, name = "n") %>%
    rename(Variant=driver_alias) %>%
    mutate(mutated = resolve_mutated_value(Variant, n, encoding_policy$driver, genes_use_counts)) %>%
    dplyr::select(-n) # will fill missing with 0 at the join stage

  # maximal genes: count all driver_alias values (_driver and _other)
  maximal_driver_count = NULL
  if(length(maximal_genes) > 0){
    maximal_driver_count = coding_maf %>%
      dplyr::filter(Hugo_Symbol %in% maximal_genes,
                    !is.na(driver_alias), driver_alias != "") %>%
      dplyr::select(driver_alias, Tumor_Sample_Barcode) %>%
      dplyr::count(driver_alias, Tumor_Sample_Barcode, name = "n") %>%
      rename(Variant = driver_alias) %>%
      mutate(mutated = resolve_mutated_value(Variant, n, encoding_policy$driver, genes_use_counts)) %>%
      dplyr::select(-n)
  }

  hotspot_count = NULL
  if(include_hotspot_alias && "hotspot_alias" %in% colnames(annotated_maf)){
    hotspot_count = annotated_maf %>%
      dplyr::filter(!is.na(hotspot_alias), hotspot_alias != "",
                    !Hugo_Symbol %in% detailed_genes,
                    !Hugo_Symbol %in% maximal_genes) %>%
      dplyr::select(hotspot_alias, Tumor_Sample_Barcode) %>%
      dplyr::count(hotspot_alias, Tumor_Sample_Barcode, name = "n") %>%
      rename(Variant=hotspot_alias) %>%
      mutate(mutated = resolve_mutated_value(Variant, n, encoding_policy$driver, genes_use_counts)) %>% # use driver encoding for hotspot
      dplyr::select(-n)
  }

  # only exclude detailed_genes that actually have alias coverage; genes with no
  # mutation_alias or hotspot_alias in the data fall back to _coding reporting.
  # use the same column that detailed_alias_long was built from: mutation_alias_all
  # when use_all_aliases is TRUE (mutation_alias may contain catch-all "_other"
  # values that produce no rows in detailed_alias_long)
  if (use_all_aliases && "mutation_alias_all" %in% colnames(detailed_maf)) {
    has_mutation_alias <- !is.na(detailed_maf$mutation_alias_all) &
      detailed_maf$mutation_alias_all != ""
  } else {
    has_mutation_alias <- !is.na(detailed_maf$mutation_alias) &
      detailed_maf$mutation_alias != ""
  }
  if (include_hotspot_alias && "hotspot_alias" %in% colnames(detailed_maf)) {
    has_hotspot_alias <- !is.na(detailed_maf$hotspot_alias) &
      detailed_maf$hotspot_alias != ""
  } else {
    has_hotspot_alias <- logical(nrow(detailed_maf))
  }
  detailed_genes_with_aliases <- unique(
    detailed_maf$Hugo_Symbol[has_mutation_alias | has_hotspot_alias]
  )
  detailed_genes_fallback <- setdiff(detailed_genes, detailed_genes_with_aliases)
  if (length(detailed_genes_fallback) > 0) {
    warning(
      "No aliases found for detailed genes; falling back to _coding for: ",
      paste(detailed_genes_fallback, collapse = ", ")
    )
  }
  skip_genes_coding <- union(driver_report_genes, detailed_genes_with_aliases)

  # catch-all for detailed genes: mutations with no alias coverage become GENE_other.
  # Without this, samples with unannotated mutations in a detailed gene are silently
  # dropped (they appear 0 in all alias columns and don't get a _coding fallback).
  detailed_other_count = NULL
  if (length(detailed_genes_with_aliases) > 0) {
    detailed_sub = alias_maf %>%
      dplyr::filter(Hugo_Symbol %in% detailed_genes_with_aliases)

    if (use_all_aliases && "mutation_alias_all" %in% colnames(detailed_sub)) {
      no_alias = is.na(detailed_sub$mutation_alias_all) | detailed_sub$mutation_alias_all == ""
    } else {
      no_alias = is.na(detailed_sub$mutation_alias) |
                 detailed_sub$mutation_alias == "" |
                 grepl("_other$", detailed_sub$mutation_alias)
    }
    if (include_hotspot_alias && "hotspot_alias" %in% colnames(detailed_sub)) {
      has_hotspot = !is.na(detailed_sub$hotspot_alias) & detailed_sub$hotspot_alias != ""
      no_alias = no_alias & !has_hotspot
    }

    unaliased = detailed_sub[no_alias, ]
    if (nrow(unaliased) > 0) {
      detailed_other_count = unaliased %>%
        dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
        dplyr::count(Hugo_Symbol, Tumor_Sample_Barcode, name = "n") %>%
        dplyr::mutate(Variant = paste0(Hugo_Symbol, "_other")) %>%
        dplyr::select(-Hugo_Symbol) %>%
        dplyr::mutate(mutated = resolve_mutated_value(Variant, n, encoding_policy$coding, genes_use_counts)) %>%
        dplyr::select(-n)
    }
  }

  coding_count = coding_maf %>%
    dplyr::filter(!Hugo_Symbol %in% skip_genes_coding) %>%
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    dplyr::count(Hugo_Symbol, Tumor_Sample_Barcode, name = "n") %>%
    mutate(Variant=paste0(Hugo_Symbol,"_coding")) %>%
    dplyr::select(-Hugo_Symbol) %>%
    mutate(mutated = resolve_mutated_value(Variant, n, encoding_policy$coding, genes_use_counts)) %>%
    dplyr::select(-n) # will fill missing with 0 at the join stage
  
  if(!is.null(genes_noncoding) && length(genes_noncoding) > 0){
    long = bind_rows(
      noncoding_count,
      coding_count,
      detailed_driver_count,
      detailed_other_count,
      maximal_driver_count,
      driver_count,
      hotspot_count
    )
  }else{
    long = bind_rows(
      coding_count,
      detailed_driver_count,
      detailed_other_count,
      maximal_driver_count,
      driver_count,
      hotspot_count
    )
  }


  all_samples = if (!is.null(meta_sample_ids)) {
    tibble::tibble(Tumor_Sample_Barcode = meta_sample_ids)
  } else {
    dplyr::distinct(coding_maf, Tumor_Sample_Barcode)
  }
  all_variants = distinct(long, Variant)

  wide = long %>%
    tidyr::complete(all_samples, all_variants, fill = list(mutated = 0)) %>%
    tidyr::pivot_wider(names_from = Variant, values_from = mutated, values_fill = 0) %>%
    column_to_rownames("Tumor_Sample_Barcode")

  # maximal policy: collapse all GENE_* columns into GENE_any; _noncoding is
  # included automatically if the gene is also in genes_noncoding (it then has
  # a GENE_noncoding column that matches the regex), otherwise it is excluded
  if(length(maximal_genes) > 0){
    maximal_present = intersect(
      maximal_genes,
      unique(sub("_[^_]+$", "", colnames(wide)))
    )
    if(length(maximal_present) > 0){
      wide = collapse_columns_by_regex(
        wide,
        genes        = maximal_present,
        suffix       = "_any",
        drop_sources = FALSE,
        agg_fn       = max
      )
    }
  }

  if(length(gene_drop_policy) > 0){
    drop_cols = character(0)
    if(!is.null(gene_drop_policy$unannotated)){
      drop_cols = c(drop_cols, unlist(lapply(gene_drop_policy$unannotated, function(g){
        c(paste0(g, "_other"))
      })))
      if(verbose && length(gene_drop_policy$unannotated) > 0){
        message("Will drop _other columns for genes: ", paste(gene_drop_policy$unannotated, collapse = ", "))
      }
    }
    if(!is.null(gene_drop_policy$all)){
      drop_cols = c(drop_cols, unlist(lapply(gene_drop_policy$all, function(g){
        grep(paste0("^", g, "_"), colnames(wide), value = TRUE)
      })))
      if(verbose && length(gene_drop_policy$all) > 0){
        message("Will drop all columns for genes: ", paste(gene_drop_policy$all, collapse = ", "))
      }
    }
    if(length(drop_cols) > 0){
      wide = wide[, setdiff(colnames(wide), unique(drop_cols)), drop = FALSE]
    }
  }

  return(wide)

}



#' @title Add binary status columns from metadata
#'
#' @description Add one or more indicator columns to a matrix/data frame based on
#' metadata fields (e.g., BCL2/BCL6/MYC rearrangement status).
#'
#' @param mat A matrix or data frame with sample IDs as row names.
#' @param metadata A data frame containing sample-level metadata.
#' @param mapping Named character vector or list mapping new column names to
#' metadata field names. For example, \code{c(BCL6_SV = "bcl6_ba")}.
#' @param positive_values Vector of values in the metadata field that indicate
#' a positive status (default \code{"POS"}). Any non-matching or missing values
#' are encoded as \code{encoding_policy$negative}.
#' @param encoding_policy Named list controlling the numeric values used for
#' positive and negative samples. Default is
#' \code{list(positive = 2, negative = 0)}.
#' @param sample_id_col Column name in \code{metadata} containing sample IDs
#' (default \code{"sample_id"}).
#' @param verbose Logical. If TRUE, prints matching diagnostics.
#'
#' @return A data.frame with the input \code{mat} columns and new SV columns
#' appended. All non-positive samples are encoded as
#' \code{encoding_policy$negative}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dlbcl_mat = encode_sv_status(
#'   dlbcl_mat,
#'   dlbcl_meta,
#'   mapping = c(BCL6_SV = "bcl6_ba", BCL2_SV = "bcl2_ba", MYC_SV = "myc_ba"),
#'   positive_values = "POS",
#'   encoding_policy = list(positive = 2, negative = 0)
#' )
#' }
#'
encode_sv_status = function(mat,
                                            metadata,
                                            mapping = c(BCL6_SV = "bcl6_ba", BCL2_SV = "bcl2_ba", MYC_SV = "myc_ba"),
                                            positive_values = "POS",
                                            encoding_policy = list(positive = 2, negative = 0),
                                            sample_id_col = "sample_id",
                                            verbose = FALSE){
  if(is.null(rownames(mat))){
    stop("mat must have row names corresponding to sample IDs")
  }
  if(is.null(mapping) || length(mapping) == 0){
    stop("mapping must be a named vector or list")
  }
  if(is.list(mapping) && is.null(names(mapping))){
    stop("mapping must be a named vector or named list")
  }
  if(!sample_id_col %in% colnames(metadata)){
    stop(paste0("metadata must contain column: ", sample_id_col))
  }

  mat_df = as.data.frame(mat)
  sample_ids = rownames(mat_df)

  if(verbose){
    message("encode_sv_status: rows in mat = ", length(sample_ids))
    message("encode_sv_status: rows in metadata = ", nrow(metadata))
    message("encode_sv_status: unique sample ids in metadata = ", length(unique(metadata[[sample_id_col]])))
    message("encode_sv_status: overlap with mat rownames = ", length(base::intersect(sample_ids, metadata[[sample_id_col]])))
  }

  for(new_col in names(mapping)){
    field = mapping[[new_col]]
    if(!field %in% colnames(metadata)){
      stop(paste0("metadata missing field: ", field))
    }
    pos_ids = metadata[[sample_id_col]][metadata[[field]] %in% positive_values]
    pos_ids = as.character(pos_ids)
    pos_ids = base::intersect(pos_ids, as.character(sample_ids))
    if(verbose){
      message("encode_sv_status: ", new_col, " field=", field,
              " positives in metadata=", sum(metadata[[field]] %in% positive_values, na.rm = TRUE),
              ", matched=", length(pos_ids))
    }
    mat_df[[new_col]] = encoding_policy$negative
    mat_df[pos_ids, new_col] = encoding_policy$positive
  }

  return(mat_df)
}


#' @title Collapse columns by regex
#'
#' @description Collapse multiple columns into a single column using regex matches.
#'
#' @details This helper creates new columns based on regex patterns and computes
#' an aggregate value (default: row-wise max) across matched source columns.
#'
#' @param x A data frame or matrix with columns to collapse.
#' @param collapse_policy A named list where each name is the new column name and
#' each value is a character vector of regex patterns used to select source columns.
#' @param genes Optional character vector of gene symbols to collapse using a
#' default pattern of \code{^GENE_}. If provided, \code{collapse_policy} is
#' generated automatically unless explicitly supplied.
#' @param suffix Suffix for auto-generated columns when \code{genes} is used
#' (default \code{"_any"}, producing columns like \code{GENE_any}).
#' @param drop_sources Logical. If TRUE (default), matched source columns are dropped
#' after the new column is created.
#' @param agg_fn Function to aggregate matched columns per row. Defaults to max.
#'
#' @return A data frame with collapsed columns added (and optionally source columns removed).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # status_df created with a call to summarise_mutation_status
#' grep("MYC",colnames(status_df),value=T)
#' grep("BCL2",colnames(status_df),value=T)
#' 
#' # each of these genes has multiple columns such as MYC_SV, MYC_noncoding
#' # that we want combined into one unified column per gene
#' status_df = collapse_columns_by_regex(
#'   status_df,
#'   collapse_policy = list(MYC_any = c("^MYC_"), BCL2_any = c("^BCL2_"))
#' )
#' }
#' # status_df will have a MYC_any and BCL2_any column and the original
#' # columns matching the regular exprssions will be gone
#' 
#' \dontrun{
#' # Simpler interface by gene list
#' status_df = collapse_columns_by_regex(
#'   status_df,
#'   genes = c("MYC","BCL2","TP53"),
#'   suffix = "_any"
#' )
#' }
collapse_columns_by_regex = function(x,
                                     collapse_policy = NULL,
                                     genes = NULL,
                                     suffix = "_any",
                                     drop_sources = TRUE,
                                     agg_fn = max){
  if(is.null(collapse_policy) || length(collapse_policy) == 0){
    if(!is.null(genes) && length(genes) > 0){
      collapse_policy = setNames(lapply(genes, function(g){
        paste0("^", g, "_")
      }), paste0(genes, suffix))
    }else{
      return(x)
    }
  }
  if(is.null(names(collapse_policy)) || any(names(collapse_policy) == "")){
    stop("collapse_policy must be a named list with new column names")
  }

  x_df = as.data.frame(x)
  original_cols = colnames(x_df)
  drop_cols = character(0)

  for(new_col in names(collapse_policy)){
    patterns = collapse_policy[[new_col]]
    if(length(patterns) == 0){
      next
    }
    matched = unique(unlist(lapply(patterns, function(p){
      grep(p, original_cols, value = TRUE)
    })))
    if(length(matched) == 0){
      next
    }
    drop_cols = c(drop_cols, matched)
    vals = x_df[, matched, drop = FALSE]
    x_df[[new_col]] = apply(vals, 1, agg_fn, na.rm = TRUE)
  }

  if(drop_sources && length(drop_cols) > 0){
    drop_cols = setdiff(unique(drop_cols), names(collapse_policy))
    x_df = x_df[, setdiff(colnames(x_df), drop_cols), drop = FALSE]
  }

  return(x_df)
}

#' @title Order a set of genetic features that contain gene names based on a 
#' user-provided gene ordering, with optional suffix-based tie-breaking and 
#' metadata annotation.
#'
#' @description Order feature names of the form \code{GENE_something} to match a
#' desired gene ordering, with optional suffix priority for tie-breaking.
#'
#' @param gene_order Character vector of genes defining the desired order. If a
#' feature name (e.g. \code{MYD88_L265P}) is included, it will be placed at the
#' exact index position defined by \code{gene_order}.
#' @param features Character vector of feature names (e.g. \code{MYD88_L265P},
#' \code{EZH2_hotspot}, \code{TP53_any}, \code{CREBBP_driver}).
#' @param suffix_priority Character vector defining preferred suffix order.
#' Default is \code{c("hotspot","driver","any")}. Any non-matching suffixes are
#' ordered after these, preserving input order.
#' @param gene_metadata Optional data frame with gene- or feature-level metadata
#' to join onto the output. If a column named
#' \code{feature} or \code{name} is present, exact feature matches are joined
#' first, then remaining rows are filled using \code{gene_key}.
#' @param gene_key Column name in \code{gene_metadata} used to join on gene names
#' (fallback after feature-level matching).
#'
#' @return A data frame with one row per feature in the desired order, including
#' columns \code{feature}, \code{gene}, \code{suffix}, and \code{rank}, plus any
#' joined metadata fields.
#'
#' @export
#'
order_features_by_gene = function(gene_order,
                                  features,
                                  suffix_priority = c("hotspot","driver","any"),
                                  gene_metadata = NULL,
                                  gene_key = "gene"){
  if(is.null(features) || length(features) == 0){
    return(data.frame(feature = character(0), gene = character(0), suffix = character(0), rank = integer(0)))
  }
  if(is.null(gene_order) || length(gene_order) == 0){
    gene = sub("_.*$", "", features)
    suffix = sub("^[^_]*_?", "", features)
    suffix = ifelse(suffix == features, "", suffix)
    return(data.frame(
      feature = features,
      gene = gene,
      suffix = suffix,
      rank = seq_along(features),
      stringsAsFactors = FALSE
    ))
  }

  gene = sub("_.*$", "", features)
  suffix = sub("^[^_]*_?", "", features)
  suffix = ifelse(suffix == features, "", suffix)

  gene_order_genes = gene_order[!grepl("_", gene_order)]
  gene_positions = match(gene_order_genes, gene_order)
  gene_rank = match(gene, gene_order_genes)
  rank = gene_positions[gene_rank]
  rank[is.na(rank)] = length(gene_order) + 1

  explicit_rank = match(features, gene_order)
  rank[!is.na(explicit_rank)] = explicit_rank[!is.na(explicit_rank)]

  suffix_rank = match(suffix, suffix_priority)
  suffix_rank[is.na(suffix_rank)] = length(suffix_priority) + 1

  ord = order(rank, suffix_rank, seq_along(features))
  ordered = features[ord]

  out = data.frame(
    feature = ordered,
    gene = gene[ord],
    suffix = suffix[ord],
    rank = rank[ord],
    stringsAsFactors = FALSE
  )

  if(!is.null(gene_metadata)){
    feature_key = NULL
    if("feature" %in% colnames(gene_metadata)){
      feature_key = "feature"
    }else if("name" %in% colnames(gene_metadata)){
      feature_key = "name"
    }

    metadata_cols = setdiff(colnames(gene_metadata), c(feature_key, gene_key))

    if(!is.null(feature_key) && any(gene_metadata[[feature_key]] %in% out$feature)){
      feature_meta = dplyr::distinct(gene_metadata, .data[[feature_key]], .keep_all = TRUE)
      out = dplyr::left_join(out, feature_meta, by = setNames(feature_key, "feature"))
    }

    if(gene_key %in% colnames(gene_metadata) && length(metadata_cols) > 0){
      gene_meta = dplyr::distinct(gene_metadata, .data[[gene_key]], .keep_all = TRUE)
      out = dplyr::left_join(out, gene_meta, by = setNames(gene_key, "gene"), suffix = c("", ".gene"))
      for(col in metadata_cols){
        gene_col = paste0(col, ".gene")
        if(gene_col %in% colnames(out)){
          out[[col]] = dplyr::coalesce(out[[col]], out[[gene_col]])
          out[[gene_col]] = NULL
        }
      }
    }else if(is.null(feature_key)){
      stop(paste0("gene_metadata must contain column: ", gene_key, " (or a feature/name column)"))
    }
  }

  out
}

#' @title Annotate curated drivers
#'
#' @description Annotate MAF-like data frame with curated driver/hotspot annotations.
#'
#' @details This function annotates a MAF-like data frame and will create or overwrite the
#' "hot_spot" column based on curated driver/hotspot rules. Genes for review are supplied
#' with the `genes_of_interest` parameter. Overlap is determined using a fuzzy interval
#' match so any overlap between a mutation and a hotspot region is considered a match.
#' Currently only a few sets of genes are supported, see parameter description for more information and limitations.
#' The desired genome build can be specified with `genome_build` parameter (useful if you loaded the MAF from disk).
#' Obviously, you need to specify the same genome_build as the coordinate system used in the incoming MAF.
#'
#' @param maf_data A maf_data object or data frame in MAF format.
#' @param genes_of_interest An optional vector of genes for hotspot annotation. By default, the genes present in the curated hotspot tables are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is grch37 genome build.
#' @param custom_coordinates A data frame defining custom driver/hotspot regions.
#' Follow the same structure as \code{get_driver_coordinates()}. At minimum, it must
#' include \code{Hugo_Symbol}, \code{chrom}, \code{start}, and \code{end}. Optional
#' columns include \code{alias} (identifier), \code{type} (token string), and/or
#' \code{classes} (regex for Variant_Classification). If \code{classes} is missing,
#' it will be set to all non-synonymous classes. See \code{get_driver_coordinates()}
#' for the expected column structure.
#' @param existing_values_action Character. How to handle existing columns.
#' Use "clobber" (default) to always reset "hot_spot" and "mutation_alias"
#' before re-annotating. Use "update" to only fill missing values.
#' @return The same data frame (as given to the `annotated_maf` parameter) with the reviewed columns
#' "hot_spot", "mutation_alias", and "driver_alias".
#'
#' @import dplyr
#' @export
#'
#' @examples
#'
#' # Get metadata for all DLBCLs available
#' dlbcl_meta = GAMBLR.open::get_gambl_metadata() %>% 
#'    dplyr::filter(pathology == "DLBCL")
#' # get mutations (as maf_data) for all these samples including non-coding (silent) for DLBCLone
#' dlbcl_coding = GAMBLR.open::get_all_coding_ssm(
#'                                    dlbcl_meta,
#'                                    include_silent=TRUE)
#'
#' # annotate known driver mutations using bundled annotations
#' # WARNING: currently hg38 annotations are out of sync with grch37
#' # Stick with grch37 for now
#' dlbcl_coding = annotate_curated_drivers(maf_data=dlbcl_coding)
#'
#' #look at the new columns we have
#' dplyr::filter(dlbcl_coding,
#'              !grepl("other",mutation_alias)) %>%
#'  dplyr::select(Chromosome,Start_Position,Variant_Classification,HGVSp_Short,hot_spot,mutation_alias,driver_alias)
#'
#' # Annotate TP53 driver categories across all available samples
#' all_meta = suppressMessages(GAMBLR.open::get_gambl_metadata())
#' all_maf = GAMBLR.open::get_all_coding_ssm(all_meta)
#'
#' all_maf_tp53 = annotate_curated_drivers(
#'   maf_data = all_maf,
#'   genes_of_interest = "TP53"
#' )
#'
#' # TP53 mutations are classified into truncations, named hotspot codons
#' # (R175, R248, R273, R282, G245) and structural loop regions (LSH1/2, LoopL2/3)
#' dplyr::filter(all_maf_tp53, Hugo_Symbol == "TP53") %>%
#'   dplyr::select(Tumor_Sample_Barcode, Variant_Classification,
#'                 HGVSp_Short, hot_spot, mutation_alias, driver_alias) %>%
#'   head(20)
#'
#' # Count how many mutations fall into each TP53 driver category
#' dplyr::filter(all_maf_tp53, Hugo_Symbol == "TP53") %>%
#'   dplyr::count(mutation_alias, sort = TRUE)
#'
#'
annotate_curated_drivers = function(maf_data,
                           genes_of_interest = c("MYD88", "NOTCH1", "NOTCH2","NFKBIZ"),
                           genome_build,
                           custom_coordinates,
                           existing_values_action = "clobber",
                           class_priority = c("hotspot", "truncation", "region")){
  original_has_maf_class = "maf_data" %in% class(maf_data)
  annotated_maf = strip_genomic_classes(maf_data)
  if(missing(genome_build)){
    if("maf_data" %in% class(maf_data)){
      genome_build = get_genome_build(maf_data)
      # drop our S3 classes because these additional attributes seem to
      # cause some problems when the data is subsequently munged.
      annotated_maf = strip_genomic_classes(maf_data)
    }else{
      if("genome_build" %in% colnames(maf_data)){
        gb = unique(maf_data$genome_build[!is.na(maf_data$genome_build)])
        if(length(gb) == 1){
          genome_build = gb
        }else{
          stop("genome_build is required (multiple or missing genome_build values in data)")
        }
      }else if("NCBI_Build" %in% colnames(maf_data)){
        gb = unique(maf_data$NCBI_Build[!is.na(maf_data$NCBI_Build)])
        if(length(gb) == 1){
          genome_build = gb
        }else{
          stop("genome_build is required (multiple or missing NCBI_Build values in data)")
        }
      }else{
        stop("genome_build is required")
      }
      annotated_maf = maf_data
    }
  }

  if(!missing(custom_coordinates)){
    #check for the required columns
    if(any(!"Hugo_Symbol" %in% colnames(custom_coordinates))){
      stop("coordinates data frame must have Hugo_Symbol column")
    }
    if(!"chrom" %in% colnames(custom_coordinates)){
      stop("custom_coordinates requires a column named chrom specifying the chromosome of each hot spot region")
    }
    if(!"start" %in% colnames(custom_coordinates) || !"end" %in% colnames(custom_coordinates) ){
      stop("custom_coordinates requires a start and end column specifying the boundaries of each hot spot region")
    }
    coordinates = custom_coordinates
    if(!"classes" %in% colnames(coordinates)){
      coordinates = coordinates %>%
        dplyr::mutate(classes = paste0(vc_nonSynonymous, collapse="|"))
    }
  }else{
    coordinates = get_driver_coordinates(genome_build)
  }

  if(!existing_values_action %in% c("clobber", "update")){
    stop("existing_values_action must be one of: clobber, update")
  }
  if(is.null(class_priority) || length(class_priority) == 0){
    stop("class_priority must be a non-empty character vector")
  }

  if(existing_values_action == "clobber"){
    annotated_maf = annotated_maf %>%
      dplyr::mutate(
        hot_spot = "FALSE",
        mutation_alias = paste0(Hugo_Symbol, "_other")
      )
  }else{
    if(!"hot_spot" %in% colnames(annotated_maf)){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = "FALSE")
    }
    if(!"mutation_alias" %in% colnames(annotated_maf)){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(mutation_alias = paste0(Hugo_Symbol, "_other"))
    }

    annotated_maf = annotated_maf %>%
      dplyr::mutate(
        hot_spot = ifelse(hot_spot %in% c(TRUE, "TRUE"), "TRUE", "FALSE"),
        mutation_alias = ifelse(
          hot_spot %in% c(TRUE, "TRUE"),
          mutation_alias,
          paste0(Hugo_Symbol, "_other")
        )
      )
  }

  if(!".row_id" %in% colnames(annotated_maf)){
    annotated_maf = annotated_maf %>%
      dplyr::mutate(.row_id = dplyr::row_number())
  }

  protected_cols = grep("\\.x$", colnames(annotated_maf), value = TRUE)
  protected_col_map = character(0)
  if(length(protected_cols) > 0){
    protected_names = paste0("..gamblr_protected_col_", seq_along(protected_cols), "..")
    while(any(protected_names %in% c(colnames(annotated_maf), colnames(coordinates)))){
      protected_names = paste0(".", protected_names)
    }
    protected_col_map = stats::setNames(protected_cols, protected_names)
    colnames(annotated_maf)[match(protected_cols, colnames(annotated_maf))] = protected_names
  }

  if(!missing(genes_of_interest) && !is.null(genes_of_interest)){
    # Expanded to include known aliases (e.g. old/new HGNC histone names) on
    # both sides of the comparison, since `coordinates` (bundled curated
    # data) and a caller-supplied `genes_of_interest` aren't guaranteed to
    # use the same nomenclature for a given gene.
    supported_genes = sort(unique(expand_gene_aliases(coordinates$Hugo_Symbol)))
    genes_of_interest_expanded = expand_gene_aliases(genes_of_interest)
    if (length(intersect(supported_genes, genes_of_interest_expanded))==0){
        stop(paste0("Currently only ",  paste(supported_genes, collapse=", "),
                    " are supported. Please specify one of these genes."))
    }
    if (length(setdiff(genes_of_interest_expanded, supported_genes))>0){
        message(strwrap(paste0("Currently only ", paste(supported_genes, collapse=", "),
                               " are supported. By default only these genes from the",
                               " supplied list will be reviewed. Reviewing hotspots for genes ",
                               paste(intersect(supported_genes, genes_of_interest_expanded),
                                     collapse = ", "), ", it will take a second ...")))
    }
    coordinates = coordinates %>%
      dplyr::filter(Hugo_Symbol %in% genes_of_interest_expanded)
  }

  coordinates = coordinates %>%
    dplyr::mutate(
      max = ifelse(start > end, start, end),
      min = ifelse(end < start, end, start)
    ) %>%
    dplyr::mutate(start = min, end = max) %>%
    dplyr::select(-min, -max)

  if("Chromosome" %in% colnames(annotated_maf)){
    annotated_maf = annotated_maf %>%
      dplyr::mutate(Chromosome = as.character(Chromosome))
  }
  if("chrom" %in% colnames(coordinates)){
    coordinates = coordinates %>%
      dplyr::mutate(chrom = as.character(chrom))
  }

  # Expanded to include known aliases so a maf row annotated under a gene's
  # other name (e.g. old/new HGNC histone names) isn't misclassified as
  # having no driver/hotspot coordinates.
  coordinates_genes = expand_gene_aliases(coordinates$Hugo_Symbol)
  maf_keep = annotated_maf %>%
    dplyr::filter(Hugo_Symbol %in% coordinates_genes)
  maf_skip = annotated_maf %>%
    dplyr::filter(!Hugo_Symbol %in% coordinates_genes)

  if (nrow(maf_keep) == 0){
    warning("No mutations found in genes of interest; i.e no hotspot or driver genes. Skipping hotspot and driver annotations.")

    # Just return annotated_maf with default columns already set
    annotated_maf = annotated_maf %>%
      dplyr::mutate(
        hot_spot = "FALSE",
        driver = "FALSE",
        mutation_alias = paste0(Hugo_Symbol, "_other"),
        mutation_alias_all = NA_character_,
        hotspot_alias = NA_character_,
        driver_alias = paste0(Hugo_Symbol, "_other")
      )

    if(original_has_maf_class && exists("create_maf_data", mode = "function")){
      annotated_maf = create_maf_data(annotated_maf, genome_build)
    }

    return(annotated_maf)
  }

  reviewed_maf = suppressMessages(cool_overlaps(
    maf_keep,
    coordinates,
    columns1 = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position"),
    columns2 = c("Hugo_Symbol", "chrom", "start", "end"),
    type = "fuzzy",
    nomatch = TRUE
  ))

  original_cols = setdiff(colnames(annotated_maf), ".row_id")

  reviewed_maf = reviewed_maf %>%
    dplyr::mutate(
      .is_hot = !is.na(start) & !is.na(classes) &
        mapply(grepl, pattern = classes, x = Variant_Classification),
      .class_rank = match(class, class_priority)
    ) %>%
    dplyr::group_by(.row_id) %>%
    dplyr::summarise(
      dplyr::across(all_of(original_cols), ~ dplyr::first(.x)),
      hot_spot = ifelse(any(.is_hot), "TRUE", dplyr::first(hot_spot)),
      mutation_alias = {
        matched_idx = which(.is_hot & !is.na(alias))
        if(length(matched_idx) > 0){
          order_idx = order(.class_rank[matched_idx], alias[matched_idx], na.last = TRUE)
          alias[matched_idx][order_idx][1]
        }else{
          dplyr::first(mutation_alias)
        }
      },
      mutation_alias_all = {
        matched_idx = which(.is_hot & !is.na(alias))
        if(length(matched_idx) > 0){
          order_idx = order(.class_rank[matched_idx], alias[matched_idx], na.last = TRUE)
          paste(unique(alias[matched_idx][order_idx]), collapse = ";")
        }else{
          NA_character_
        }
      },
      hotspot_alias = {
        if(any(.is_hot & class == "hotspot", na.rm = TRUE)){
          paste0(dplyr::first(Hugo_Symbol), "_hotspot")
        }else{
          NA_character_
        }
      },
      .groups = "drop"
    )

  reviewed_maf = reviewed_maf %>%
    dplyr::mutate(
      hot_spot = ifelse(!is.na(hotspot_alias), "TRUE", "FALSE"),
      driver = ifelse(!is.na(mutation_alias_all), "TRUE", "FALSE")
    )

  reviewed_maf = reviewed_maf %>%
    dplyr::arrange(.row_id) %>%
    dplyr::select(-.row_id)

  if(nrow(maf_skip) > 0){
    reviewed_maf = dplyr::bind_rows(reviewed_maf, maf_skip) %>%
      dplyr::arrange(.row_id) %>%
      dplyr::select(-.row_id)
  }

  reviewed_maf = reviewed_maf %>%
    dplyr::mutate(
      driver_alias = ifelse(
        mutation_alias == paste0(Hugo_Symbol, "_other"),
        paste0(Hugo_Symbol, "_other"),
        paste0(Hugo_Symbol, "_driver")
      )
    )

  if(length(protected_col_map) > 0){
    colnames(reviewed_maf) = dplyr::recode(colnames(reviewed_maf), !!!protected_col_map)
  }

  if(original_has_maf_class && exists("create_maf_data", mode = "function")){
    reviewed_maf = create_maf_data(reviewed_maf, genome_build)
  }

  return(reviewed_maf)
}



vc_truncating = c("Nonsense_Mutation","Splice_Site","Frame_Shift_Ins","Frame_Shift_Del")
vc_non_truncating = vc_nonSynonymous[!vc_nonSynonymous %in% vc_truncating]

#' @title Driver coordinate definitions
#'
#' @description Return curated driver coordinate regions for a given genome build,
#' including per-region variant class patterns and aliases.
#' The \code{type} field may include tokens such as missense, inframe, trunc,
#' startloss, stoploss, splicesite, spliceregion, and utr_3 (comma-separated).
#' This table can be used as a template for \code{custom_coordinates} in
#' \code{annotate_curated_drivers()}.
#'
#' @param genome_build Reference genome build for the coordinates in the MAF file.
#' Supported values include hg19/grch37/hs37d5/GRCh37 and hg38/grch38/GRCh38.
#' @return A data frame with columns including "Hugo_Symbol", "chrom", "start",
#' "end", and "classes", plus any columns present in the underlying curated
#' coordinate table (e.g. "strand", "type", "alias"). The \code{size} and
#' \code{classes} columns are added by this function.
#'
#' @import dplyr tidyr
#' @export
#'
get_driver_coordinates = function(genome_build){
  if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates = GAMBLR.utils::hotspot_regions_grch37
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates = GAMBLR.utils::hotspot_regions_hg38
  }else{
    stop("The genome build specified is not currently supported. Please provide MAF file in one of the following coordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38")
  }
  coordinates = coordinates %>%
    dplyr::rename(Hugo_Symbol = gene)

  #for some reason sometimes this data frame has the higher number as the start.
  coordinates = dplyr::mutate(coordinates,
    max = ifelse(start > end, start, end),
    min = ifelse(end < start, end, start)
  ) %>%
    dplyr::mutate(start = min, end = max) %>%
    dplyr::select(-min, -max)

  coordinates = dplyr::mutate(coordinates, size = end - start + 1)

  type_to_classes = function(type_value){
    if(is.na(type_value) || type_value == "" || type_value == "all"){
      return(paste0(vc_nonSynonymous, collapse = "|"))
    }
    tokens = unlist(strsplit(type_value, ","))
    tokens = trimws(tokens)
    classes = character(0)
    if("trunc" %in% tokens){
      classes = c(classes, vc_truncating)
    }
    if("utr_3" %in% tokens){
      classes = c(classes, "3'UTR")
    }
    if("missense" %in% tokens){
      classes = c(classes, "Missense_Mutation")
    }
    if("inframe" %in% tokens){
      classes = c(classes, "In_Frame_Del", "In_Frame_Ins")
    }
    if("stoploss" %in% tokens){
      classes = c(classes, "Nonstop_Mutation")
    }
    if("startloss" %in% tokens){
      classes = c(classes, "Translation_Start_Site")
    }
    if("splicesite" %in% tokens){
      classes = c(classes, "Splice_Site")
    }
    if("spliceregion" %in% tokens){
      classes = c(classes, "Splice_Region")
    }
    if(length(classes) == 0){
      classes = vc_nonSynonymous
    }
    paste0(unique(classes), collapse = "|")
  }

  coordinates = coordinates %>%
    dplyr::mutate(classes = vapply(type, type_to_classes, character(1)))

  return(coordinates)

}

