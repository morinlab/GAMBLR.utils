#' @title Summarize mutation status
#'
#' @description Create a wide sample-by-variant matrix summarizing coding,
#' driver, and hotspot-level mutation status from an annotated MAF.
#'
#' @details This function expects output from \code{annotate_curated_hotspots()}
#' (or equivalent) with \code{mutation_alias} and \code{mutation_annotation}
#' columns. It produces a wide matrix with 0/1 indicators for each
#' sample-variant combination, including per-gene coding status, driver status,
#' and granular hotspot aliases. If you want to collapse multiple columns for a
#' gene into a single summary column (e.g., \code{MYC_any}), see
#' \code{collapse_columns_by_regex()}.
#'
#' @param annotated_maf A data frame in MAF format that includes
#' \code{Hugo_Symbol}, \code{Tumor_Sample_Barcode}, \code{Variant_Classification},
#' \code{mutation_alias}, and \code{mutation_annotation}.
#' @param genes_coding Optional character vector of genes to include in the
#' coding summary. Defaults to Tier 1 lymphoma genes.
#' @param genes_noncoding Optional character vector of genes to include for
#' non-coding mutation summaries. If NULL, non-coding mutations are not included.
#' @param gene_reporting_policy Optional named list defining per-gene reporting
#' mode. Valid keys are \code{detailed} (granular alias columns) and \code{coding}
#' (keep only \code{GENE_coding}). Genes not listed default to driver/unknown.
#' @param genes_drop_unannotated Optional character vector of genes for which
#' catch-all columns (e.g. \code{GENE_other} or \code{GENE_unknown}) should
#' be dropped from the output. Mutually exclusive with \code{gene_drop_policy}.
#' @param gene_drop_policy Optional named list defining gene drop policies, e.g.
#' \code{list(all = c("IGLL5"), unannotated = c("KMT2D"))}. Valid keys are
#' \code{all} (drop all columns for those genes) and \code{unannotated} (drop only
#' catch-all columns). Mutually exclusive with \code{genes_drop_unannotated}.
#' @param encoding_policy Named list controlling the numeric values used for
#' coding, non-coding, and driver indicators. Default is
#' \code{list(coding = 2, noncoding = 1, driver = 1)}.
#' @param verbose Logical. If TRUE, prints messages describing which columns
#' will be dropped based on reporting and drop policies.
#'
#' @return A wide data frame (matrix-like) with rows as samples and columns as
#' variant categories, containing 0/1 indicators.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' maf = GAMBLR.open::get_coding_ssm(
#'   these_samples_metadata = GAMBLR.open::get_gambl_metadata(),
#'   this_seq_type = "capture"
#' )
#' maf_anno = annotate_curated_drivers(maf)
#' summary_mat = summarise_mutation_status(maf_anno)
#' summary_mat[1:5, 1:5]
#' }
#'
#' \dontrun{
#' maf_anno = annotate_curated_drivers(GAMBLR.open::get_coding_ssm())
#' summary_mat = summarise_mutation_status(
#'   maf_anno,
#'   genes_coding = c("MYD88","NOTCH1","FOXO1"),
#'   gene_reporting_policy = list(detailed = c("FOXO1"), coding = c("EZH2")),
#'   genes_drop_unannotated = c("NOTCH1")
#' )
#' }
#'
#' \dontrun{
#' # Include selected non-coding genes and custom indicator values
#' maf_anno = annotate_curated_drivers(GAMBLR.open::get_coding_ssm(include_silent=TRUE))
#' summary_mat = summarise_mutation_status(
#'   maf_anno,
#'   genes_noncoding = c("PIM1","SOCS1"),
#'   encoding_policy = list(coding = 2, noncoding = 1, driver = 1)
#' )
#' }
#'
summarise_mutation_status = function(annotated_maf,
                                     genes_coding = NULL,
                                     genes_noncoding = NULL,
                                     gene_reporting_policy = NULL,
                                     genes_drop_unannotated = NULL,
                                     gene_drop_policy = NULL,
                                     encoding_policy = list(coding = 2, noncoding = 1, driver = 1),
                                     verbose = FALSE){
  if(is.null(genes_coding)){
    genes_coding = lymphoma_genes %>%
    dplyr::filter(FL_Tier==1 | MCL_Tier== 1 | BL_Tier==1 | DLBCL_Tier==1) %>%
    pull(Gene)
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
  overlap = intersect(detailed_genes, coding_policy_genes)
  if(length(overlap) > 0){
    stop(paste0("Genes cannot be in both gene_reporting_policy$detailed and $coding: ", paste(overlap, collapse=", ")))
  }

  if(verbose){
    if(length(detailed_genes) > 0){
      message("Will keep detailed (alias-level) columns for genes: ", paste(detailed_genes, collapse = ", "))
      message("Will drop _driver/_unknown columns for genes: ", paste(detailed_genes, collapse = ", "))
    }
    if(length(coding_policy_genes) > 0){
      message("Will keep _coding columns (and drop driver/unknown and detailed) for genes: ", paste(coding_policy_genes, collapse = ", "))
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


  #expects the output of Annotate curated drivers
  #drop Silent and hang onto them for later
  if(!is.null(genes_noncoding)){

    noncoding_maf = dplyr::filter(annotated_maf, ! Variant_Classification %in% vc_nonSynonymous,Hugo_Symbol %in% genes_noncoding)

    noncoding_count = noncoding_maf %>%
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      unique() %>%
      mutate(Variant=paste0(Hugo_Symbol,"_noncoding")) %>%
      ungroup() %>%
      dplyr::select(-Hugo_Symbol) %>%
      mutate(mutated = encoding_policy$noncoding) # will fill missing with 0 at the join stage

  }


  coding_maf = dplyr::filter(annotated_maf, Variant_Classification %in% vc_nonSynonymous)

  #separately count:
  # - coding variants
  # - drivers
  # - non-drivers
  # - each distinct hotspot

  coding_maf = coding_maf %>%
    dplyr::filter(Hugo_Symbol %in% genes_coding) %>%
    strip_genomic_classes()

  

  detailed_driver_count = coding_maf %>%
    dplyr::filter(Hugo_Symbol %in% detailed_genes) %>%
    dplyr::select(mutation_alias, Tumor_Sample_Barcode) %>%
    group_by(mutation_alias, Tumor_Sample_Barcode) %>%
    unique() %>%
    rename(Variant=mutation_alias) %>%
    mutate(mutated = encoding_policy$driver) # will fill missing with 0 at the join stage
  driver_report_genes = dplyr::filter(coding_maf, grepl("_driver", mutation_annotation)) %>%
    pull(Hugo_Symbol) %>%
    unique()
  if(length(detailed_genes) > 0){
    driver_report_genes = setdiff(driver_report_genes, detailed_genes)
  }
  if(length(coding_policy_genes) > 0){
    driver_report_genes = setdiff(driver_report_genes, coding_policy_genes)
  }
  driver_count = coding_maf %>%
    dplyr::filter(Hugo_Symbol %in% driver_report_genes) %>%
    dplyr::select(mutation_annotation, Tumor_Sample_Barcode) %>%
    group_by(mutation_annotation, Tumor_Sample_Barcode) %>%
    unique() %>%
    rename(Variant=mutation_annotation) %>%
    mutate(mutated = encoding_policy$driver) # will fill missing with 0 at the join stage
  
  skip_genes_coding = union(driver_report_genes, detailed_genes)
  coding_count = coding_maf %>%
    dplyr::filter(!Hugo_Symbol %in% skip_genes_coding) %>%
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    unique() %>%
    mutate(Variant=paste0(Hugo_Symbol,"_coding")) %>%
    ungroup() %>%
    dplyr::select(-Hugo_Symbol) %>%
    mutate(mutated = encoding_policy$coding) # will fill missing with 0 at the join stage
  
  if(!is.null(genes_noncoding)){
    long = bind_rows(
      noncoding_count,
      coding_count,
      detailed_driver_count,
      driver_count
    )
  }else{
    long = bind_rows(
      coding_count,
      detailed_driver_count,
      driver_count
    )
  }


  all_samples = distinct(coding_maf, Tumor_Sample_Barcode)
  all_variants = distinct(long, Variant)

  wide = long %>%
    tidyr::complete(all_samples, all_variants, fill = list(mutated = 0)) %>%
    tidyr::pivot_wider(names_from = Variant, values_from = mutated, values_fill = 0) %>%
    column_to_rownames("Tumor_Sample_Barcode")

  if(length(gene_drop_policy) > 0){
    drop_cols = character(0)
    if(!is.null(gene_drop_policy$unannotated)){
      drop_cols = c(drop_cols, unlist(lapply(gene_drop_policy$unannotated, function(g){
        c(paste0(g, "_other"), paste0(g, "_unknown"))
      })))
      if(verbose && length(gene_drop_policy$unannotated) > 0){
        message("Will drop _other/_unknown columns for genes: ", paste(gene_drop_policy$unannotated, collapse = ", "))
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
#' a positive status (default \code{"POS"}).
#' @param encoding_policy Named list controlling the numeric values used for
#' positive and negative samples. Default is
#' \code{list(positive = 2, negative = 0)}.
#' @param sample_id_col Column name in \code{metadata} containing sample IDs
#' (default \code{"sample_id"}).
#' @param debug Logical. If TRUE, prints matching diagnostics.
#'
#' @return The input \code{mat} with new columns appended.
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
                                            debug = FALSE){
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

  if(debug){
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
    if(debug){
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
#' (default \code{"any"}, producing columns like \code{GENE_any}).
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
#'   suffix = "any"
#' )
#' }
collapse_columns_by_regex = function(x,
                                     collapse_policy = NULL,
                                     genes = NULL,
                                     suffix = "any",
                                     drop_sources = TRUE,
                                     agg_fn = max){
  if(is.null(collapse_policy) || length(collapse_policy) == 0){
    if(!is.null(genes) && length(genes) > 0){
      collapse_policy = setNames(lapply(genes, function(g){
        paste0("^", g, "_")
      }), paste0(genes, "_", suffix))
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

#' @title Annotate curated drivers
#'
#' @description Annotate MAF-like data frame with a hot_spot column indicating recurrent mutations.
#'
#' @details This function annotates a MAF-like data frame and will create or overwrite the
#' "hot_spot" column based on curated hotspot rules. Genes for hotspot review are supplied
#' with the `genes_of_interest` parameter. Overlap is determined using a fuzzy interval
#' match so any overlap between a mutation and a hotspot region is considered a match.
#' Currently only a few sets of genes are supported, see parameter description for more information and limitations.
#' The desired genome build can be specified with `genome_build` parameter (useful if you loaded the MAF from disk).
#' Obviously, you need to specify the same genome_build as the coordinate system used in the incoming MAF.
#'
#' @param annotated_maf A maf_data object or data frame in MAF format that has hotspots annotated using the function annotate_hotspots().
#' @param genes_of_interest An optional vector of genes for hotspot annotation. By default, the genes present in the curated hotspot tables are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is grch37 genome build.
#' @param custom_coordinates A data frame with custom coordinates for the hot spots.
#' All mutations in any of the regions specified in the data frame will be marked as hot spots.
#' The data frame must have the following columns: "Hugo_Symbol", "chrom", "start", and "end".
#' Optional columns are: "classes" (regex for Variant_Classification) and "alias" (identifier for the hotspot).
#' The "type" column can include tokens such as missense, inframe, trunc, startloss,
#' stoploss, splicesite, and spliceregion (comma-separated).
#' @param existing_values_action Character. How to handle existing columns.
#' Use "clobber" (default) to always reset "hot_spot" and "mutation_alias"
#' before re-annotating. Use "update" to only fill missing values.
#' @return The same data frame (as given to the `annotated_maf` parameter) with the reviewed columns
#' "hot_spot", "mutation_alias", and "mutation_annotation".
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
#'  dplyr::select(Chromosome,Start_Position,Variant_Classification,HGVSp_Short,hot_spot,mutation_alias,mutation_annotation)
#'
#'
annotate_curated_drivers = function(maf_data,
                           genes_of_interest = c("MYD88", "NOTCH1", "NOTCH2"),
                           genome_build,
                           custom_coordinates,
                           existing_values_action = "clobber"){
  original_has_maf_class = "maf_data" %in% class(maf_data)
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
    coordinates = get_hotspot_coordinates(genome_build)
  }

  if(!existing_values_action %in% c("clobber", "update")){
    stop("existing_values_action must be one of: clobber, update")
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

  if(!missing(genes_of_interest) && !is.null(genes_of_interest)){
    supported_genes = sort(unique(coordinates$Hugo_Symbol))
    if (length(intersect(supported_genes, genes_of_interest))==0){
        stop(paste0("Currently only ",  paste(supported_genes, collapse=", "),
                    " are supported. Please specify one of these genes."))
    }
    if (length(setdiff(genes_of_interest, supported_genes))>0){
        message(strwrap(paste0("Currently only ", paste(supported_genes, collapse=", "),
                               " are supported. By default only these genes from the",
                               " supplied list will be reviewed. Reviewing hotspots for genes ",
                               paste(intersect(supported_genes, genes_of_interest),
                                     collapse = ", "), ", it will take a second ...")))
    }
    coordinates = coordinates %>%
      dplyr::filter(Hugo_Symbol %in% genes_of_interest)
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

  maf_keep = annotated_maf %>%
    dplyr::filter(Hugo_Symbol %in% coordinates$Hugo_Symbol)
  maf_skip = annotated_maf %>%
    dplyr::filter(!Hugo_Symbol %in% coordinates$Hugo_Symbol)

  reviewed_maf = cool_overlaps(
    maf_keep,
    coordinates,
    columns1 = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position"),
    columns2 = c("Hugo_Symbol", "chrom", "start", "end"),
    type = "fuzzy",
    nomatch = TRUE
  )

  original_cols = setdiff(colnames(annotated_maf), ".row_id")

  reviewed_maf = reviewed_maf %>%
    dplyr::mutate(
      .is_hot = !is.na(start) & !is.na(classes) &
        mapply(grepl, pattern = classes, x = Variant_Classification)
    ) %>%
    dplyr::group_by(.row_id) %>%
    dplyr::summarise(
      dplyr::across(all_of(original_cols), ~ dplyr::first(.x)),
      hot_spot = ifelse(any(.is_hot), "TRUE", dplyr::first(hot_spot)),
      mutation_alias = {
        aliases = alias[.is_hot & !is.na(alias)]
        if(length(aliases) > 0) aliases[1] else dplyr::first(mutation_alias)
      },
      .groups = "drop"
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
      mutation_annotation = ifelse(
        mutation_alias == paste0(Hugo_Symbol, "_other"),
        paste0(Hugo_Symbol, "_unknown"),
        paste0(Hugo_Symbol, "_driver")
      )
    )

  if(original_has_maf_class && exists("create_maf_data", mode = "function")){
    reviewed_maf = create_maf_data(reviewed_maf, genome_build)
  }

  return(reviewed_maf)
}



vc_truncating = c("Nonsense_Mutation","Splice_Site","Frame_Shift_Ins","Frame_Shift_Del")
vc_non_truncating = vc_nonSynonymous[!vc_nonSynonymous %in% vc_truncating]

#' @title Hotspot coordinate definitions
#'
#' @description Return curated hotspot coordinate regions for a given genome build,
#' including per-region variant class patterns and aliases.
#'
#' @param genome_build Reference genome build for the coordinates in the MAF file.
#' Supported values include hg19/grch37/hs37d5/GRCh37 and hg38/grch38/GRCh38.
#' @return A data frame with columns including "Hugo_Symbol", "chrom", "start",
#' "end", "strand", "type", "size", "alias", and "classes".
#'
#' @import dplyr tidyr
#' @export
#'
get_hotspot_coordinates = function(genome_build){
  if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates = GAMBLR.utils::hotspot_regions_grch37
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates = GAMBLR.utils::hotspot_regions_hg38
  }else{
    stop("The genome build specified is not currently supported. Please provide MAF file in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38")
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
