#' Shared genetic feature definitions (grch37)
#'
#' Definitions of the mutation, structural variant and copy number features
#' used by DLBCL classifiers, one row per feature. Each feature is defined
#' once and joined to the classifiers that use it through
#' [classifier_features]. Built by `data-raw/classifier_features.R` from
#' `inst/extdata/genetic_features.grch37.tsv`. Use [get_classifier_features]
#' rather than this object directly.
#'
#' @format A data frame with one row per feature:
#' \describe{
#'   \item{feature}{Unique feature name, e.g. "CREBBP_nonsilent",
#'     "NOTCH1_trunc_PEST" or "BCL2_SV".}
#'   \item{data_type}{"ssm", "sv" or "cnv".}
#'   \item{gene}{Gene symbol (as in GAMBL grch37 MAFs, e.g. HIST1H1D).}
#'   \item{entrez}{NCBI Entrez gene ID, where known.}
#'   \item{variant_class}{For ssm features, comma-separated mutation classes
#'     that count (see [classify_mutations]).}
#'   \item{protein_change, protein_change_exclude}{For ssm features, a
#'     protein change that is required (e.g. "p.L265P") or excluded.}
#'   \item{chrom, start, end}{For ssm features, an optional grch37 range the
#'     mutation start must fall in (e.g. the LymphGen NOTCH1 PEST rule).}
#'   \item{position_rule}{The named rule the range implements, so it can be
#'     switched off (see `position_rules` in [assemble_classifier_features]):
#'     "PEST" (NOTCH1/NOTCH2), "ITAM" (CD79B), "SET" (EZH2) or "tss4kb"
#'     (silent and 5' UTR mutations in the first 4 kb of the gene from its
#'     TSS, as in [subset_ssm_by_tss] with its defaults; the LymphGen
#'     Synon definition).}
#'   \item{cnv_feature, cnv_feature_set}{For cnv features, the `name` and
#'     `feature_set` of the [cnv_features_grch37] row that defines it.}
#'   \item{requires}{The input the feature is called from: "maf_data",
#'     "noncoding_maf" (e.g. 5' UTR features), "seg_data", or
#'     "sv_data or sv_from_metadata".}
#'   \item{summary}{Notes.}
#' }
"genetic_features_grch37"

#' Classifier feature map
#'
#' Which [genetic_features_grch37] features each classifier uses and how it
#' encodes them. Built by `data-raw/classifier_features.R` from
#' `inst/extdata/classifier_features.tsv`.
#'
#' A classifier feature can combine several shared features, taking the
#' largest encoded value: e.g. the DLBClass "CREBBP" feature is 2 for a
#' non-silent mutation (CREBBP_nonsilent) and otherwise 1 for a silent one
#' (CREBBP_silent).
#'
#' The LymphGen features are those its R code (Model9) selects, with
#' complete gene lists, in each of its nine submodels: the full model (copy
#' number and both translocations available), the A53 model, and the
#' submodels LymphGen uses for samples without copy number and/or BCL2 or
#' BCL6 translocation status, which select some other features (e.g. BCL2
#' Synon, a stand-in for the BCL2 translocation). `submodels` records which
#' use each feature. LymphGen refits its model to the genes available, so
#' smaller panels can use still other features. Its feature types map to shared features as follows ("Sub"
#' types count the same as the plain ones for new data): TRUNC, truncating
#' mutations; MUTATION, any truncating or non-truncating coding mutation;
#' Synon, any of these or a silent or 5' UTR mutation near the TSS; HETLOSS
#' (which includes HOMDEL), GAIN (which includes AMP), HOMDEL and AMP, copy
#' number levels from the "LymphGen_model" rows of [cnv_features_grch37].
#'
#' @format A data frame with one row per classifier feature and shared
#' feature:
#' \describe{
#'   \item{classifier}{"DLBClass" or "LymphGen".}
#'   \item{classifier_feature}{The classifier's own feature name, e.g.
#'     "MYD88.L265P" or "X17P.DEL" (DLBClass). LymphGen features use the
#'     LymphGen code's own names, "<gene or arm>:<types>", where several types
#'     separated by "|" are alternatives and "/" joins types that are combined
#'     into one feature, e.g. "TNFAIP3:TRUNC|HOMDEL" or "BTG1:HETLOSS/SubSyn".}
#'   \item{group}{For LymphGen, the subtype whose model uses the feature
#'     (BN2, EZB, MCD, N1, ST2 or A53).}
#'   \item{feature_order}{Position of the classifier feature in the
#'     classifier's order (the DLBClass GSM order for DLBClass).}
#'   \item{feature}{The shared feature (joins to [genetic_features_grch37]).}
#'   \item{value}{The value the classifier gives the feature when present: a
#'     number, or "level" to pass the copy number level (1 or 2) through.}
#'   \item{min_level}{For copy number features, the minimum level (1 or 2)
#'     that counts; empty means 1.}
#'   \item{submodels}{For LymphGen, the submodels whose model uses the
#'     feature, comma-separated: "Full" (copy number and both translocations
#'     available), "NoCGH" (no copy number), "NoBCL2Fus", "NoBCL6Fus" and
#'     "NoFus" (BCL2, BCL6 or both translocation statuses missing), the same
#'     three ending in "CGH" (also no copy number), and "A53".}
#' }
"classifier_features"

# mutation classes used in genetic_features_grch37$variant_class
.mutation_classes = c("trunc", "missense", "inframe", "splice_region",
                      "other_nonsilent", "silent", "utr5")

#' @title Get classifier feature definitions
#'
#' @description Return the shared genetic feature definitions, joined to the
#' features of one or more classifiers.
#'
#' @param classifier Optional classifier name(s), e.g. "DLBClass" or
#' "LymphGen". Default (missing) returns every shared feature definition
#' without classifier columns.
#' @param genome_build Genome build. Default "grch37", the only one
#' currently available.
#' @param submodel Optional LymphGen submodel name(s) (see `submodels` in
#' [classifier_features]), e.g. "Full" or c("NoCGH", "A53"), to keep only
#' the features those submodels use. Default NULL keeps every feature.
#' Features of classifiers without submodels (e.g. DLBClass) are kept.
#'
#' @return A data frame: [genetic_features_grch37] rows, or with
#' `classifier`, one row per classifier feature and shared feature with the
#' [classifier_features] columns first, in each classifier's feature order.
#'
#' @export
#'
#' @examples
#' # every DLBClass feature, with its definition
#' dlbclass = get_classifier_features("DLBClass")
#' head(dlbclass[, c("classifier_feature", "feature", "value", "data_type")])
#'
#' # shared definitions used by both classifiers
#' both = intersect(get_classifier_features("DLBClass")$feature,
#'                  get_classifier_features("LymphGen")$feature)
get_classifier_features = function(classifier, genome_build = "grch37", submodel = NULL){
  if(!genome_build %in% c("grch37", "hg19")){
    stop("only grch37 classifier features are currently available")
  }
  features = GAMBLR.utils::genetic_features_grch37
  if(missing(classifier)) return(features)
  feature_map = GAMBLR.utils::classifier_features
  unknown = setdiff(classifier, feature_map$classifier)
  if(length(unknown) > 0){
    stop("unknown classifier: ", paste(unknown, collapse = ", "),
         ". Available: ", paste(unique(feature_map$classifier), collapse = ", "))
  }
  feature_map = feature_map[feature_map$classifier %in% classifier, , drop = FALSE]
  if(!is.null(submodel)){
    known = unique(trimws(unlist(strsplit(stats::na.omit(feature_map$submodels), ","))))
    unknown = setdiff(submodel, known)
    if(length(unknown) > 0){
      stop("unknown submodel: ", paste(unknown, collapse = ", "),
           ". Available: ", paste(known, collapse = ", "))
    }
    uses = vapply(strsplit(ifelse(is.na(feature_map$submodels), "", feature_map$submodels), ","),
                  function(x) any(trimws(x) %in% submodel), logical(1))
    feature_map = feature_map[is.na(feature_map$submodels) | uses, , drop = FALSE]
  }
  feature_map = feature_map[order(match(feature_map$classifier, classifier),
                                  feature_map$feature_order), , drop = FALSE]
  out = cbind(feature_map,
              features[match(feature_map$feature, features$feature),
                       setdiff(colnames(features), "feature"), drop = FALSE])
  rownames(out) = NULL
  out
}

#' @title Classify mutations for classifier features
#'
#' @description Assign each MAF row one of the mutation classes used by
#' [genetic_features_grch37]: "trunc" (nonsense, frameshift, splice site or
#' start codon loss), "missense", "inframe" (in-frame indels),
#' "splice_region" (splice region variants with a protein change),
#' "other_nonsilent" (any other mutation with a protein change, e.g. nonstop,
#' or a de novo start), "silent" or "utr5" (5' UTR). Other rows (e.g. intronic
#' or 3' UTR) get NA.
#'
#' @param maf_data MAF data frame (Variant_Classification and a protein
#' change column).
#' @param protein_change_column Column with the protein change, e.g.
#' "p.L265P". Default NULL uses Protein_Change, or HGVSp_Short (as in GAMBL
#' MAFs) if there is no Protein_Change column. Without either, every
#' splice region variant is left unclassified.
#'
#' @return `maf_data` with added columns mutation_class and protein_change.
#'
#' @export
#'
#' @examples
#' maf = data.frame(Variant_Classification = c("Nonsense_Mutation",
#'                                             "Splice_Region", "Silent"),
#'                  HGVSp_Short = c("p.Q12*", "", "p.L5L"))
#' classify_mutations(maf)
classify_mutations = function(maf_data, protein_change_column = NULL){
  maf_data = as.data.frame(maf_data)
  if(is.null(protein_change_column)){
    protein_change_column = intersect(c("Protein_Change", "HGVSp_Short"),
                                      colnames(maf_data))[1]
  }
  if(!is.na(protein_change_column) && protein_change_column %in% colnames(maf_data)){
    protein_change = as.character(maf_data[[protein_change_column]])
  }else{
    protein_change = rep(NA_character_, nrow(maf_data))
  }
  has_change = !is.na(protein_change) & protein_change != ""
  vc = as.character(maf_data$Variant_Classification)
  maf_data$mutation_class = dplyr::case_when(
    vc %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
              "Splice_Site", "Translation_Start_Site",
              "Start_Codon_Del") ~ "trunc",
    vc == "Missense_Mutation" ~ "missense",
    vc %in% c("In_Frame_Del", "In_Frame_Ins") ~ "inframe",
    vc == "Silent" ~ "silent",
    vc == "5'UTR" ~ "utr5",
    vc == "Splice_Region" & has_change ~ "splice_region",
    vc %in% c("De_novo_Start_InFrame", "De_novo_Start_OutOfFrame") ~ "other_nonsilent",
    has_change ~ "other_nonsilent",
    TRUE ~ NA_character_
  )
  maf_data$protein_change = protein_change
  maf_data
}

#' @title Assemble a sample by feature table of classifier features
#'
#' @description Call the shared genetic features ([genetic_features_grch37])
#' for every sample, from mutations, copy number segments and structural
#' variants, as the common input for classifier input files (e.g. the
#' DLBClass GSM, via [classifier_feature_matrix]) and for developing new
#' models. Every sample in `these_samples_metadata` is included, with each
#' feature marked as assessed or not.
#'
#' @details
#' \itemize{
#'   \item Mutations: each MAF row is classified with [classify_mutations]
#'     and a feature is present (1) when a mutation in its gene has one of
#'     its classes and passes its protein change and position rules.
#'     Features that need `noncoding_maf` (e.g. 5' UTR) are only called when
#'     it is given; `maf_data` should include silent mutations (e.g.
#'     `get_coding_ssm(include_silent = TRUE)`).
#'   \item Structural variants: present (1) when either `sv_data` has an SV
#'     rearranging the gene (its `gene` column, with an annotated partner,
#'     from [annotate_sv], which is run if `sv_data` has no gene column) or
#'     the metadata column in `sv_from_metadata` is positive. A gene named
#'     only as the `partner` (the locus driving the rearranged gene, e.g.
#'     BCL6 in a BCL6-MYC) doesn't count. Assessed when either source has
#'     information for the sample.
#'   \item Copy number: the feature's level (0, 1 or 2) from
#'     [segmented_data_to_cnv_status], called separately for each
#'     [cnv_features_grch37] set used (so that, e.g., DLBClass features are
#'     measured as DLBClass measures them).
#' }
#' A sample counts as assessed for a data type when the matching
#' availability column of `these_samples_metadata` (ssm_available,
#' sv_available, cn_available; see GAMBLR.results annotate_data_availability)
#' is TRUE, or, without that column, when the data was given (mutations) or
#' has rows for the sample (segments).
#'
#' @param these_samples_metadata Metadata with one row per sample_id; defines
#' the samples. Data for other samples is dropped.
#' @param classifier Optional classifier name(s), e.g. "DLBClass", to call
#' only the features those classifiers use. Default NULL calls every feature
#' whose input is given.
#' @param features Optional feature definitions to use instead (a subset of
#' [get_classifier_features], with at least the [genetic_features_grch37]
#' columns).
#' @param maf_data MAF data frame of coding mutations, including silent ones.
#' @param noncoding_maf Optional MAF of non-coding mutations (e.g. 5' UTR),
#' for the features that need it. Not used by default.
#' @param seg_data Optional segmented copy number data, as for
#' [segmented_data_to_cnv_status].
#' @param sv_data Optional SV data (bedpe, as for [annotate_sv], or already
#' annotated with gene and partner columns).
#' @param sv_from_metadata Named vector mapping SV genes to metadata columns
#' (e.g. break-apart FISH). Default `c(BCL2 = "bcl2_ba", BCL6 = "bcl6_ba",
#' MYC = "myc_ba")`. NULL to use only `sv_data`.
#' @param sv_positive_values Metadata values that indicate an SV. Default
#' "POS".
#' @param genome_build Genome build. Default "grch37" (the only one
#' currently supported).
#' @param protein_change_column Passed to [classify_mutations].
#' @param cnv_args Optional named list of further arguments for
#' [segmented_data_to_cnv_status], e.g. `value_scale`, `centering` or
#' `blacklist`.
#' @param position_rules The named position rules (`position_rule` in
#' [genetic_features_grch37]) to apply. Default all of them: "PEST", "ITAM",
#' "SET" and "tss4kb". Features whose rule is left out count mutations at any
#' position, e.g. `position_rules = c("PEST", "ITAM", "SET")` counts silent
#' and 5' UTR mutations anywhere in the gene for LymphGen Synon features, and
#' `position_rules = "PEST"` matches what the LymphGen code applies.
#'
#' @return A data frame with one row per sample and feature: sample_id,
#' feature, data_type, value (0/1, or the copy number level), assessed
#' (logical) and source (where a positive call came from: "maf_data",
#' "noncoding_maf", "seg_data", "sv_data", "metadata" or
#' "sv_data,metadata"; NA when absent).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta = GAMBLR.results::get_gambl_metadata() %>%
#'   dplyr::filter(pathology == "DLBCL") %>%
#'   GAMBLR.results::annotate_data_availability()
#' maf = GAMBLR.results::get_coding_ssm(these_samples_metadata = meta,
#'                                      include_silent = TRUE)
#' segs = GAMBLR.results::get_cn_segments(these_samples_metadata = meta)
#' sv = GAMBLR.results::get_manta_sv(these_samples_metadata = meta)
#' feature_table = assemble_classifier_features(meta, maf_data = maf,
#'                                              seg_data = segs, sv_data = sv)
#' gsm = classifier_feature_matrix(feature_table, "DLBClass",
#'                                 output_format = "gsm")
#' }
assemble_classifier_features = function(these_samples_metadata,
                                        classifier = NULL,
                                        features = NULL,
                                        maf_data = NULL,
                                        noncoding_maf = NULL,
                                        seg_data = NULL,
                                        sv_data = NULL,
                                        sv_from_metadata = c(BCL2 = "bcl2_ba",
                                                             BCL6 = "bcl6_ba",
                                                             MYC = "myc_ba"),
                                        sv_positive_values = "POS",
                                        genome_build = "grch37",
                                        protein_change_column = NULL,
                                        cnv_args = list(),
                                        position_rules = c("PEST", "ITAM", "SET", "tss4kb")){
  if(!"sample_id" %in% colnames(these_samples_metadata)){
    stop("these_samples_metadata needs a sample_id column")
  }
  these_samples_metadata = as.data.frame(these_samples_metadata)
  samples = unique(as.character(these_samples_metadata$sample_id))
  if(is.null(features)){
    features = if(is.null(classifier)) get_classifier_features(genome_build = genome_build)
               else get_classifier_features(classifier, genome_build)
  }
  features = features[!duplicated(features$feature), , drop = FALSE]
  # per-sample availability from the metadata, or NULL without the column
  available = function(column){
    if(!column %in% colnames(these_samples_metadata)) return(NULL)
    stats::setNames(these_samples_metadata[[column]] %in% TRUE,
                    as.character(these_samples_metadata$sample_id))[samples]
  }

  calls = list()
  add_calls = function(feature, present, assessed, source){
    calls[[length(calls) + 1]] <<- data.frame(
      sample_id = samples,
      feature = feature,
      value = as.numeric(present),
      assessed = assessed,
      source = ifelse(present > 0, source, NA_character_),
      stringsAsFactors = FALSE
    )
  }

  # mutations ---------------------------------------------------------------
  ssm = features[features$data_type == "ssm", , drop = FALSE]
  if(nrow(ssm) > 0){
    ssm_assessed = available("ssm_available")
    inputs = list(maf_data = maf_data, noncoding_maf = noncoding_maf)
    for(input in names(inputs)){
      these = ssm[ssm$requires == input, , drop = FALSE]
      if(nrow(these) == 0) next
      if(is.null(inputs[[input]])){
        message(nrow(these), " mutation features need ", input,
                " and are not assessed")
        for(f in these$feature) add_calls(f, rep(0, length(samples)), FALSE, input)
        next
      }
      maf = classify_mutations(inputs[[input]], protein_change_column)
      maf = maf[maf$Tumor_Sample_Barcode %in% samples &
                  maf$Hugo_Symbol %in% these$gene &
                  !is.na(maf$mutation_class), , drop = FALSE]
      maf_chrom = sub("^chr", "", as.character(maf$Chromosome))
      assessed = if(is.null(ssm_assessed)) rep(TRUE, length(samples)) else ssm_assessed
      for(i in seq_len(nrow(these))){
        row = these[i, ]
        hit = maf$Hugo_Symbol == row$gene &
          maf$mutation_class %in% trimws(strsplit(row$variant_class, ",")[[1]])
        if(!is.na(row$protein_change)){
          hit = hit & maf$protein_change %in% row$protein_change
        }
        if(!is.na(row$protein_change_exclude)){
          hit = hit & !maf$protein_change %in% row$protein_change_exclude
        }
        if(!is.na(row$start) &&
           (is.na(row$position_rule) || row$position_rule %in% position_rules)){
          hit = hit & maf_chrom == sub("^chr", "", row$chrom) &
            maf$Start_Position >= row$start & maf$Start_Position <= row$end
        }
        add_calls(row$feature, samples %in% maf$Tumor_Sample_Barcode[hit],
                  assessed, input)
      }
    }
  }

  # structural variants -----------------------------------------------------
  sv = features[features$data_type == "sv", , drop = FALSE]
  if(nrow(sv) > 0){
    sv_genes = NULL
    if(!is.null(sv_data)){
      sv_data = as.data.frame(sv_data)
      if(!"gene" %in% colnames(sv_data)){
        sv_data = annotate_sv(sv_data, genome_build = genome_build)
      }
      # only rearranged genes with an annotated partner count; the partner is
      # the locus driving the gene (e.g. BCL6 in a BCL6-MYC), not a target
      sv_data = sv_data[sv_data$tumour_sample_id %in% samples &
                          !is.na(sv_data$partner), , drop = FALSE]
      sv_genes = split(sv_data$gene, sv_data$tumour_sample_id)
      sv_assessed = available("sv_available")
      if(is.null(sv_assessed)) sv_assessed = rep(TRUE, length(samples))
    }
    meta_rows = match(samples, as.character(these_samples_metadata$sample_id))
    for(i in seq_len(nrow(sv))){
      gene = sv$gene[i]
      from_sv = rep(FALSE, length(samples))
      assessed = rep(FALSE, length(samples))
      if(!is.null(sv_genes)){
        from_sv = vapply(samples, function(s) gene %in% sv_genes[[s]], logical(1))
        assessed = sv_assessed
      }
      from_meta = rep(FALSE, length(samples))
      column = if(!is.null(sv_from_metadata) && gene %in% names(sv_from_metadata))
        sv_from_metadata[[gene]] else NA
      if(!is.na(column)){
        if(column %in% colnames(these_samples_metadata)){
          status = as.character(these_samples_metadata[[column]][meta_rows])
          from_meta = status %in% sv_positive_values
          assessed = assessed | (!is.na(status) & status != "" & status != "NA")
        }else{
          warning("SV metadata column ", column, " not found")
        }
      }
      source = ifelse(from_sv & from_meta, "sv_data,metadata",
                      ifelse(from_sv, "sv_data", "metadata"))
      add_calls(sv$feature[i], from_sv | from_meta, assessed, source)
    }
  }

  # copy number -------------------------------------------------------------
  cnv = features[features$data_type == "cnv", , drop = FALSE]
  if(nrow(cnv) > 0){
    if(is.null(seg_data)){
      message(nrow(cnv), " copy number features need seg_data and are not assessed")
      for(f in cnv$feature) add_calls(f, rep(0, length(samples)), FALSE, "seg_data")
    }else{
      cn_assessed = available("cn_available")
      if(is.null(cn_assessed)) cn_assessed = samples %in% seg_data$ID
      # call each cnv_features set on its own so that its features are
      # measured as that set defines (e.g. DLBClass arm backgrounds)
      for(set in unique(cnv$cnv_feature_set)){
        these = cnv[cnv$cnv_feature_set == set, , drop = FALSE]
        set_features = get_cnv_features(genome_build, feature_set = set)
        set_features = set_features[set_features$name %in% these$cnv_feature, , drop = FALSE]
        cnv_call = utils::modifyList(cnv_args, list(
          seg_data = seg_data,
          cnv_features = set_features,
          these_samples_metadata = these_samples_metadata,
          genome_build = genome_build,
          output_format = "matrix"
        ))
        level = do.call(segmented_data_to_cnv_status, cnv_call)
        level = as.matrix(level)[samples, , drop = FALSE]
        for(i in seq_len(nrow(these))){
          add_calls(these$feature[i], level[, these$cnv_feature[i]],
                    cn_assessed, "seg_data")
        }
      }
    }
  }

  out = do.call(rbind, calls)
  out$data_type = features$data_type[match(out$feature, features$feature)]
  out = out[, c("sample_id", "feature", "data_type", "value", "assessed", "source")]
  rownames(out) = NULL
  out
}

#' @title Encode classifier features for one classifier
#'
#' @description Convert a feature table from [assemble_classifier_features]
#' into one classifier's features and encoding, e.g. the DLBClass GSM.
#'
#' @param feature_table Output of [assemble_classifier_features].
#' @param classifier Classifier name, e.g. "DLBClass" or "LymphGen".
#' @param output_format "matrix" (default; samples as rows, classifier
#' features as columns, in the classifier's order), "long" (one row per
#' sample and classifier feature, with an assessed column) or "gsm" (the
#' DLBClass layout: a classifier_name column, then one column per sample).
#' @param unassessed_value Value for a classifier feature none of whose
#' features was assessed in a sample. Default 0 (as DLBClass expects); NA
#' marks them.
#' @param genome_build Genome build. Default "grch37".
#' @param submodel Optional LymphGen submodel name(s), as for
#' [get_classifier_features], to encode only the features those submodels
#' use, e.g. "NoCGH" for samples without copy number. Default NULL encodes
#' every feature.
#'
#' @return A data frame in the requested format.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gsm = classifier_feature_matrix(feature_table, "DLBClass",
#'                                 output_format = "gsm")
#' }
classifier_feature_matrix = function(feature_table,
                                     classifier,
                                     output_format = "matrix",
                                     unassessed_value = 0,
                                     genome_build = "grch37",
                                     submodel = NULL){
  if(!output_format %in% c("matrix", "long", "gsm")){
    stop("output_format must be \"matrix\", \"long\" or \"gsm\"")
  }
  if(length(classifier) != 1) stop("give one classifier")
  feature_map = get_classifier_features(classifier, genome_build, submodel = submodel)
  samples = unique(feature_table$sample_id)
  # features the table doesn't have (e.g. one assembled from mutations only)
  # count as not assessed
  missing_features = setdiff(feature_map$feature, feature_table$feature)
  if(length(missing_features) > 0){
    message(length(missing_features), " ", classifier, " features are not in ",
            "feature_table and count as not assessed (e.g. ",
            paste(utils::head(missing_features, 3), collapse = ", "), ")")
    feature_table = rbind(
      feature_table[, c("sample_id", "feature", "value", "assessed")],
      data.frame(sample_id = rep(samples, length(missing_features)),
                 feature = rep(missing_features, each = length(samples)),
                 value = 0, assessed = FALSE))
  }
  classifier_names = unique(feature_map$classifier_feature)
  joined = merge(feature_map[, c("classifier_feature", "feature", "value", "min_level")],
                 feature_table[, c("sample_id", "feature", "value", "assessed")],
                 by = "feature", suffixes = c("_map", ""))
  min_level = ifelse(is.na(joined$min_level), 1, joined$min_level)
  present = !is.na(joined$value) & joined$value >= min_level
  joined$encoded = ifelse(!present, 0,
                          ifelse(joined$value_map == "level", joined$value,
                                 suppressWarnings(as.numeric(joined$value_map))))
  key = list(sample_id = joined$sample_id,
             classifier_feature = joined$classifier_feature)
  encoded = stats::aggregate(joined$encoded, key, max)
  assessed = stats::aggregate(joined$assessed, key, any)
  long = data.frame(sample_id = encoded$sample_id,
                    classifier_feature = encoded$classifier_feature,
                    value = encoded$x,
                    assessed = assessed$x[match(paste(encoded$sample_id, encoded$classifier_feature),
                                                paste(assessed$sample_id, assessed$classifier_feature))],
                    stringsAsFactors = FALSE)
  long$value[!long$assessed] = unassessed_value
  long = long[order(match(long$sample_id, samples),
                    match(long$classifier_feature, classifier_names)), , drop = FALSE]
  rownames(long) = NULL
  if(output_format == "long") return(long)

  mat = matrix(unassessed_value, nrow = length(samples), ncol = length(classifier_names),
               dimnames = list(samples, classifier_names))
  mat[cbind(match(long$sample_id, samples),
            match(long$classifier_feature, classifier_names))] = long$value
  if(output_format == "gsm"){
    out = data.frame(classifier_name = classifier_names, t(mat), check.names = FALSE)
    rownames(out) = NULL
    return(out)
  }
  as.data.frame(mat)
}
