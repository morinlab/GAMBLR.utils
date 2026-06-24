#' @title Test feature associations between two or more groups
#'
#' @description Run per-feature statistical tests comparing mutation status
#' across two or more groups. Accepts either a MAF data frame or a pre-built
#' sample-by-feature matrix. When \code{test = "fisher"} and more than two
#' \code{comparison_values} are supplied, all pairwise Fisher tests are run
#' automatically and the results include a \code{comparison} column identifying
#' each pair. For a single omnibus test across all groups, use
#' \code{test = "chi_square"}.
#'
#' @details
#' \strong{MAF input:} Each unique value in \code{maf_column} for a given gene
#' becomes one feature, and one test is run per feature. For example, with
#' \code{maf_column = "mutation_alias"} and genes \code{c("TP53", "EZH2")},
#' features such as \code{TP53_trunc}, \code{TP53_R175}, and \code{EZH2_Y641}
#' each receive their own contingency table and test. NA values in
#' \code{maf_column} are silently dropped.
#'
#' When a gene produces no features with at least \code{min_samples} mutated
#' samples under \code{maf_column}, the gene is retested using
#' \code{fallback_column} instead, and \code{used_fallback} is set to TRUE in
#' the output.
#'
#' \strong{Matrix input:} Each column of \code{feature_matrix} is tested as a
#' feature. Values greater than zero are treated as mutated (i.e. counts are
#' acceptable; they are coerced to binary internally). Samples present in the
#' matrix but absent from \code{metadata} are silently excluded. Features with
#' fewer than \code{min_samples} non-zero rows are skipped. The \code{gene}
#' column in the output is set equal to \code{feature} (the column name).
#'
#' @param maf A MAF-format data frame. Must contain at least
#' \code{Hugo_Symbol} and \code{Tumor_Sample_Barcode}. Mutually exclusive with
#' \code{feature_matrix}.
#' @param metadata A sample metadata data frame. Must contain
#' \code{Tumor_Sample_Barcode} and \code{comparison_column}.
#' @param comparison_column Name of the metadata column that defines groups.
#' @param feature_matrix Optional sample-by-feature matrix or data frame where
#' rows are samples and columns are features (e.g. genes or mutation categories).
#' Values should be numeric; anything greater than zero is treated as mutated.
#' Either supply a \code{sample_id_column} column or set
#' \code{sample_id_column = NULL} to use rownames. Mutually exclusive with
#' \code{maf}.
#' @param sample_id_column Name of the column in \code{feature_matrix} that
#' holds sample identifiers. Defaults to \code{"Tumor_Sample_Barcode"}. Set to
#' \code{NULL} to use rownames instead.
#' @param comparison_values Optional character vector of group values to
#' include. For \code{test = "fisher"} with exactly two values, one comparison
#' is run. With three or more values all C(k, 2) pairwise comparisons are run
#' and a \code{comparison} column is added to the output. Defaults to all
#' unique values in \code{comparison_column}.
#' @param genes Optional character vector of \code{Hugo_Symbol} values to
#' restrict testing to. Ignored (with a warning) when \code{feature_matrix} is
#' supplied — subset the matrix columns directly instead.
#' @param maf_column The MAF column whose unique values define features.
#' Default is \code{"Hugo_Symbol"} (one test per gene). Set to
#' \code{"mutation_alias"}, \code{"driver_alias"}, or \code{"hotspot_alias"}
#' to use sub-gene categories from
#' \code{\link[GAMBLR.utils]{annotate_curated_drivers}}. Setting
#' \code{"hot_spot"} derives labels as \code{<Hugo_Symbol>HOTSPOT} for TRUE
#' rows and \code{Hugo_Symbol} for FALSE rows. Ignored when
#' \code{feature_matrix} is supplied.
#' @param fallback_column MAF column to use per-gene when \code{maf_column}
#' yields no features with at least \code{min_samples} mutated samples.
#' Default is \code{"Hugo_Symbol"}. Ignored when \code{feature_matrix} is
#' supplied.
#' @param min_samples Minimum number of mutated samples required for a
#' feature to be tested. Default is 2.
#' @param test Statistical test to apply. One of \code{"fisher"} (default;
#' runs all pairwise comparisons when more than 2 comparison values are
#' supplied) or \code{"chi_square"} (single omnibus test across all groups).
#' @param contrast How to define comparisons. One of \code{"pairwise"}
#' (default; all C(k, 2) group pairs) or \code{"one_vs_rest"} (each group
#' tested against all remaining samples pooled together). In
#' \code{"one_vs_rest"} mode the output contains k comparisons labelled
#' \code{"<group> vs rest"} and the OR represents enrichment in \code{<group>}
#' relative to the pooled rest. This produces a single signed estimate per
#' feature per group with no directionality ambiguity, making the results
#' particularly suitable for heatmap visualisation.
#' @param p_adjust_method Multiple testing correction method passed to
#' \code{\link[stats]{p.adjust}}. Default is \code{"BH"}.
#' @param restrict_to_maf Logical. If TRUE (default), the denominator for each
#' group is restricted to samples present in both \code{maf} and
#' \code{metadata}. Set to FALSE to use all metadata samples regardless of
#' MAF presence; a warning is emitted because samples absent from the MAF are
#' treated as wild-type, which may be unsafe if they were not sequenced.
#' Ignored when \code{feature_matrix} is supplied.
#' @param verbose If TRUE, messages which genes fell back to
#' \code{fallback_column}. Ignored when \code{feature_matrix} is supplied.
#'
#' @return A tibble with one row per tested feature (or one row per feature
#' per pairwise comparison when more than two groups are supplied with
#' \code{test = "fisher"}). Always present: \code{gene} (Hugo_Symbol, or
#' equal to \code{feature} in matrix mode), \code{feature},
#' \code{used_fallback} (always FALSE in matrix mode), \code{p_value},
#' \code{q_value}, \code{n_mutated}.
#' Fisher-only columns: \code{OR}, \code{conf_low}, \code{conf_high}, and
#' per-group count columns \code{n_mutated_<group>} and denominator columns
#' \code{n_total_<group>} (e.g. \code{n_mutated_FL}, \code{n_total_FL}).
#' In pairwise or one-vs-rest mode a \code{comparison} column
#' (e.g. \code{"FL vs DLBCL"} or \code{"FL vs rest"}) identifies the
#' comparison; \code{n_mutated_<group>} columns cover all original groups.
#' Chi-square-only column: \code{statistic}.
#'
#' @import dplyr purrr
#' @importFrom stats fisher.test chisq.test p.adjust
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.open))
#'
#' meta = get_gambl_metadata()
#' fl_dlbcl_meta = dplyr::filter(meta, pathology %in% c("FL", "DLBCL","BL")) %>%
#'   GAMBLR.helpers::check_and_clean_metadata(duplicate_action = "keep_first")
#'
#' maf = get_all_coding_ssm(fl_dlbcl_meta)
#' maf = annotate_curated_drivers(
#'   maf,
#'   genes_of_interest = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2","DDX3X","MYC")
#' )
#'
#' # Fine-grained: one test per mutation alias category
#' results = test_feature_associations(
#'   maf = maf,
#'   metadata = fl_dlbcl_meta,
#'   comparison_column = "pathology",
#'   comparison_values = c("FL", "DLBCL","BL"),
#'   genes = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2","DDX3X","MYC"),
#'   maf_column = "mutation_alias"
#' )
#' dplyr::arrange(results, q_value)
#'
#' # With fallback: use hotspot_alias where available, Hugo_Symbol otherwise
#' results_fb = test_feature_associations(
#'   maf = maf,
#'   metadata = fl_dlbcl_meta,
#'   comparison_column = "pathology",
#'   comparison_values = c("FL", "DLBCL"),
#'   genes = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2"),
#'   maf_column = "hotspot_alias",
#'   fallback_column = "mutation_alias",
#'   verbose = TRUE
#' )
#'
#' # Pre-built feature matrix via assemble_genetic_features
#' # assemble_genetic_features returns a sample x feature matrix with rownames
#' # as sample IDs and integer-encoded mutation status (0 = absent, >0 = present).
#' coding_genes = c("EZH2", "CREBBP", "KMT2D", "TP53", "MYD88",
#'                  "NOTCH1", "NOTCH2", "DDX3X", "MYC", "BCL2")
#' nc_genes = intersect(coding_genes,
#'                      unique(GAMBLR.data::grch37_ashm_regions$gene))
#'
#' raw_maf = get_all_coding_ssm(fl_dlbcl_meta)
#'
#' # Create a GENE_alias column for each distinct annotation category (e.g. "mutation_alias", "hotspot_alias", "driver_alias")
#' feat_mat  = summarise_mutation_status(
#'    annotated_maf
#'  )
#' head(colnames(feat_mat))
#'
#' # sample IDs are rownames — pass sample_id_column = NULL
#' results_mat = test_feature_associations(
#'   feature_matrix     = feat_mat,
#'   sample_id_column   = NULL,
#'   metadata           = fl_dlbcl_meta,
#'   comparison_column  = "pathology",
#'   comparison_values  = c("FL", "DLBCL", "BL")
#' )
#' dplyr::arrange(results_mat, q_value)
#' }
#'
test_feature_associations = function(maf = NULL,
                                     metadata,
                                     comparison_column,
                                     feature_matrix = NULL,
                                     sample_id_column = "Tumor_Sample_Barcode",
                                     comparison_values = NULL,
                                     genes = NULL,
                                     maf_column = "Hugo_Symbol",
                                     fallback_column = "Hugo_Symbol",
                                     min_samples = 2,
                                     test = "fisher",
                                     contrast = "pairwise",
                                     p_adjust_method = "BH",
                                     restrict_to_maf = TRUE,
                                     verbose = FALSE) {

  test     = match.arg(test,     c("fisher", "chi_square"))
  contrast = match.arg(contrast, c("pairwise", "one_vs_rest"))

  if (is.null(maf) == is.null(feature_matrix)) {
    stop("Provide exactly one of 'maf' or 'feature_matrix'.")
  }
  matrix_mode = !is.null(feature_matrix)

  if (matrix_mode && !is.null(genes)) {
    warning("'genes' is ignored when 'feature_matrix' is supplied. ",
            "Subset the matrix columns directly to restrict features.")
  }

  if (is.null(comparison_values)) {
    if ("factor" %in% class(metadata[[comparison_column]])) {
      comparison_values = levels(metadata[[comparison_column]])
    } else {
      comparison_values = unique(metadata[[comparison_column]])
    }
  }

  pairwise_mode = test == "fisher" &&
    (length(comparison_values) > 2 || contrast == "one_vs_rest")
  if (test == "fisher" && length(comparison_values) < 2) {
    stop("Fisher's exact test requires at least 2 comparison_values.")
  }

  metadata = dplyr::filter(metadata, .data[[comparison_column]] %in% comparison_values)
  metadata$.group = factor(metadata[[comparison_column]], levels = comparison_values)

  # shared sample universe (drives denominators and contingency tables)
  all_meta_samples = dplyr::select(metadata, Tumor_Sample_Barcode, .group)

  if (!matrix_mode) {
    if ("maf_data" %in% class(maf)) {
      maf = strip_genomic_classes(maf)
    }
    maf = dplyr::filter(maf, Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)

    if (is.null(genes)) {
      genes = unique(maf$Hugo_Symbol)
    } else {
      maf = dplyr::filter(maf, Hugo_Symbol %in% genes)
    }

    if (restrict_to_maf) {
      maf_samples = unique(maf$Tumor_Sample_Barcode)
      all_samples = dplyr::filter(all_meta_samples,
                                  Tumor_Sample_Barcode %in% maf_samples)
    } else {
      warning(
        "restrict_to_maf = FALSE: all metadata samples are used as the ",
        "denominator, including those absent from the MAF. Samples without MAF ",
        "entries are treated as wild-type, which may be unsafe if they were ",
        "not sequenced."
      )
      all_samples = all_meta_samples
    }
  } else {
    # matrix mode: use all metadata samples as denominator
    all_samples = all_meta_samples
  }

  group_levels = levels(metadata$.group)
  group_totals = table(all_samples$.group)

  # ── shared test helpers ────────────────────────────────────────────────────

  run_test = function(feature_label, mutated_ids) {
    sample_status = all_samples %>%
      dplyr::mutate(is_mutated = as.integer(Tumor_Sample_Barcode %in% mutated_ids))

    n_mut = sum(sample_status$is_mutated)
    if (n_mut < min_samples) return(NULL)

    ct = table(sample_status$is_mutated, sample_status$.group)
    for (lvl in c("0", "1")) {
      if (!lvl %in% rownames(ct)) {
        extra = matrix(0L, nrow = 1, ncol = ncol(ct),
                       dimnames = list(lvl, colnames(ct)))
        ct = rbind(ct, extra)
      }
    }
    ct = ct[c("1", "0"), , drop = FALSE]

    per_group = sample_status %>%
      dplyr::group_by(.group) %>%
      dplyr::summarise(n_mut = sum(is_mutated), .groups = "drop")

    if (test == "fisher") {
      res = tryCatch(fisher.test(ct), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      row = tibble::tibble(
        OR        = unname(res$estimate),
        conf_low  = res$conf.int[1],
        conf_high = res$conf.int[2],
        p_value   = res$p.value,
        n_mutated = n_mut
      )
      for (grp in group_levels) {
        row[[paste0("n_mutated_", grp)]] =
          per_group$n_mut[per_group$.group == grp]
        row[[paste0("n_total_", grp)]] = as.integer(group_totals[grp])
      }
    } else {
      res = tryCatch(chisq.test(ct), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      row = tibble::tibble(
        statistic = unname(res$statistic),
        p_value   = res$p.value,
        n_mutated = n_mut
      )
      for (grp in group_levels) {
        row[[paste0("n_total_", grp)]] = as.integer(group_totals[grp])
      }
    }
    row
  }

  run_test_pair = function(feature_label, mutated_ids, pair) {
    pair_samples = dplyr::filter(all_samples, .group %in% pair)
    pair_samples = dplyr::mutate(pair_samples,
      is_mutated = as.integer(Tumor_Sample_Barcode %in% mutated_ids))

    n_mut_pair = sum(pair_samples$is_mutated)
    if (n_mut_pair < min_samples) return(NULL)

    ct = table(pair_samples$is_mutated, pair_samples$.group)
    for (lvl in c("0", "1")) {
      if (!lvl %in% rownames(ct)) {
        extra = matrix(0L, nrow = 1, ncol = ncol(ct),
                       dimnames = list(lvl, colnames(ct)))
        ct = rbind(ct, extra)
      }
    }
    ct = ct[c("1", "0"), pair, drop = FALSE]

    res = tryCatch(fisher.test(ct), error = function(e) NULL)
    if (is.null(res)) return(NULL)

    full_status = dplyr::mutate(all_samples,
      is_mutated = as.integer(Tumor_Sample_Barcode %in% mutated_ids))
    per_group = full_status %>%
      dplyr::group_by(.group) %>%
      dplyr::summarise(n_mut = sum(is_mutated), .groups = "drop")

    row = tibble::tibble(
      OR        = unname(res$estimate),
      conf_low  = res$conf.int[1],
      conf_high = res$conf.int[2],
      p_value   = res$p.value,
      n_mutated = n_mut_pair
    )
    for (grp in group_levels) {
      row[[paste0("n_mutated_", grp)]] =
        per_group$n_mut[per_group$.group == grp]
      row[[paste0("n_total_", grp)]] = as.integer(group_totals[grp])
    }
    row
  }

  run_test_ovr = function(feature_label, mutated_ids, g) {
    ovr_samples = all_samples %>%
      dplyr::mutate(
        .group_bin = factor(
          ifelse(as.character(.group) == g, g, "rest"),
          levels = c(g, "rest")
        ),
        is_mutated = as.integer(Tumor_Sample_Barcode %in% mutated_ids)
      )

    n_mut_ovr = sum(ovr_samples$is_mutated)
    if (n_mut_ovr < min_samples) return(NULL)

    ct = table(ovr_samples$is_mutated, ovr_samples$.group_bin)
    for (lvl in c("0", "1")) {
      if (!lvl %in% rownames(ct)) {
        extra = matrix(0L, nrow = 1, ncol = ncol(ct),
                       dimnames = list(lvl, colnames(ct)))
        ct = rbind(ct, extra)
      }
    }
    ct = ct[c("1", "0"), c(g, "rest"), drop = FALSE]

    res = tryCatch(fisher.test(ct), error = function(e) NULL)
    if (is.null(res)) return(NULL)

    per_group = all_samples %>%
      dplyr::mutate(is_mutated = as.integer(Tumor_Sample_Barcode %in% mutated_ids)) %>%
      dplyr::group_by(.group) %>%
      dplyr::summarise(n_mut = sum(is_mutated), .groups = "drop")

    row = tibble::tibble(
      OR        = unname(res$estimate),
      conf_low  = res$conf.int[1],
      conf_high = res$conf.int[2],
      p_value   = res$p.value,
      n_mutated = n_mut_ovr
    )
    for (grp in group_levels) {
      row[[paste0("n_mutated_", grp)]] =
        per_group$n_mut[per_group$.group == grp]
      row[[paste0("n_total_", grp)]] = as.integer(group_totals[grp])
    }
    row
  }

  # ── MAF path ──────────────────────────────────────────────────────────────

  if (!matrix_mode) {
    derive_feature = function(m, col) {
      if (col == "hot_spot") {
        if (!"hot_spot" %in% colnames(m))
          stop("Column 'hot_spot' not found in maf. Run annotate_curated_drivers() first.")
        return(ifelse(m$hot_spot %in% c(TRUE, "TRUE"),
                      paste0(m$Hugo_Symbol, "HOTSPOT"),
                      m$Hugo_Symbol))
      }
      if (col == "Hugo_Symbol") return(m$Hugo_Symbol)
      if (!col %in% colnames(m))
        stop("Column '", col, "' not found in maf. Run annotate_curated_drivers() first.")
      m[[col]]
    }

    n_testable = function(gene_maf, col) {
      feat = derive_feature(gene_maf, col)
      valid = !is.na(feat) & nchar(as.character(feat)) > 0
      if (!any(valid)) return(0L)
      gene_maf$.feat = feat
      gene_maf[valid, ] %>%
        dplyr::filter(Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode) %>%
        dplyr::distinct(Tumor_Sample_Barcode, .feat) %>%
        dplyr::count(.feat) %>%
        dplyr::filter(n >= min_samples) %>%
        nrow()
    }

    run_gene = function(g, pair = NULL, group_ovr = NULL, verbose_fallback = FALSE) {
      gene_maf = dplyr::filter(maf, Hugo_Symbol == g)
      if (nrow(gene_maf) == 0) return(NULL)

      used_fallback = FALSE
      active_col = maf_column

      if (n_testable(gene_maf, maf_column) == 0 && maf_column != fallback_column) {
        if (verbose_fallback) {
          message(g, ": no testable features under '", maf_column,
                  "', falling back to '", fallback_column, "'")
        }
        active_col = fallback_column
        used_fallback = TRUE
      }

      feat = derive_feature(gene_maf, active_col)
      valid = !is.na(feat) & nchar(as.character(feat)) > 0
      gene_maf$.feat = feat
      gene_maf = gene_maf[valid, ]
      if (nrow(gene_maf) == 0) return(NULL)

      features = unique(gene_maf$.feat)

      purrr::map_dfr(features, function(f) {
        mutated_ids = gene_maf$Tumor_Sample_Barcode[gene_maf$.feat == f]
        row = if (!is.null(pair))      run_test_pair(f, mutated_ids, pair) else
              if (!is.null(group_ovr)) run_test_ovr(f, mutated_ids, group_ovr) else
                                        run_test(f, mutated_ids)
        if (is.null(row)) return(NULL)
        dplyr::bind_cols(
          tibble::tibble(gene = g, feature = f, used_fallback = used_fallback),
          row
        )
      })
    }

    if (contrast == "one_vs_rest") {
      results = dplyr::bind_rows(lapply(comparison_values, function(g) {
        res = dplyr::bind_rows(lapply(genes, run_gene, group_ovr = g))
        if (nrow(res) > 0) res$comparison = paste(g, "vs rest")
        res
      }))
    } else if (pairwise_mode) {
      pairs = combn(comparison_values, 2, simplify = FALSE)

      if (verbose && maf_column != fallback_column) {
        fallen = Filter(function(g) {
          gm = dplyr::filter(maf, Hugo_Symbol == g)
          nrow(gm) > 0 && n_testable(gm, maf_column) == 0
        }, genes)
        if (length(fallen) > 0) {
          message(paste(fallen, collapse = ", "),
                  ": no testable features under '", maf_column,
                  "', falling back to '", fallback_column, "'")
        }
      }

      results = dplyr::bind_rows(lapply(pairs, function(pair) {
        res = dplyr::bind_rows(lapply(genes, run_gene, pair = pair))
        if (nrow(res) > 0) res$comparison = paste(pair, collapse = " vs ")
        res
      }))
    } else {
      results = dplyr::bind_rows(lapply(genes, run_gene, verbose_fallback = verbose))
    }

  # ── matrix path ───────────────────────────────────────────────────────────

  } else {
    feature_matrix = as.data.frame(feature_matrix)

    if (!is.null(sample_id_column)) {
      if (!sample_id_column %in% colnames(feature_matrix)) {
        stop("'", sample_id_column, "' not found in feature_matrix. ",
             "Set sample_id_column = NULL to use rownames instead.")
      }
      mat_ids  = feature_matrix[[sample_id_column]]
      feat_cols = setdiff(colnames(feature_matrix), sample_id_column)
    } else {
      if (is.null(rownames(feature_matrix)) ||
          identical(rownames(feature_matrix), as.character(seq_len(nrow(feature_matrix))))) {
        stop("feature_matrix has no meaningful rownames. ",
             "Add rownames or supply a sample_id_column.")
      }
      mat_ids   = rownames(feature_matrix)
      feat_cols = colnames(feature_matrix)
    }

    # restrict to intersection of matrix rows and metadata (silently in both directions)
    keep      = mat_ids %in% all_samples$Tumor_Sample_Barcode
    mat_ids   = mat_ids[keep]
    feature_matrix = feature_matrix[keep, , drop = FALSE]
    # update all_samples/group_totals so denominators only cover matrix-represented samples
    all_samples  = dplyr::filter(all_samples, Tumor_Sample_Barcode %in% mat_ids)
    group_totals = table(all_samples$.group)

    run_matrix_feature = function(f, pair = NULL, group_ovr = NULL) {
      mutated_ids = mat_ids[feature_matrix[[f]] > 0]
      row = if (!is.null(pair))      run_test_pair(f, mutated_ids, pair) else
            if (!is.null(group_ovr)) run_test_ovr(f, mutated_ids, group_ovr) else
                                      run_test(f, mutated_ids)
      if (is.null(row)) return(NULL)
      dplyr::bind_cols(
        tibble::tibble(gene = sub("_.*$", "", f), feature = f, used_fallback = FALSE),
        row
      )
    }

    if (contrast == "one_vs_rest") {
      results = dplyr::bind_rows(lapply(comparison_values, function(g) {
        res = dplyr::bind_rows(lapply(feat_cols, run_matrix_feature, group_ovr = g))
        if (nrow(res) > 0) res$comparison = paste(g, "vs rest")
        res
      }))
    } else if (pairwise_mode) {
      pairs = combn(comparison_values, 2, simplify = FALSE)
      results = dplyr::bind_rows(lapply(pairs, function(pair) {
        res = dplyr::bind_rows(lapply(feat_cols, run_matrix_feature, pair = pair))
        if (nrow(res) > 0) res$comparison = paste(pair, collapse = " vs ")
        res
      }))
    } else {
      results = dplyr::bind_rows(lapply(feat_cols, run_matrix_feature))
    }
  }

  # ── shared post-processing ────────────────────────────────────────────────

  if (nrow(results) == 0) {
    warning("No testable features found. Try reducing min_samples or checking maf_column.")
    return(results)
  }

  results = dplyr::mutate(results,
                          q_value = stats::p.adjust(p_value, method = p_adjust_method))
  results
}


#' @title Select Informative Features Using Elastic Net
#'
#' @description Uses elastic net regularization to identify non-redundant,
#' informative features from a sample-by-feature matrix with respect to group
#' membership. Addresses the correlated/hierarchical feature problem common in
#' mutation matrices (e.g. \code{GENE_coding}, \code{GENE_driver}, and
#' \code{GENE_any} all measure similar biology). Features are penalized
#' according to specificity: catch-all columns (\code{_any}) receive the
#' highest penalty, intermediate summaries (\code{_coding}, \code{_driver},
#' \code{_other}, \code{_noncoding}) receive a moderate penalty, and specific
#' aliases (e.g. \code{EZH2_Y641}) receive the baseline penalty. This biases
#' the model toward selecting the most granular informative representation of
#' each gene.
#'
#' @details One \code{cv.glmnet} model (binomial family) is fit per comparison.
#' Features with a non-zero coefficient at the selected lambda are retained.
#' Returned odds ratios (\code{OR}) are exponentiated regularized logistic
#' regression coefficients, not Fisher test ORs, and reflect the direction and
#' relative magnitude of each feature's contribution after adjusting for the
#' others. Set \code{run_fisher = TRUE} to also return Fisher-based ORs,
#' confidence intervals, and adjusted p-values for use with
#' \code{plot_feature_associations}.
#'
#' @param feature_matrix A data frame with samples as rows and binary features
#'   as columns. Either include a sample ID column (named by
#'   \code{sample_id_column}) or set \code{sample_id_column = NULL} and supply
#'   meaningful row names.
#' @param metadata A data frame with at least a sample ID column and
#'   \code{comparison_column}.
#' @param comparison_column Column name in \code{metadata} defining groups.
#' @param sample_id_column Column in \code{feature_matrix} containing sample
#'   IDs. Set to \code{NULL} to use row names instead. Default
#'   \code{"Tumor_Sample_Barcode"}.
#' @param comparison_values Character vector of group values to include.
#'   Defaults to all unique values in \code{comparison_column}.
#' @param contrast One of \code{"one_vs_rest"} or \code{"pairwise"}. Default
#'   \code{"one_vs_rest"}.
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = LASSO). Values
#'   around 0.5–0.8 are recommended for correlated features. Default \code{0.8}.
#' @param lambda_select Which cross-validated lambda to use: \code{"1se"}
#'   (sparser, more regularized) or \code{"min"} (minimum CV error, more
#'   features). Default \code{"1se"}.
#' @param nfolds Number of cross-validation folds. Automatically capped at the
#'   smaller class size for each comparison. Default \code{10}.
#' @param run_fisher Logical. If \code{TRUE}, runs Fisher's exact test post-hoc
#'   on selected features and adds \code{OR}, \code{conf_low},
#'   \code{conf_high}, \code{p_value}, and \code{q_value} columns that are
#'   compatible with \code{plot_feature_associations}. Default \code{FALSE}.
#' @param p_adjust_method P-value adjustment method passed to
#'   \code{\link[stats]{p.adjust}}. Used only when \code{run_fisher = TRUE}.
#'   Default \code{"BH"}.
#' @param min_samples Minimum number of samples with a non-zero value required
#'   to include a feature. Default \code{2}.
#' @param positive_only Logical. If \code{TRUE}, restricts the returned features
#'   to those with a positive coefficient (\code{coef > 0}) in at least one
#'   comparison. Features whose coefficients are negative in every comparison
#'   are informative only as depletions and may reflect collinearity rather than
#'   group-specific enrichment. Default \code{FALSE}.
#' @param seed Random seed for reproducibility. Default \code{42}.
#' @param balance_classes Logical. If \code{TRUE}, per-sample observation
#'   weights are applied inside each \code{cv.glmnet} call to compensate for
#'   class imbalance. Positive-class samples receive weight
#'   \code{n_negative / n_positive} and negative-class samples receive weight
#'   1, so both classes contribute equally to the loss regardless of their
#'   sizes. Recommended when one class dominates the training data (e.g. DLBCL
#'   outnumbering minority classes several-fold). Default \code{FALSE}.
#' @param priority_features Optional character vector of feature names that
#'   should always enter the model regardless of regularisation. These features
#'   receive a \code{penalty.factor} of 0 in \code{glmnet}, making them
#'   unpenalised — lambda shrinks all other coefficients but leaves these
#'   unconstrained. Features not found in the filtered feature matrix are
#'   ignored with a warning. In the pairwise contrast, priority features are
#'   also exempt from the \code{min_samples} variance filter. Default
#'   \code{NULL}.
#' @param return_models Logical. If \code{TRUE}, returns a list with elements
#'   \code{results} (the usual tibble), \code{models} (a named list of fitted
#'   model objects, one per comparison), \code{contrast}, and
#'   \code{comparison_values}. Each element of \code{models} contains
#'   \code{fit} (the \code{cv.glmnet} object), \code{lambda} (the selected
#'   lambda), \code{features} (character vector of feature names in the order
#'   the model was trained on, after collinear-column dropping), \code{label}
#'   (the comparison label string), and \code{pos_class} (the group encoded as
#'   1 in the binary outcome). Pass the returned list to
#'   \code{GAMBLR.predict::predict_from_glmnet_models()} to classify new
#'   samples. Default \code{FALSE}.
#'
#' @return When \code{return_models = FALSE} (default): a tibble. When
#'   \code{run_fisher = FALSE}: columns \code{gene}, \code{feature},
#'   \code{comparison}, \code{coef}, \code{OR} (from \code{exp(coef)}).
#'   When \code{run_fisher = TRUE}: additionally \code{OR} (Fisher),
#'   \code{conf_low}, \code{conf_high}, \code{p_value}, \code{q_value},
#'   \code{n_mutated}, and per-group count columns — matching
#'   \code{test_feature_associations} output for use with
#'   \code{plot_feature_associations}. When \code{return_models = TRUE}: a
#'   named list with \code{results} (the tibble described above),
#'   \code{models}, \code{contrast}, and \code{comparison_values}.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' mat = get_binary_matrix(
#'   these_samples_metadata = get_gambl_metadata(),
#'   maf_data = get_coding_ssm()
#' )
#' meta = get_gambl_metadata()
#' selected = select_informative_features(
#'   feature_matrix    = mat,
#'   metadata          = meta,
#'   comparison_column = "pathology",
#'   comparison_values = c("DLBCL", "FL", "BL"),
#'   run_fisher        = TRUE
#' )
#' plot_feature_associations(selected, max_q = 0.1)
#' }
select_informative_features = function(
  feature_matrix,
  metadata,
  comparison_column,
  sample_id_column  = "Tumor_Sample_Barcode",
  comparison_values = NULL,
  contrast          = "one_vs_rest",
  alpha             = 0.8,
  lambda_select     = "1se",
  nfolds            = 10,
  run_fisher        = FALSE,
  p_adjust_method   = "BH",
  min_samples       = 2,
  positive_only     = FALSE,
  seed              = 42,
  return_models     = FALSE,
  priority_features = NULL,
  balance_classes   = FALSE
) {
  contrast      = match.arg(contrast, c("one_vs_rest", "pairwise"))
  lambda_select = match.arg(lambda_select, c("1se", "min"))

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required. Install with: install.packages('glmnet')")
  }

  # ── resolve sample IDs and feature columns ────────────────────────────────

  feature_matrix = as.data.frame(feature_matrix)

  if (!is.null(sample_id_column)) {
    if (!sample_id_column %in% colnames(feature_matrix)) {
      stop("'", sample_id_column, "' not found in feature_matrix.")
    }
    mat_ids   = feature_matrix[[sample_id_column]]
    feat_cols = setdiff(colnames(feature_matrix), sample_id_column)
  } else {
    mat_ids   = rownames(feature_matrix)
    feat_cols = colnames(feature_matrix)
    if (is.null(mat_ids) ||
        identical(mat_ids, as.character(seq_len(nrow(feature_matrix))))) {
      stop("feature_matrix has no meaningful row names. ",
           "Add row names or supply sample_id_column.")
    }
  }

  # ── align metadata ────────────────────────────────────────────────────────

  meta_id_col = if (!is.null(sample_id_column)) sample_id_column else
                  "Tumor_Sample_Barcode"

  if (!meta_id_col %in% colnames(metadata)) {
    stop("'", meta_id_col, "' not found in metadata.")
  }

  if (is.null(comparison_values)) {
    comparison_values = sort(unique(as.character(
      metadata[[comparison_column]]
    )))
  }

  meta_sub = dplyr::filter(
    metadata,
    .data[[comparison_column]] %in% comparison_values,
    .data[[meta_id_col]] %in% mat_ids
  )
  meta_sub$.group = as.character(meta_sub[[comparison_column]])
  keep_ids = meta_sub[[meta_id_col]]

  # align feature matrix rows to metadata (same order as keep_ids)
  row_idx  = match(keep_ids, mat_ids)
  feat_mat = as.matrix(feature_matrix[row_idx, feat_cols, drop = FALSE])
  rownames(feat_mat) = keep_ids
  feat_mat[is.na(feat_mat)] = 0

  # ── filter low-prevalence features ───────────────────────────────────────

  n_nonzero = colSums(feat_mat > 0)
  feat_mat  = feat_mat[, n_nonzero >= min_samples |
                           colnames(feat_mat) %in% priority_features,
                       drop = FALSE]
  feat_cols_kept = colnames(feat_mat)

  if (ncol(feat_mat) == 0) {
    stop("No features remain after min_samples = ", min_samples, " filtering.")
  }

  # ── penalty factors (lower = less penalized = more likely selected) ───────
  # _any columns are the most redundant → highest penalty
  # intermediate summaries (_coding, _driver, _other, _noncoding) → moderate
  # specific aliases (e.g. EZH2_Y641) → baseline

  penalty_factors = dplyr::case_when(
    grepl("_any$", feat_cols_kept)                                   ~ 3.0,
    grepl("_coding$|_driver$|_hotspot$|_other$|_noncoding$", feat_cols_kept)  ~ 1.5,
    TRUE                                                             ~ 1.0
  )

  if (!is.null(priority_features)) {
    missing_pf = setdiff(priority_features, feat_cols_kept)
    if (length(missing_pf) > 0) {
      warning(
        "priority_features not found in feature matrix (ignored): ",
        paste(missing_pf, collapse = ", ")
      )
    }
    penalty_factors[feat_cols_kept %in% priority_features] = 0
  }

  group_vec = meta_sub$.group

  # ── inner helper: fit one glmnet model and extract non-zero features ──────

  fit_one = function(x, y, label, pf, feat_names, pos_class = NULL) {
    pos_n = sum(y == 1L)
    neg_n = sum(y == 0L)
    if (pos_n < 2L || neg_n < 2L) {
      warning("Skipping comparison '", label, "': too few samples in one class.")
      return(NULL)
    }

    # ── drop perfectly collinear columns ─────────────────────────────────────
    # Identical columns cause numerical instability in glmnet and arise
    # naturally when a gene has only one hotspot alias (GENE_hotspot ==
    # GENE_alias) or all coding mutations are drivers (GENE_coding ==
    # GENE_driver). Within each group of identical columns, keep the one with
    # the lowest penalty factor (most specific feature); break ties by position.
    col_sig    = apply(x, 2, paste, collapse = "\r")
    dup_groups = split(seq_along(col_sig), col_sig)
    keep_idx   = sort(vapply(dup_groups, function(idx) idx[which.min(pf[idx])],
                             integer(1L)))
    if (length(keep_idx) < ncol(x)) {
      n_dropped = ncol(x) - length(keep_idx)
      message(n_dropped, " perfectly collinear feature(s) removed for '",
              label, "'.")
      x          = x[, keep_idx, drop = FALSE]
      pf         = pf[keep_idx]
      feat_names = feat_names[keep_idx]
    }

    # ── within-gene correlation-aware penalty scaling ─────────────────────────
    # Summary features (e.g. GENE_hotspot, GENE_coding) are by construction
    # supersets of more-specific features of the same gene, creating partial
    # collinearity. For each gene with multiple features, detect which features
    # are binary supersets of others via crossprod, then scale each parent's
    # penalty upward by (1 + max |cor| with any child feature).
    gene_labels = sub("_[^_]+$", "", feat_names)
    for (g in unique(gene_labels)) {
      idx = which(gene_labels == g)
      if (length(idx) < 2L) next
      sub_x = x[, idx, drop = FALSE]
      cs    = colSums(sub_x)
      cp    = crossprod(sub_x)          # cp[j,k] = dot(col_j, col_k)
      # is_child[j,k] = TRUE when col_j is a strict binary subset of col_k
      # (every 1 in col_j also appears in col_k, but col_k has more 1s)
      cs_mat   = matrix(cs, length(idx), length(idx))  # row j = cs[j]
      is_child = (cs_mat == cp) & (cs_mat < matrix(cs, length(idx), length(idx),
                                                    byrow = TRUE))
      cor_mat  = suppressWarnings(stats::cor(sub_x))
      cor_mat[is.na(cor_mat)] = 0
      diag(cor_mat) = 0
      # for each parent k, take the max |cor| with any child j
      max_child_r  = apply(is_child * abs(cor_mat), 2, max)
      pf[idx]      = pf[idx] * (1 + max_child_r)
    }

    safe_folds  = min(nfolds, pos_n, neg_n)
    obs_weights = if (balance_classes) {
      ifelse(y == 1L, neg_n / pos_n, 1.0)
    } else {
      NULL
    }

    set.seed(seed)
    fit = tryCatch(
      glmnet::cv.glmnet(
        x              = x,
        y              = y,
        family         = "binomial",
        alpha          = alpha,
        nfolds         = safe_folds,
        penalty.factor = pf,
        standardize    = TRUE,
        weights        = obs_weights
      ),
      error = function(e) {
        warning("cv.glmnet failed for '", label, "': ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) return(NULL)

    lam      = if (lambda_select == "1se") fit$lambda.1se else fit$lambda.min
    coefs    = glmnet::coef.glmnet(fit, s = lam)
    coef_vec = as.numeric(coefs)[-1]    # drop intercept; length == ncol(x)
    names(coef_vec) = feat_names

    nonzero = which(coef_vec != 0)
    if (length(nonzero) == 0) {
      message("No features selected for '", label,
              "'. Try lambda_select = 'min' or lower alpha.")
      return(NULL)
    }

    coefs = tibble::tibble(
      gene       = sub("_.*$", "", feat_names[nonzero]),
      feature    = feat_names[nonzero],
      comparison = label,
      coef       = coef_vec[nonzero],
      OR_glmnet  = exp(coef_vec[nonzero])
    )

    if (!return_models) return(coefs)

    list(
      coefs = coefs,
      model = list(
        fit       = fit,
        lambda    = lam,
        features  = feat_names,
        label     = label,
        pos_class = pos_class
      )
    )
  }

  # ── run comparisons ───────────────────────────────────────────────────────

  model_list = list()

  if (contrast == "one_vs_rest") {
    raw = lapply(comparison_values, function(g) {
      y = as.integer(group_vec == g)
      fit_one(feat_mat, y, paste(g, "vs rest"),
              penalty_factors, feat_cols_kept, pos_class = g)
    })
    if (return_models) {
      results    = dplyr::bind_rows(lapply(raw, `[[`, "coefs"))
      model_list = stats::setNames(
        lapply(raw, `[[`, "model"),
        sapply(raw, function(r) if (is.null(r)) NA_character_ else r$model$label)
      )
      model_list = model_list[!vapply(model_list, is.null, logical(1L))]
    } else {
      results = dplyr::bind_rows(raw)
    }
  } else {
    pairs = combn(comparison_values, 2, simplify = FALSE)
    raw   = lapply(pairs, function(pair) {
      idx   = group_vec %in% pair
      y     = as.integer(group_vec[idx] == pair[[1]])
      x     = feat_mat[idx, , drop = FALSE]
      label = paste(pair, collapse = " vs ")

      keep  = colSums(x > 0) >= min_samples |
              feat_cols_kept %in% priority_features
      if (!any(keep)) return(NULL)

      fit_one(x[, keep, drop = FALSE], y, label,
              penalty_factors[keep], feat_cols_kept[keep],
              pos_class = pair[[1]])
    })
    if (return_models) {
      results    = dplyr::bind_rows(lapply(raw, `[[`, "coefs"))
      model_list = stats::setNames(
        lapply(raw, `[[`, "model"),
        sapply(raw, function(r) if (is.null(r)) NA_character_ else r$model$label)
      )
      model_list = model_list[!vapply(model_list, is.null, logical(1L))]
    } else {
      results = dplyr::bind_rows(raw)
    }
  }

  if (positive_only && !is.null(results) && nrow(results) > 0) {
    keep_features = results |>
      dplyr::group_by(.data$feature) |>
      dplyr::summarise(any_positive = any(.data$coef > 0), .groups = "drop") |>
      dplyr::filter(.data$any_positive) |>
      dplyr::pull(.data$feature)
    results = dplyr::filter(results, .data$feature %in% keep_features)
  }

  if (is.null(results) || nrow(results) == 0) {
    warning("No features selected across any comparison.")
    empty = tibble::tibble(
      gene       = character(),
      feature    = character(),
      comparison = character(),
      coef       = numeric(),
      OR_glmnet  = numeric()
    )
    if (run_fisher) {
      empty$OR        = numeric()
      empty$conf_low  = numeric()
      empty$conf_high = numeric()
      empty$p_value   = numeric()
      empty$q_value   = numeric()
    }
    if (return_models) {
      return(list(
        results           = empty,
        models            = model_list,
        contrast          = contrast,
        comparison_values = comparison_values
      ))
    }
    return(empty)
  }

  # ── optional post-hoc Fisher test ─────────────────────────────────────────

  if (run_fisher) {
    selected_feats = unique(results$feature)

    # build feature_matrix subset: sample_id_column + selected features only
    sub_cols = if (!is.null(sample_id_column)) {
      c(sample_id_column, intersect(selected_feats, colnames(feature_matrix)))
    } else {
      intersect(selected_feats, colnames(feature_matrix))
    }
    sub_mat = feature_matrix[, sub_cols, drop = FALSE]
    if (!is.null(sample_id_column) && is.null(rownames(sub_mat))) {
      rownames(sub_mat) = seq_len(nrow(sub_mat))
    }

    fisher_res = test_feature_associations(
      feature_matrix    = sub_mat,
      metadata          = metadata,
      comparison_column = comparison_column,
      sample_id_column  = sample_id_column,
      comparison_values = comparison_values,
      contrast          = contrast,
      test              = "fisher",
      p_adjust_method   = p_adjust_method,
      min_samples       = min_samples
    )

    # join Fisher columns onto glmnet results; glmnet OR kept as OR_glmnet
    join_cols = intersect(
      c("feature", "comparison", "OR", "conf_low", "conf_high",
        "p_value", "q_value", "n_mutated",
        grep("^n_mutated_|^n_total_", colnames(fisher_res), value = TRUE)),
      colnames(fisher_res)
    )
    results = dplyr::left_join(
      results,
      dplyr::select(fisher_res, dplyr::all_of(join_cols)),
      by = c("feature", "comparison")
    )
  }

  if (return_models) {
    return(list(
      results           = results,
      models            = model_list,
      contrast          = contrast,
      comparison_values = comparison_values
    ))
  }

  results
}
