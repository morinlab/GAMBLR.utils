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
