#' @title Test feature associations between two or more groups
#'
#' @description Run per-feature statistical tests comparing mutation status
#' across two or more groups. When \code{test = "fisher"} and more than two
#' \code{comparison_values} are supplied, all pairwise Fisher tests are run
#' automatically and the results include a \code{comparison} column identifying
#' each pair. For a single omnibus test across all groups, use
#' \code{test = "chi_square"}.
#'
#' @details Each unique value in \code{maf_column} for a given gene becomes
#' one feature, and one test is run per feature. For example, with
#' \code{maf_column = "mutation_alias"} and genes \code{c("TP53", "EZH2")},
#' features such as \code{TP53_trunc}, \code{TP53_R175}, and \code{EZH2_Y641}
#' each receive their own contingency table and test. NA values in
#' \code{maf_column} are silently dropped (e.g. \code{hotspot_alias} is NA
#' for non-hotspot mutations).
#'
#' When a gene produces no features with at least \code{min_samples} mutated
#' samples under \code{maf_column} (e.g. TP53 has no \code{hotspot_alias}
#' entries because all its drivers are truncations), the gene is retested
#' using \code{fallback_column} instead, and \code{used_fallback} is set to
#' TRUE in the output.
#'
#' @param maf A MAF-format data frame. Must contain at least
#' \code{Hugo_Symbol} and \code{Tumor_Sample_Barcode}.
#' @param metadata A sample metadata data frame. Must contain
#' \code{Tumor_Sample_Barcode} and \code{comparison_column}.
#' @param comparison_column Name of the metadata column that defines groups.
#' @param comparison_values Optional character vector of group values to
#' include. For \code{test = "fisher"} with exactly two values, one comparison
#' is run. With three or more values all C(k, 2) pairwise comparisons are run
#' and a \code{comparison} column is added to the output. Defaults to all
#' unique values in \code{comparison_column}.
#' @param genes Optional character vector of \code{Hugo_Symbol} values to
#' restrict testing to. Defaults to all genes present in \code{maf}.
#' @param maf_column The MAF column whose unique values define features.
#' Default is \code{"Hugo_Symbol"} (one test per gene). Set to
#' \code{"mutation_alias"}, \code{"driver_alias"}, or \code{"hotspot_alias"}
#' to use sub-gene categories from
#' \code{\link[GAMBLR.utils]{annotate_curated_drivers}}. Setting
#' \code{"hot_spot"} derives labels as \code{<Hugo_Symbol>HOTSPOT} for TRUE
#' rows and \code{Hugo_Symbol} for FALSE rows.
#' @param fallback_column MAF column to use per-gene when \code{maf_column}
#' yields no features with at least \code{min_samples} mutated samples.
#' Default is \code{"Hugo_Symbol"}.
#' @param min_samples Minimum number of mutated samples required for a
#' feature to be tested. Default is 2.
#' @param test Statistical test to apply. One of \code{"fisher"} (default;
#' runs all pairwise comparisons when more than 2 comparison values are
#' supplied) or \code{"chi_square"} (single omnibus test across all groups).
#' @param p_adjust_method Multiple testing correction method passed to
#' \code{\link[stats]{p.adjust}}. Default is \code{"BH"}.
#' @param restrict_to_maf Logical. If TRUE (default), the denominator for each
#' group is restricted to samples present in both \code{maf} and
#' \code{metadata}. Set to FALSE to use all metadata samples regardless of
#' MAF presence; a warning is emitted because samples absent from the MAF are
#' treated as wild-type, which may be unsafe if they were not sequenced.
#' @param verbose If TRUE, messages which genes fell back to
#' \code{fallback_column}.
#'
#' @return A tibble with one row per tested feature (or one row per feature
#' per pairwise comparison when more than two groups are supplied with
#' \code{test = "fisher"}). Always present: \code{gene} (Hugo_Symbol),
#' \code{feature} (value from \code{maf_column} or \code{fallback_column}),
#' \code{used_fallback} (logical), \code{p_value}, \code{q_value},
#' \code{n_mutated}.
#' Fisher-only columns: \code{OR}, \code{conf_low}, \code{conf_high}, and
#' per-group count columns \code{n_mutated_<group>} and denominator columns
#' \code{n_total_<group>} (e.g. \code{n_mutated_FL}, \code{n_total_FL}).
#' In pairwise mode a \code{comparison} column (e.g. \code{"FL vs DLBCL"})
#' identifies which pair each row belongs to; \code{n_mutated_<group>} columns
#' cover all groups, not only the two being compared in that row.
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
#' # (e.g. TP53 falls back because it has no hotspot_alias entries)
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
#' dplyr::filter(results_fb, used_fallback)
#' }
#'
test_feature_associations = function(maf,
                                     metadata,
                                     comparison_column,
                                     comparison_values = NULL,
                                     genes = NULL,
                                     maf_column = "Hugo_Symbol",
                                     fallback_column = "Hugo_Symbol",
                                     min_samples = 2,
                                     test = "fisher",
                                     p_adjust_method = "BH",
                                     restrict_to_maf = TRUE,
                                     verbose = FALSE) {

  test = match.arg(test, c("fisher", "chi_square"))

  if (is.null(comparison_values)) {
    if ("factor" %in% class(metadata[[comparison_column]])) {
      comparison_values = levels(metadata[[comparison_column]])
    } else {
      comparison_values = unique(metadata[[comparison_column]])
    }
  }

  pairwise_mode <- test == "fisher" && length(comparison_values) > 2
  if (test == "fisher" && length(comparison_values) < 2) {
    stop("Fisher's exact test requires at least 2 comparison_values.")
  }

  metadata = dplyr::filter(metadata, .data[[comparison_column]] %in% comparison_values)
  metadata$.group = factor(metadata[[comparison_column]], levels = comparison_values)

  if ("maf_data" %in% class(maf)) {
    maf = strip_genomic_classes(maf)
  }
  maf = dplyr::filter(maf, Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)

  if (is.null(genes)) {
    genes = unique(maf$Hugo_Symbol)
  } else {
    maf = dplyr::filter(maf, Hugo_Symbol %in% genes)
  }

  all_meta_samples = dplyr::select(metadata, Tumor_Sample_Barcode, .group)

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

  group_levels = levels(metadata$.group)
  group_totals = table(all_samples$.group)

  # derive feature labels from a given column, with hot_spot special handling
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

  # count features with enough mutated samples under a given column
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

  # run one test for a single feature column vector against the comparison
  run_test = function(feature_label, mutated_ids) {
    sample_status = all_samples %>%
      dplyr::mutate(is_mutated = as.integer(Tumor_Sample_Barcode %in% mutated_ids))

    n_mut = sum(sample_status$is_mutated)
    if (n_mut < min_samples) return(NULL)

    ct = table(sample_status$is_mutated, sample_status$.group)

    # ensure both 0 and 1 rows are present
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

  # variant of run_test used in pairwise mode:
  # Fisher test is restricted to `pair` samples; n_mutated_* covers all k groups
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

  run_gene = function(g, pair = NULL, verbose_fallback = FALSE) {
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
      row = if (is.null(pair)) run_test(f, mutated_ids) else
                                run_test_pair(f, mutated_ids, pair)
      if (is.null(row)) return(NULL)
      dplyr::bind_cols(
        tibble::tibble(gene = g, feature = f, used_fallback = used_fallback),
        row
      )
    })
  }

  if (pairwise_mode) {
    pairs = combn(comparison_values, 2, simplify = FALSE)

    # report fallbacks once (before looping over pairs)
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

  if (nrow(results) == 0) {
    warning("No testable features found. Try reducing min_samples or checking maf_column.")
    return(results)
  }

  results = dplyr::mutate(results,
                          q_value = stats::p.adjust(p_value, method = p_adjust_method))
  results
}
