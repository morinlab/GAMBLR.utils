utils::globalVariables(c(
  # column names from LymphGen input files
  "ENTREZ.ID", "Type", "Sample", "Sample.ID", "Location",
  # columns created inside dplyr verbs
  "Hugo_Symbol", "direction", "col_name", "score",
  "feature_col", "is_named_hotspot",
  # columns in the result tibble used inside mutate
  "gene_requested", "status", "suggested_symbol",
  # columns from annotate_curated_drivers output
  "driver_alias"
))

#' @title Check gene symbols against a LymphGen feature matrix
#'
#' @description Validate a set of HGNC gene symbols against the genes
#' available in a \code{lymphgen_to_dlbclone()} feature matrix, and
#' identify any outdated symbols using the \pkg{org.Hs.eg.db} alias table.
#'
#' @details
#' For each requested gene the function reports one of three statuses:
#'
#' * **found** - the symbol matches at least one mutation-related column
#'   in `feature_matrix`. Columns are matched by the pattern
#'   `^GENE(_|$)`, so `MYD88_hotspot`, `MYD88_other`, and plain `MYD88`
#'   all count as a match for the gene `MYD88`.
#' * **alias** - the symbol is not a current HGNC name but is a known
#'   alias; a suggested replacement is returned
#'   (e.g. `HIST1H1C` -> `H1-2`).
#' * **unknown** - the symbol was not found as either a current HGNC
#'   name or any known alias in \pkg{org.Hs.eg.db}.
#'
#' The function always prints a one-line summary and, for `alias` and
#' `unknown` cases, a detailed breakdown. The return value is invisibly
#' returned so the function can be used for its side-effect messages or
#' captured for programmatic use.
#'
#' @param genes Character vector of HGNC gene symbols to validate.
#' @param feature_matrix Data frame returned by
#'   \code{lymphgen_to_dlbclone()}. Mutation-related columns are those
#'   lacking a `_gain` or `_loss` suffix and not named `sample_id`.
#'   A gene is considered present if any column matches `^GENE(_|$)`.
#'
#' @return A data frame (invisibly) with columns:
#'   \describe{
#'     \item{gene_requested}{Symbol as supplied by the user.}
#'     \item{status}{\code{"found"}, \code{"alias"}, or
#'       \code{"unknown"}.}
#'     \item{suggested_symbol}{Current HGNC symbol when
#'       `status == "alias"`; \code{NA} otherwise. Note: the suggested
#'       symbol may itself not be present in `feature_matrix` if no
#'       samples in the dataset carried that gene's mutation.}
#'   }
#'
#' @importFrom AnnotationDbi mapIds
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- lymphgen_to_dlbclone(
#'   mutation_flat_path     = "lymphgen_inputs/mutation_flat.txt",
#'   cnv_path               = "lymphgen_inputs/cnv_flat.txt",
#'   sample_annotation_path = "lymphgen_inputs/sample_annotation.txt"
#' )
#'
#' check_lymphgen_genes(
#'   genes          = c("MYD88", "HIST1H1C", "NOTCH1", "FAKEGENE"),
#'   feature_matrix = mat
#' )
#' }
#'
check_lymphgen_genes <- function(genes, feature_matrix) {
  stopifnot(
    "`genes` must be a character vector" =
      is.character(genes),
    "`feature_matrix` must be a data frame" =
      is.data.frame(feature_matrix),
    "`feature_matrix` must contain 'sample_id'" =
      "sample_id" %in% colnames(feature_matrix)
  )

  # Mutation-related columns: exclude CNV (_gain/_loss) and sample_id.
  # SV (_SV) cols are intentionally kept here because a gene like BCL2
  # could be "found" via BCL2_SV even if it has no mutation column.
  mut_cols <- colnames(feature_matrix)
  mut_cols <- mut_cols[!grepl("_gain$|_loss$", mut_cols)]
  mut_cols <- setdiff(mut_cols, "sample_id")

  # A gene is "found" if ANY column matches ^GENE(_|$),
  # catching plain GENE, GENE_hotspot, GENE_other, GENE_noncoding, GENE_SV
  gene_found <- vapply(genes, function(g) {
    any(grepl(paste0("^", g, "(_|$)"), mut_cols))
  }, logical(1))

  found   <- genes[gene_found]
  missing <- genes[!gene_found]

  result <- dplyr::tibble(
    gene_requested   = genes,
    status           = dplyr::if_else(
      gene_found, "found", "pending"
    ),
    suggested_symbol = NA_character_
  )

  if (length(missing) > 0) {
    # Step 1: resolve each missing symbol as a known alias -> ENTREZ ID
    entrez_by_alias <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys      = missing,
        column    = "ENTREZID",
        keytype   = "ALIAS",
        multiVals = "first"
      )
    )

    # Step 2: for resolved aliases, look up the current HGNC symbol
    alias_to_entrez <- entrez_by_alias[!is.na(entrez_by_alias)]

    if (length(alias_to_entrez) > 0) {
      entrez_to_symbol <- suppressMessages(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys      = unname(alias_to_entrez),
          column    = "SYMBOL",
          keytype   = "ENTREZID",
          multiVals = "first"
        )
      )
      # alias_map: old symbol -> current symbol
      alias_map <- setNames(
        entrez_to_symbol[alias_to_entrez],
        names(alias_to_entrez)
      )
    } else {
      alias_map <- character(0)
    }

    result <- result |>
      dplyr::mutate(
        status = dplyr::case_when(
          status == "found" ~ "found",
          gene_requested %in% names(alias_map) ~ "alias",
          TRUE ~ "unknown"
        ),
        suggested_symbol = dplyr::case_when(
          status == "alias" ~
            unname(alias_map[gene_requested]),
          TRUE ~ NA_character_
        )
      )
  } else {
    result <- dplyr::mutate(
      result,
      status = dplyr::if_else(
        status == "pending", "unknown", status
      )
    )
  }

  # ── Summary messages ──────────────────────────────────────────────────
  n_found   <- sum(result$status == "found")
  n_alias   <- sum(result$status == "alias")
  n_unknown <- sum(result$status == "unknown")

  message(sprintf(
    "%d gene(s) found | %d outdated symbol(s) | %d unknown",
    n_found, n_alias, n_unknown
  ))

  if (n_alias > 0) {
    aliases <- dplyr::filter(result, status == "alias")
    message(
      "Outdated symbols - use the suggested current name instead:"
    )
    for (i in seq_len(nrow(aliases))) {
      in_mat <- any(grepl(
        paste0("^", aliases$suggested_symbol[i], "(_|$)"),
        mut_cols
      ))
      note <- if (in_mat) {
        " [available in matrix]"
      } else {
        " [not in matrix - gene may have no observations]"
      }
      message(sprintf(
        "  %-20s ->  %s%s",
        aliases$gene_requested[i],
        aliases$suggested_symbol[i],
        note
      ))
    }
  }

  if (n_unknown > 0) {
    unknowns <- dplyr::filter(result, status == "unknown")
    message(
      "Unknown symbols (not found as current HGNC name or alias):"
    )
    message("  ", paste(unknowns$gene_requested, collapse = ", "))
  }

  invisible(result)
}


# ──────────────────────────────────────────────────────────────────────────────

#' @title Convert LymphGen inputs to a feature matrix compatible with DLBCLone
#'
#' @description Read the standard LymphGen input files (mutation flat
#' file, copy-number file, and sample annotation file) and assemble them
#' into a wide sample-by-feature numeric matrix compatible with DLBCLone.
#'
#' @details
#' Four data types are extracted and combined:
#'
#' ## Mutations
#' Column naming mirrors the logic of \code{assemble_genetic_features()}:
#'
#' * **Hotspot genes** (`genes_hotspot`, default `"MYD88"`) produce two
#'   mutation columns per gene:
#'   - `<GENE>_hotspot` - named hotspot mutations (any LymphGen
#'     `Type` that is not `MUTATION`, `TRUNC`, or `Synon`, e.g. `L265P`).
#'   - `<GENE>_other` - coding mutations outside the hotspot
#'     (`MUTATION` or `TRUNC`).
#' * **Non-coding / aSHM genes** (`genes_noncoding`, default `NULL`)
#'   produce an additional `<GENE>_noncoding` column for `Synon`
#'   mutations. These genes also retain their regular mutation column(s).
#' * **All other genes** produce a single `<GENE>` column. Only
#'   non-synonymous mutations (`MUTATION`, `TRUNC`, named hotspots)
#'   contribute to this column. `Synon` mutations are included only for
#'   genes explicitly listed in `genes_synonymous`; for every other gene
#'   a sample with only a `Synon` record scores 0.
#'
#' ## Per-gene mutation type filtering
#' `per_gene_type_filter` is a named list that restricts which LymphGen
#' `Type` values are kept for specific genes before any column assignment.
#' This mirrors LymphGen's own per-gene rules — e.g. NOTCH1/NOTCH2 count
#' only truncating mutations. `Synon` is always exempt from this filter
#' (it is governed independently by `genes_noncoding` / `genes_synonymous`).
#' Genes absent from the list are unrestricted.
#'
#' ## Copy-number variants
#' Each gene contributes two columns: `<GENE>_gain` (GAIN / AMP) and
#' `<GENE>_loss` (HETLOSS / HOMDEL). Default gain encoding: GAIN = 1,
#' AMP = 2. Default loss encoding: HETLOSS = 1, HOMDEL = 2.
#'
#' ## Structural variants / translocations
#' Any column ending in `.transloc` in the sample annotation file is
#' mapped to a `<GENE>_SV` column (e.g. `BCL2.transloc` -> `BCL2_SV`).
#'
#' ## ENTREZ ID mapping
#' All ENTREZ IDs are converted to HGNC gene symbols via
#' \pkg{org.Hs.eg.db}. Genes that cannot be mapped are silently dropped.
#'
#' When `genes` is provided the output is further subset. Any symbol not
#' found in the matrix is automatically checked via
#' \code{\link{check_lymphgen_genes}} so outdated names (e.g.
#' `HIST1H1C`) trigger a helpful suggestion.
#'
#' @param mutation_flat_path Path to the LymphGen mutation flat file
#'   (tab-separated; columns: Sample, ENTREZ.ID, Type, Location).
#' @param cnv_path Path to the LymphGen copy-number file
#'   (tab-separated; columns: Sample, ENTREZ.ID, Type).
#' @param sample_annotation_path Path to the LymphGen sample annotation
#'   file (tab-separated; must contain `Sample.ID`; translocation columns
#'   are detected automatically by the `.transloc` suffix).
#' @param mutation_gene_list_path Optional path to the LymphGen gene list
#'   file (single column of ENTREZ IDs, one per row, with a header). When
#'   supplied with `restrict_to_gene_list = TRUE`, only genes from this
#'   list appear in the output.
#' @param genes Optional character vector of HGNC gene symbols. When
#'   provided, the output matrix is subset to all columns belonging to
#'   those genes (mutation, hotspot, other, noncoding, CNV, and SV
#'   columns alike). Any symbol not found is checked via
#'   \code{\link{check_lymphgen_genes}}, which messages the user with
#'   alias suggestions for outdated names.
#' @param genes_hotspot Character vector of genes whose mutations should
#'   be split into `<GENE>_hotspot` (named hotspot types, e.g. L265P)
#'   and `<GENE>_other` (MUTATION / TRUNC) columns instead of a single
#'   `<GENE>` column. Default: `c("MYD88")`.
#' @param genes_noncoding Optional character vector of genes for which
#'   `Synon` mutations are tracked in a dedicated `<GENE>_noncoding`
#'   column (modelling aSHM signal). For genes also in `genes_hotspot`,
#'   noncoding tracking is applied on top of the hotspot/other split.
#'   Default: `NULL` (no noncoding columns).
#' @param genes_synonymous Optional character vector of genes for which
#'   `Synon` mutations are included in the plain `<GENE>` column (scored
#'   with `mutation_encoding$synonymous`). For genes **not** in this list
#'   (and not in `genes_noncoding`), a sample with only a `Synon`
#'   mutation scores 0, matching the behaviour of
#'   \code{assemble_genetic_features()} where synonymous signal is only
#'   carried for genes explicitly listed in `genes_noncoding`. This
#'   parameter has no effect on genes in `genes_hotspot` (those never
#'   produce a plain `<GENE>` column). Default: `NULL`.
#' @param per_gene_type_filter Named list restricting which LymphGen
#'   `Type` values are counted for specific genes. Each name is an HGNC
#'   gene symbol; each value is a character vector of the `Type` strings
#'   that are **kept** (e.g. `c("TRUNC")` to accept truncating mutations
#'   only). `Synon` rows are always exempt — they are governed by
#'   `genes_noncoding` / `genes_synonymous` and are never removed by this
#'   filter. Genes not present in the list are unrestricted. To disable
#'   all defaults pass `list()`. Default encodes the two clearest
#'   LymphGen rules:
#'   \code{list(NOTCH1 = c("TRUNC"), NOTCH2 = c("TRUNC"))}.
#' @param mutation_encoding Named list controlling the numeric score for
#'   each LymphGen mutation class. Recognised keys: `none`, `synonymous`,
#'   `coding`, `truncating`, `hotspot`. Default:
#'   `list(none = 0, synonymous = 1, coding = 2, truncating = 2,
#'   hotspot = 2)`.
#' @param cnv_gain_encoding Named list controlling scores for copy gains.
#'   Recognised keys: `none`, `gain`, `amp`.
#'   Default: `list(none = 0, gain = 1, amp = 2)`.
#' @param cnv_loss_encoding Named list controlling scores for copy losses.
#'   Recognised keys: `none`, `hetloss`, `homdel`.
#'   Default: `list(none = 0, hetloss = 1, homdel = 2)`.
#' @param sv_encoding Named list controlling scores for translocation
#'   status. Recognised keys: `negative`, `positive`.
#'   Default: `list(negative = 0, positive = 1)`.
#' @param restrict_to_gene_list Logical. When `TRUE` (default) and
#'   `mutation_gene_list_path` is provided, output columns are restricted
#'   to genes present in that list.
#' @param include_cnv Logical. When `TRUE` (default), CNV gain/loss
#'   columns are appended to the output.
#' @param include_sv Logical. When `TRUE` (default), SV/translocation
#'   columns from the sample annotation file are appended to the output.
#' @param genome_build Reference genome build string used by the
#'   coordinate-based driver filter and passed to
#'   \code{\link{get_driver_coordinates}}. Supported values:
#'   `"grch37"` (default), `"hg19"`, `"hs37d5"`, `"GRCh37"`,
#'   `"hg38"`, `"grch38"`, `"GRCh38"`. Ignored when
#'   `genes_coordinate_filter` is `NULL` or empty.
#' @param genes_coordinate_filter Optional character vector of HGNC gene
#'   symbols to filter by curated driver coordinates via
#'   \code{\link{annotate_curated_drivers}}. Only mutations that fall
#'   inside the defined driver region for each gene **and** match the
#'   expected variant class are retained; all others are dropped.
#'   `Synon` rows are always exempt (governed by `genes_noncoding` /
#'   `genes_synonymous`). Genes listed in `genes_hotspot` are silently
#'   skipped (hotspot splitting already handles them). For genes handled
#'   by this filter, any entry in `per_gene_type_filter` is automatically
#'   skipped to avoid double-filtering. A warning is issued for any gene
#'   not present in the curated coordinate table returned by
#'   \code{\link{get_driver_coordinates}}. Set to `NULL` or
#'   `character(0)` to disable coordinate filtering entirely and rely on
#'   `per_gene_type_filter` alone. Default: `c("NOTCH1", "NOTCH2")`,
#'   which mirrors the coordinate restrictions LymphGen itself applies.
#' @param verbose Logical. When `TRUE`, progress messages are printed.
#'
#' @return A data frame with one row per sample and one column per
#'   genomic feature. The first column is `sample_id`. Mutation columns
#'   are named by HGNC symbol with optional `_hotspot`, `_other`, or
#'   `_noncoding` suffixes; CNV columns carry `_gain` or `_loss`; SV
#'   columns carry `_SV`. All values are numeric and reflect the supplied
#'   encoding lists.
#'
#' @import dplyr tidyr readr
#' @importFrom AnnotationDbi mapIds
#' @export
#'
#' @examples
#' \dontrun{
#' # Full matrix; MYD88 is split into MYD88_hotspot + MYD88_other
#' feature_matrix <- lymphgen_to_dlbclone(
#'   mutation_flat_path      = "lymphgen_inputs/mutation_flat.txt",
#'   cnv_path                = "lymphgen_inputs/cnv_flat.txt",
#'   sample_annotation_path  = "lymphgen_inputs/sample_annotation.txt",
#'   mutation_gene_list_path = "lymphgen_inputs/mutation_gene_list.txt"
#' )
#'
#' # Subset to a specific gene panel; includes _hotspot/_other for MYD88
#' dlbcl_panel <- lymphgen_to_dlbclone(
#'   mutation_flat_path      = "lymphgen_inputs/mutation_flat.txt",
#'   cnv_path                = "lymphgen_inputs/cnv_flat.txt",
#'   sample_annotation_path  = "lymphgen_inputs/sample_annotation.txt",
#'   mutation_gene_list_path = "lymphgen_inputs/mutation_gene_list.txt",
#'   genes = c("MYD88", "CD79B", "EZH2", "CREBBP", "BCL2", "TP53"),
#'   genes_noncoding = c("MYD88", "CREBBP", "PIM1", "SOCS1")
#' )
#' }
#'
lymphgen_to_dlbclone <- function(
    mutation_flat_path,
    cnv_path,
    sample_annotation_path,
    mutation_gene_list_path = NULL,
    genes                   = NULL,
    genes_hotspot           = c("MYD88"),
    genes_noncoding         = NULL,
    genes_synonymous        = NULL,
    per_gene_type_filter    = list(
      NOTCH1 = c("TRUNC"),
      NOTCH2 = c("TRUNC")
    ),
    mutation_encoding = list(
      none = 0, synonymous = 1, coding = 2,
      truncating = 2, hotspot = 2
    ),
    cnv_gain_encoding = list(none = 0, gain = 1, amp = 2),
    cnv_loss_encoding = list(none = 0, hetloss = 1, homdel = 2),
    sv_encoding       = list(negative = 0, positive = 1),
    restrict_to_gene_list   = TRUE,
    include_cnv             = TRUE,
    include_sv              = TRUE,
    genome_build            = "grch37",
    genes_coordinate_filter = c("NOTCH1", "NOTCH2"),
    verbose                 = FALSE
) {

  # ── 1. Read input files ──────────────────────────────────────────────
  if (verbose) message("Reading input files...")

  mutations <- readr::read_tsv(
    mutation_flat_path,
    col_types      = readr::cols(),
    show_col_types = FALSE
  )
  cnv_data <- readr::read_tsv(
    cnv_path,
    col_types      = readr::cols(),
    show_col_types = FALSE
  )
  sample_annot <- readr::read_tsv(
    sample_annotation_path,
    col_types      = readr::cols(),
    show_col_types = FALSE
  )

  stopifnot(
    "mutation_flat_path must have columns: Sample, ENTREZ.ID, Type" =
      all(c("Sample", "ENTREZ.ID", "Type") %in% colnames(mutations)),
    "cnv_path must have columns: Sample, ENTREZ.ID, Type" =
      all(c("Sample", "ENTREZ.ID", "Type") %in% colnames(cnv_data)),
    "sample_annotation_path must contain a 'Sample.ID' column" =
      "Sample.ID" %in% colnames(sample_annot)
  )

  all_samples <- as.character(sample_annot$Sample.ID)

  # Coerce to character(0) if NULL for safe %in% checks below
  genes_hotspot_vec    <- genes_hotspot    %||% character(0)
  genes_noncoding_vec  <- genes_noncoding  %||% character(0)
  genes_synonymous_vec <- genes_synonymous %||% character(0)

  # ── 2. Optional gene-list restriction ───────────────────────────────
  gene_list_entrez <- NULL
  if (!is.null(mutation_gene_list_path) && restrict_to_gene_list) {
    if (verbose) message("Reading mutation gene list...")
    gene_list_entrez <- readr::read_tsv(
      mutation_gene_list_path,
      col_names      = TRUE,
      col_types      = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    ) |>
      dplyr::pull(1)
  }

  # ── 3. ENTREZ ID -> HGNC symbol mapping ─────────────────────────────
  all_entrez <- unique(c(
    as.character(mutations$ENTREZ.ID),
    as.character(cnv_data$ENTREZ.ID)
  ))

  if (!is.null(gene_list_entrez)) {
    all_entrez <- intersect(all_entrez, gene_list_entrez)
  }

  if (verbose) {
    message(
      "Mapping ", length(all_entrez), " ENTREZ IDs to HGNC symbols..."
    )
  }

  symbol_map <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys      = all_entrez,
      column    = "SYMBOL",
      keytype   = "ENTREZID",
      multiVals = "first"
    )
  )
  symbol_map <- symbol_map[!is.na(symbol_map)]

  if (verbose) {
    message(sprintf(
      "  Mapped: %d  |  Unmapped (dropped): %d",
      length(symbol_map),
      length(all_entrez) - length(symbol_map)
    ))
  }

  # Annotate a data frame with Hugo_Symbol from ENTREZ.ID
  add_symbol <- function(df) {
    df |>
      dplyr::mutate(ENTREZ.ID = as.character(ENTREZ.ID)) |>
      dplyr::filter(ENTREZ.ID %in% names(symbol_map)) |>
      dplyr::mutate(Hugo_Symbol = unname(symbol_map[ENTREZ.ID]))
  }

  # ── 4. Mutation matrix ───────────────────────────────────────────────
  # LymphGen mutation types that indicate a named hotspot (e.g. L265P):
  # anything that is NOT one of the three standard generic types.
  std_types <- c("MUTATION", "TRUNC", "Synon")

  if (verbose) message("Processing mutations...")

  mutations_annotated <- mutations |> add_symbol()

  # ── 4a. Coordinate-based driver filter setup ─────────────────────────
  # genes_coordinate_filter names the genes to be filtered by coordinate
  # range. Hotspot genes are silently excluded (their splitting already
  # handles them). If NULL/empty, coordinate filtering is disabled.
  coord_filter_genes <- character(0)
  driver_coords      <- NULL

  genes_coord_vec <- setdiff(
    genes_coordinate_filter %||% character(0),
    genes_hotspot_vec
  )

  if (length(genes_coord_vec) > 0) {
    driver_coords <- tryCatch(
      get_driver_coordinates(genome_build),
      error = function(e) {
        if (verbose) {
          message(
            "  Could not load driver coordinates: ", e$message,
            " -- coordinate filter disabled."
          )
        }
        NULL
      }
    )
    if (!is.null(driver_coords)) {
      supported       <- unique(driver_coords$Hugo_Symbol)
      not_in_table    <- setdiff(genes_coord_vec, supported)
      if (length(not_in_table) > 0) {
        warning(
          "genes_coordinate_filter: gene(s) not found in the curated ",
          "coordinate table and will be skipped: ",
          paste(not_in_table, collapse = ", ")
        )
      }
      coord_filter_genes <- intersect(genes_coord_vec, supported)
      if (verbose && length(coord_filter_genes) > 0) {
        message(
          "  Coordinate filter will apply to: ",
          paste(coord_filter_genes, collapse = ", ")
        )
      }
    }
  }

  # ── 4b. Per-gene type filter (skips coordinate-handled genes) ────────
  # Build the effective filter: genes already governed by the coordinate
  # filter are excluded so they are not double-filtered.
  effective_type_filter <- per_gene_type_filter[
    !names(per_gene_type_filter) %in% coord_filter_genes
  ]

  if (length(effective_type_filter) > 0) {
    # Build a long table of gene + allowed (non-Synon) types
    allowed_long   <- dplyr::bind_rows(lapply(
      names(effective_type_filter),
      function(g) dplyr::tibble(
        Hugo_Symbol = g,
        Type        = effective_type_filter[[g]]
      )
    ))
    filtered_genes <- names(effective_type_filter)

    if (verbose) {
      rules <- paste(
        vapply(filtered_genes, function(g) {
          paste0(
            g, ": ",
            paste(effective_type_filter[[g]], collapse = "/")
          )
        }, character(1)),
        collapse = ", "
      )
      message("  Per-gene type filter applied: ", rules)
    }

    mutations_annotated <- dplyr::bind_rows(
      # Genes with no filter: keep all rows unchanged
      dplyr::filter(
        mutations_annotated, !Hugo_Symbol %in% filtered_genes
      ),
      # Filtered genes: Synon rows always pass through;
      # non-Synon rows kept only when Type is in the allowed list
      dplyr::filter(
        mutations_annotated,
        Hugo_Symbol %in% filtered_genes, Type == "Synon"
      ),
      dplyr::semi_join(
        dplyr::filter(
          mutations_annotated,
          Hugo_Symbol %in% filtered_genes, Type != "Synon"
        ),
        allowed_long,
        by = c("Hugo_Symbol", "Type")
      )
    )
  }

  # ── 4c. Coordinate-based driver filter (apply) ───────────────────────
  # For each gene in coord_filter_genes, build a pseudo-MAF and run
  # annotate_curated_drivers() to keep only mutations that fall inside
  # the curated driver region with a matching variant class.
  # Synon rows are always exempt (governed by genes_noncoding /
  # genes_synonymous).
  if (length(coord_filter_genes) > 0) {
    coord_mutations <- dplyr::filter(
      mutations_annotated,
      Hugo_Symbol %in% coord_filter_genes
    )

    if (nrow(coord_mutations) > 0) {
      active_coord_genes <- intersect(
        coord_filter_genes,
        unique(coord_mutations$Hugo_Symbol)
      )

      if (verbose) {
        message(
          "  Applying coordinate driver filter to: ",
          paste(active_coord_genes, collapse = ", ")
        )
      }

      # Look up chromosome for each ENTREZ ID (CHR column is deprecated
      # in org.Hs.eg.db but still functional; suppress warning)
      coord_entrez <- unique(as.character(coord_mutations$ENTREZ.ID))
      chr_map <- suppressWarnings(suppressMessages(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys      = coord_entrez,
          column    = "CHR",
          keytype   = "ENTREZID",
          multiVals = "first"
        )
      ))

      # genome_build local alias avoids ambiguity inside mutate
      .gb <- genome_build

      # Build pseudo-MAF; embed genome_build as a column so
      # annotate_curated_drivers() can auto-detect it (the function
      # has a known issue where explicit genome_build skips setting
      # annotated_maf internally)
      pseudo_maf <- coord_mutations |>
        dplyr::mutate(
          .lymphgen_row = dplyr::row_number(),
          genome_build  = .gb,
          Chromosome    = unname(chr_map[as.character(ENTREZ.ID)]),
          Start_Position = suppressWarnings(as.integer(Location)),
          End_Position   = suppressWarnings(as.integer(Location)),
          Variant_Classification = dplyr::case_when(
            Type == "MUTATION" ~ "Missense_Mutation",
            Type == "TRUNC"    ~ "Nonsense_Mutation",
            Type == "Synon"    ~ "Silent",
            TRUE               ~ "Missense_Mutation"
          ),
          Tumor_Sample_Barcode = Sample
        )

      # annotate_curated_drivers auto-detects genome_build from the
      # column; do NOT pass genome_build explicitly to avoid the bug
      # where annotated_maf is left unset
      annotated_pseudo <- suppressMessages(
        annotate_curated_drivers(
          pseudo_maf,
          genes_of_interest = active_coord_genes
        )
      )

      # driver_alias == GENE_other → outside driver region → drop
      # driver_alias == GENE_driver → confirmed driver → keep
      # Synon is always exempt regardless of driver_alias
      driver_verdict <- dplyr::select(
        annotated_pseudo,
        dplyr::all_of(c(".lymphgen_row", "driver_alias"))
      )

      coord_filtered <- coord_mutations |>
        dplyr::mutate(.lymphgen_row = dplyr::row_number()) |>
        dplyr::left_join(driver_verdict, by = ".lymphgen_row") |>
        dplyr::filter(
          Type == "Synon" |
            driver_alias != paste0(Hugo_Symbol, "_other")
        ) |>
        dplyr::select(
          -dplyr::all_of(c(".lymphgen_row", "driver_alias"))
        )

      if (verbose) {
        n_before <- nrow(coord_mutations)
        n_after  <- nrow(coord_filtered)
        message(sprintf(
          "  Coordinate filter: %d -> %d mutations (%d removed)",
          n_before, n_after, n_before - n_after
        ))
      }

      # Reassemble: non-coordinate genes + filtered coordinate genes
      mutations_annotated <- dplyr::bind_rows(
        dplyr::filter(
          mutations_annotated, !Hugo_Symbol %in% coord_filter_genes
        ),
        coord_filtered
      )
    }
  }

  mut_matrix <- mutations_annotated |>
    dplyr::mutate(
      is_named_hotspot = !Type %in% std_types,
      # Numeric score
      score = dplyr::case_when(
        Type == "Synon"    ~ mutation_encoding$synonymous,
        Type == "MUTATION" ~ mutation_encoding$coding,
        Type == "TRUNC"    ~ mutation_encoding$truncating,
        is_named_hotspot   ~ mutation_encoding$hotspot,
        TRUE               ~ mutation_encoding$none
      ),
      # Destination column name
      feature_col = dplyr::case_when(
        # Hotspot gene + named hotspot type -> _hotspot
        Hugo_Symbol %in% genes_hotspot_vec & is_named_hotspot ~
          paste0(Hugo_Symbol, "_hotspot"),
        # Hotspot gene + coding mutation -> _other
        Hugo_Symbol %in% genes_hotspot_vec &
          Type %in% c("MUTATION", "TRUNC") ~
          paste0(Hugo_Symbol, "_other"),
        # Any gene + Synon + noncoding tracking -> _noncoding
        Hugo_Symbol %in% genes_noncoding_vec & Type == "Synon" ~
          paste0(Hugo_Symbol, "_noncoding"),
        # Hotspot gene + Synon without noncoding tracking -> skip
        Hugo_Symbol %in% genes_hotspot_vec & Type == "Synon" ~
          NA_character_,
        # Non-hotspot gene + Synon + explicit synonymous list -> plain col
        Hugo_Symbol %in% genes_synonymous_vec & Type == "Synon" ~
          Hugo_Symbol,
        # All remaining Synon (not in any tracking list) -> skip
        Type == "Synon" ~ NA_character_,
        # All other (non-Synon) mutations -> plain gene name column
        TRUE ~ Hugo_Symbol
      )
    ) |>
    dplyr::filter(!is.na(feature_col)) |>
    dplyr::group_by(Sample, feature_col) |>
    dplyr::summarise(
      score = max(score, na.rm = TRUE), .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols     = "Sample",
      names_from  = "feature_col",
      values_from = "score",
      values_fill = mutation_encoding$none
    ) |>
    dplyr::right_join(
      dplyr::tibble(Sample = all_samples), by = "Sample"
    ) |>
    dplyr::mutate(dplyr::across(
      -Sample,
      ~dplyr::coalesce(., as.numeric(mutation_encoding$none))
    ))

  if (verbose) message("  Mutation columns: ", ncol(mut_matrix) - 1)

  # ── 5. CNV matrix ────────────────────────────────────────────────────
  cnv_matrix <- NULL

  if (include_cnv) {
    if (verbose) message("Processing CNVs...")

    cnv_matrix <- cnv_data |>
      add_symbol() |>
      dplyr::filter(
        Type %in% c("GAIN", "AMP", "HETLOSS", "HOMDEL")
      ) |>
      dplyr::mutate(
        direction = dplyr::case_when(
          Type %in% c("GAIN", "AMP")       ~ "gain",
          Type %in% c("HETLOSS", "HOMDEL") ~ "loss"
        ),
        score = dplyr::case_when(
          Type == "GAIN"    ~ cnv_gain_encoding$gain,
          Type == "AMP"     ~ cnv_gain_encoding$amp,
          Type == "HETLOSS" ~ cnv_loss_encoding$hetloss,
          Type == "HOMDEL"  ~ cnv_loss_encoding$homdel
        ),
        col_name = paste0(Hugo_Symbol, "_", direction)
      ) |>
      dplyr::group_by(Sample, col_name) |>
      dplyr::summarise(
        score = max(score, na.rm = TRUE), .groups = "drop"
      ) |>
      tidyr::pivot_wider(
        id_cols     = "Sample",
        names_from  = "col_name",
        values_from = "score",
        values_fill = 0
      ) |>
      dplyr::right_join(
        dplyr::tibble(Sample = all_samples), by = "Sample"
      ) |>
      dplyr::mutate(dplyr::across(-Sample, ~dplyr::coalesce(., 0)))

    if (verbose) message("  CNV columns: ", ncol(cnv_matrix) - 1)
  }

  # ── 6. SV / translocation matrix ────────────────────────────────────
  sv_matrix <- NULL

  if (include_sv) {
    transloc_cols <- grep(
      "\\.transloc$", colnames(sample_annot), value = TRUE
    )

    if (length(transloc_cols) > 0) {
      if (verbose) {
        message(
          "Processing SVs (",
          length(transloc_cols),
          " translocation column(s))..."
        )
      }
      sv_col_rename <- setNames(
        sub("\\.transloc$", "_SV", transloc_cols),
        transloc_cols
      )

      sv_matrix <- sample_annot |>
        dplyr::select(
          Sample.ID, dplyr::all_of(transloc_cols)
        ) |>
        dplyr::rename(Sample = Sample.ID) |>
        dplyr::rename(
          dplyr::all_of(setNames(transloc_cols, sv_col_rename))
        ) |>
        dplyr::mutate(dplyr::across(
          -Sample,
          ~dplyr::if_else(
            !is.na(.) & . == 1L,
            as.numeric(sv_encoding$positive),
            as.numeric(sv_encoding$negative)
          )
        ))
    } else if (verbose) {
      message(
        "No '.transloc' columns found in sample annotation",
        " -- skipping SVs."
      )
    }
  }

  # ── 7. Combine all feature blocks ───────────────────────────────────
  if (verbose) message("Combining feature blocks...")

  result <- dplyr::rename(mut_matrix, sample_id = Sample)

  if (!is.null(cnv_matrix)) {
    result <- dplyr::left_join(
      result,
      dplyr::rename(cnv_matrix, sample_id = Sample),
      by = "sample_id"
    )
  }

  if (!is.null(sv_matrix)) {
    result <- dplyr::left_join(
      result,
      dplyr::rename(sv_matrix, sample_id = Sample),
      by = "sample_id"
    )
  }

  # ── 8. Optional gene subset ──────────────────────────────────────────
  if (!is.null(genes)) {
    if (verbose) {
      message(
        "Subsetting to ", length(genes), " requested gene(s)..."
      )
    }

    # Warn about outdated or missing symbols before subsetting
    check_lymphgen_genes(genes, result)

    # Match ALL column types for each requested gene using ^GENE(_|$)
    all_cols  <- setdiff(colnames(result), "sample_id")
    keep_cols <- unique(unlist(lapply(genes, function(g) {
      grep(paste0("^", g, "(_|$)"), all_cols, value = TRUE)
    })))

    result <- dplyr::select(
      result, sample_id, dplyr::all_of(keep_cols)
    )
  }

  if (verbose) {
    message(sprintf(
      "Done. Output: %d samples x %d features.",
      nrow(result), ncol(result) - 1
    ))
  }

  result
}

# Null-coalescing operator (avoids depending on rlang just for this)
`%||%` <- function(x, y) if (is.null(x)) y else x
