#' @title Tidy Lymphgen.
#'
#' @description Consolidate a column of LymphGen data in the original Subtype.Prediction output format to the GAMBLR tidy format.
#'
#' @details This function takes an incoming data frame (`df`) and consolidates a column of LymphGen data.
#' Specify the column with the lymphgen data to be processed with `lymphgen_column_in` and
#' what column to write the tidied data to with `lymphgen_column_out`.
#' In addition, the user can also run this function with `relevel = TRUE` (default is FALSE),
#' to return the output column as a factor with plot friendly levels.
#'
#' Composite calls (e.g. "EZB/ST2") are labelled "X-COMP" after one of their
#' classes. By default (`composite_rule = "precedence"`) that class follows a
#' fixed order, EZB > MCD > N1 > BN2 > ST2, so ST2 is the parent only when no
#' other class is present (e.g. "ST2/A53"). With `composite_rule = "confidence"`
#' the parent is the member class with the highest LymphGen confidence. A53 is
#' never the parent.
#'
#' For models that report confidences from more than one pass (GAMBL-U's
#' union of a copy-number and a no-copy-number pass: `Confidence.<class>` and
#' `NoCN.Confidence.<class>`), a class's confidence is the highest across its
#' columns, i.e. its strongest support from either pass. Ties go to the
#' precedence order. Composites with no confidence for any of their classes
#' keep the precedence label (with a message).
#'
#' @param df Input data frame.
#' @param lymphgen_column_in The name of the column with lymphgen data to be processed.
#' @param lymphgen_column_out The name of the column to write the tidied results (optional).
#' @param relevel If TRUE, will return the output column as a factor with plot-friendly levels.
#' @param composite_rule How the parent class of a composite is chosen:
#'   "precedence" (default; EZB > MCD > N1 > BN2 > ST2) or "confidence" (the
#'   member class with the highest confidence; see details).
#' @param confidence_columns With `composite_rule = "confidence"`: optional
#'   named list giving, for each class (BN2, EZB, MCD, N1, ST2), the column(s)
#'   of `df` holding its confidence. Default NULL: found from the column names
#'   `Confidence.<class>`, `NoCN.Confidence.<class>` (LymphGen output, GAMBL-U
#'   union), `CN_conf_<class>` and `noCN_conf_<class>`.
#'
#' @return A data frame with a tidied lymphGen column
#'
#' @import dplyr
#' @export
#'
#' @examples
#' metadata = GAMBLR.open::get_gambl_metadata()
#' lymphgen = tidy_lymphgen(df = metadata,
#'                          lymphgen_column_in = "lymphgen_with_cnv",
#'                          lymphgen_column_out = "lymphgen_with_cnv_tidy",
#'                          relevel = TRUE)
#'
#' # parent class of composites by confidence
#' calls = data.frame(Subtype.Prediction = c("EZB/ST2", "BN2/ST2", "MCD"),
#'                    Confidence.EZB = c(0.55, 0.02, 0.01),
#'                    Confidence.ST2 = c(0.91, 0.62, 0.00),
#'                    Confidence.BN2 = c(0.01, 0.95, 0.03),
#'                    Confidence.MCD = c(0.00, 0.01, 0.98))
#' tidy_lymphgen(calls, composite_rule = "confidence")
#'
tidy_lymphgen = function(df,
                         lymphgen_column_in = "Subtype.Prediction",
                         lymphgen_column_out = "Subtype.Prediction",
                         relevel = FALSE,
                         composite_rule = "precedence",
                         confidence_columns = NULL){

  composite_rule = match.arg(composite_rule, c("precedence", "confidence"))
  calls_in = as.character(df[[lymphgen_column_in]])
  df = mutate(df, {{ lymphgen_column_out }} := case_when(
    !grepl("/", .data[[lymphgen_column_in]])~.data[[lymphgen_column_in]],
    grepl("EZB", .data[[lymphgen_column_in]])~"EZB-COMP",
    grepl("MCD", .data[[lymphgen_column_in]])~"MCD-COMP",
    grepl("N1", .data[[lymphgen_column_in]])~"N1-COMP",
    grepl("BN2", .data[[lymphgen_column_in]])~"BN2-COMP",
    grepl("ST2", .data[[lymphgen_column_in]])~"ST2-COMP"
  ))
  if(composite_rule == "confidence"){
    classes = c("EZB", "MCD", "N1", "BN2", "ST2")   # precedence order, for ties
    if(is.null(confidence_columns)){
      confidence_columns = lapply(classes, function(k){
        intersect(c(paste0("Confidence.", k), paste0("NoCN.Confidence.", k), paste0("CN_conf_", k), paste0("noCN_conf_", k)), colnames(df))
      })
      names(confidence_columns) = classes
    }
    missing_cols = setdiff(unlist(confidence_columns), colnames(df))
    if(length(missing_cols) > 0){
      stop("confidence_columns not found in df: ", paste(missing_cols, collapse = ", "))
    }
    if(all(lengths(confidence_columns[intersect(classes, names(confidence_columns))]) == 0)){
      stop("composite_rule = \"confidence\" needs per-class confidence columns (e.g. Confidence.EZB); none found")
    }
    # each class's confidence: the highest across its columns (e.g. both passes of a union model)
    conf = sapply(classes, function(k){
      cols = confidence_columns[[k]]
      if(length(cols) == 0) return(rep(NA_real_, nrow(df)))
      m = as.matrix(as.data.frame(df)[, cols, drop = FALSE]); storage.mode(m) = "numeric"
      v = suppressWarnings(apply(m, 1, max, na.rm = TRUE)); v[!is.finite(v)] = NA_real_; v
    })
    if(is.null(dim(conf))) conf = matrix(conf, nrow = nrow(df), dimnames = list(NULL, classes))
    is_comp = which(grepl("/", calls_in))
    n_fallback = 0
    out = as.character(df[[lymphgen_column_out]])
    for(i in is_comp){
      members = intersect(classes, strsplit(calls_in[i], "/")[[1]])
      if(length(members) == 0) next
      cv = conf[i, members]
      if(all(is.na(cv))){ n_fallback = n_fallback + 1; next }
      out[i] = paste0(members[which.max(replace(cv, is.na(cv), -Inf))], "-COMP")   # which.max keeps the first maximum: precedence order breaks ties
    }
    df[[lymphgen_column_out]] = out
    if(n_fallback > 0){
      message(n_fallback, " composite call(s) without confidences for their classes keep the precedence label")
    }
  }
  if(relevel){
    df <- df %>%
      mutate({
        {
          lymphgen_column_out
        }
      } := factor(
        .data[[lymphgen_column_out]],
        levels = c(
          "EZB",
          "EZB-COMP",
          "ST2",
          "ST2-COMP",
          "BN2",
          "BN2-COMP",
          "N1",
          "N1-COMP",
          "MCD",
          "MCD-COMP",
          "A53",
          "Other"
        )
      ))
  }
  return(df)
}
