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
#' @param df Input data frame.
#' @param lymphgen_column_in The name of the column with lymphgen data to be processed.
#' @param lymphgen_column_out The name of the column to write the tidied results (optional).
#' @param relevel If TRUE, will return the output column as a factor with plot-friendly levels.
#'
#' @return A data frame with a tidied lymphGen column
#'
#' @import dplyr purrr readr
#' @export
#'
#' @examples
#' metadata = GAMBLR.data::gambl_metadata
#' lymphgen = tidy_lymphgen(df = metadata,
#'                          lymphgen_column_in = "lymphgen_with_cnv",
#'                          lymphgen_column_out = "lymphgen_with_cnv_tidy",
#'                          relevel = TRUE)
#'
tidy_lymphgen = function(df,
                         lymphgen_column_in = "Subtype.Prediction",
                         lymphgen_column_out = "Subtype.Prediction",
                         relevel = FALSE){

  df = mutate(df, {{ lymphgen_column_out }} := case_when(
    !grepl("/", .data[[lymphgen_column_in]])~.data[[lymphgen_column_in]],
    grepl("EZB", .data[[lymphgen_column_in]])~"EZB-COMP",
    grepl("MCD", .data[[lymphgen_column_in]])~"MCD-COMP",
    grepl("N1", .data[[lymphgen_column_in]])~"N1-COMP",
    grepl("BN2", .data[[lymphgen_column_in]])~"BN2-COMP",
    grepl("ST2", .data[[lymphgen_column_in]])~"ST2-COMP"
  ))
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
