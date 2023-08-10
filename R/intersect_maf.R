#' @title Intersect MAF.
#'
#' @description Perform set operations on two MAFs.
#'
#' @details Perform set operations on two MAFs.
#'
#' @param maf1 First list of MAFs.
#' @param maf2 Second list of MAFs.
#' @param set_returned List of MAFs that doesn't share the same start positions as the other list of MAFs. Accepted commands are; "maf1_only" and "maf2_only", default is "maf1_only".
#'
#' @return Set of MAFs with start positions that don't match the start positions in the other supplied MAF file.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' intersected_mafs_l1 = intersect_maf(maf_list1, maf_list2, "maf1_only")
#' intersected_mafs_l2 = intersect_maf(maf_list1, maf_list2, "maf2_only")
#' }
#'
intersect_maf = function(maf1,
                         maf2,
                         set_returned = "maf1_only"){

  if(set_returned == "maf1_only"){
    maf_set = dplyr::filter(maf1, !Start_Position %in% maf2$Start_Position)
  }else if(set_returned == "maf2_only"){
    maf_set = dplyr::filter(maf2, !Start_Position %in% maf1$Start_Position)
  }
  return(maf_set)
}
