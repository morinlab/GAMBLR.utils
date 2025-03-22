#' Map region to bins
#'
#' @param query_region
#' @param regions
#' @param first
#'
#' @return a named list
#' @export
#'
#' @examples
#'
#' gene_region <- gene_to_region("TP53")
#'
#' all_bins <- colnames(cn_state_matrix)
#'
#' TP53_bin <- map_region_to_bins(gene_region, all_bins, first = TRUE)
#'
map_regions_to_bins <- function(query_regions,
                               regions,
                               query_type = "region",
                               first = TRUE) {
    if(query_type == "gene"){
        
        query_regions <- gene_to_region(query_regions)
        #print(query_regions)
    }else{
        names(query_regions) <- query_regions
    }
    region_matches = list()
    for(query_region in names(query_regions)) { 

      these_coords <- suppressMessages(region_to_chunks(query_regions[query_region]))
      these_coords$chromosome <- gsub("chr", "", these_coords$chromosome)
      all_matches <- c()
      for (r in regions) {
        region_coords <- region_to_chunks(r)
        region_coords$chromosome <- gsub("chr", "", region_coords$chromosome)
        if (these_coords$chromosome == region_coords$chromosome) {
            if (((as.integer(these_coords$start) > as.integer(region_coords$start)) &
                (as.integer(these_coords$start) < as.integer(region_coords$end))) ||
                (as.integer(these_coords$start) < as.integer(region_coords$start) &
                    as.integer(these_coords$end) > as.integer(region_coords$end))) {
                if (first) {
                    # just return the first match
                    #print("match")
                    #print(query_region)
                    #print(r)
                    region_matches[[query_region]] <- r
                    next;
                } else {
                    all_matches <- c(all_matches, r)
                }
            }
        }
      }
      if(length(all_matches)){
        region_matches[[query_region]] <- all_matches
      }
    }
    return(region_matches)
}
