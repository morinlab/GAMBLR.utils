#' Map regions to bins
#'
#' Given query gene or regions, return the overlapping bin(s)
#' from a candidate set of regions. Overlap is assessed by chromosome match and
#' whether the query region starts inside a candidate region or fully contains it.
#'
#' @param query_regions Character vector of genomic regions (e.g. "chr17:777-999")
#'   or gene symbols if `query_type = "gene"`.
#' @param regions Character vector of candidate regions (bins) to match against.
#' @param query_type One of `c("region", "gene")`. If `"gene"`, `query_regions`
#'   is interpreted as gene symbols and converted to regions via
#'   `gene_to_region()`.
#' @param first Logical; if `TRUE`, return only the first match per query region.
#'   If `FALSE`, return all matches.
#'
#' @return A named list keyed by query region name with value(s) of matched
#'   region(s). When `first = TRUE`, each element is a single region string; when
#'   `first = FALSE`, each element is a character vector of all matches.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_region <- gene_to_region("TP53")
#' 
#' all_bins <- colnames(cn_state_matrix)
#' 
#' TP53_bin <- map_regions_to_bins(gene_region, all_bins, first = TRUE)
#'
#' # multiple regions
#' my_regions <- c("chr17:7500000-7600000", "chr1:100000-120000")
#' region_bins <- map_regions_to_bins(my_regions, all_bins, first = FALSE)
#'}
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
