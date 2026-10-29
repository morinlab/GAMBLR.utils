#' @title Region To Gene.
#'
#' @description Return genes residing in defined region(s).
#'
#' @details This function takes a region as a vector of characters, or a data frame with of regions (e.g output from `gene_to_region(return_as = "df")`).
#' and returns the genes residing withing the specified region. For the other way around (i.e gene to regions, refer to [GAMBLR.utils::gene_to_region]).
#' If the user provides a region as a vector of characters there is the option to set `return_empty = TRUE`. 
#' This forces the function to return an empty data frame when there are no genes found within the specified region for the selected projection.
#' Note, if the user provides a data frame with regions, any regions that no genes are found for will be excluded automatically in the the returned object.
#'
#' @param region Regions to intersect genes with, this can be either a data frame with regions sorted in the following columns; chromosome, start, end. Or it can be a character vector in "region" format, i.e chromosome:start-end.
#' @param gene_format Parameter for specifying the format of returned genes, default is "hugo", other acceptable inputs are "ensembl".
#' @param projection Reference genome build. Default is grch37.
#' @param verbose Set to TRUE for noisier output to the console. Useful or debugging purposes. Default is FALSE.
#' @param return_empty Set to TRUE to force returning an empty data frame when no genes are found within the specified region for the selected projection. 
#' This only applies when the function is given a region as a vector of characters. See function documentation for more information. Default is FALSE.
#'
#' @return A data frame with columns detailing the gene/region.
#'
#' @import dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' #get a region, used for input in Example 1
#' regions_of_interest = gene_to_region(gene_symbol = c("MYC", "BCL2"),
#'                                    projection = "grch37",
#'                                    return_as = "df")
#' 
#' regions_of_interest
#' 
#' 
#'                         
#' \dontrun{
#' #Example 2 - Use a region specified as a vector of characters:
#' in_my_region = region_to_gene(region = "chr8:127735433-127742951", 
#'                               projection = "hg38", 
#'                               gene_format = "ensembl", 
#'                               verbose = TRUE)
#' print(in_my_region)
#'                               
#' #Example 3 - Use a region with no genes in it, and return the empty data frame anyway:
#' in_my_region = region_to_gene(region = "chr14:105862865-105863058",
#'                               projection = "hg38", 
#'                               verbose = TRUE, 
#'                               return_empty = TRUE)
#'}

region_to_gene = function(region,
                          gene_format = "hugo",
                          projection = "grch37",
                          verbose = TRUE,
                          return_empty = FALSE){

  #set mart based on selected genome projection
  if(projection == "grch37"){
    gene_list = GAMBLR.data::grch37_gene_coordinates
  }else if(projection == "hg38"){
    gene_list = GAMBLR.data::hg38_gene_coordinates
  }

  #rename columns to match downstream formats
  colnames(gene_list)[1] = "ensembl_gene_id"
  colnames(gene_list)[2] = "chromosome"
  colnames(gene_list)[3] = "start"
  colnames(gene_list)[4] = "end"
  colnames(gene_list)[5] = "gene_name"
  colnames(gene_list)[6] = "hugo_symbol"

  gene_list = as.data.frame(gene_list)

  if(is.character(region)){
    # one or more regions in "chromosome:start-end" format
    region = gsub(",", "", region)
    parsed = regmatches(region, regexec("^([^:]+):([0-9]+)-([0-9]+)$", region))
    malformed = lengths(parsed) == 0
    if(any(malformed)){
      stop(paste("regions must be in chromosome:start-end format:",
                 paste(region[malformed], collapse = ", ")))
    }
    region = data.frame(
      chromosome = vapply(parsed, `[`, character(1), 2),
      start = as.double(vapply(parsed, `[`, character(1), 3)),
      end = as.double(vapply(parsed, `[`, character(1), 4))
    )
  }else if(is.data.frame(region)){
    region = as.data.frame(region)
  }else{
    stop("region must be a character vector (chromosome:start-end) or a data frame with chromosome, start and end columns")
  }
  # match the chr-prefix convention of the gene coordinates
  region$chromosome = as.character(region$chromosome)
  if(grepl("^chr", gene_list$chromosome[1])){
    region$chromosome = ifelse(grepl("^chr", region$chromosome),
                               region$chromosome,
                               paste0("chr", region$chromosome))
  }else{
    region$chromosome = gsub("^chr", "", region$chromosome)
  }

  #intersect regions: any overlap, including genes that cross a region boundary
  region_dt = data.table::as.data.table(region)[, c("chromosome", "start", "end")]
  region_dt[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  gene_dt = data.table::as.data.table(gene_list)
  gene_dt[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  data.table::setkey(gene_dt, chromosome, start, end)
  intersect = data.table::foverlaps(region_dt, gene_dt, type = "any", nomatch = NULL)
  data.table::setnames(intersect, c("i.start", "i.end"), c("region_start", "region_end"))
  #transform object to data frame
  inter_df = as.data.frame(intersect)

  #organize columns to match the expected format
  if(gene_format == "hugo"){
    genes = dplyr::select(inter_df, chromosome, start, end, hugo_symbol, region_start, region_end)
  }else if(gene_format == "ensembl"){
    genes = dplyr::select(inter_df, chromosome, start, end, ensembl_gene_id, region_start, region_end)
  }
  
  genes = as.data.frame(genes) %>%
    dplyr::arrange(chromosome, start) %>%
    distinct(.keep_all = TRUE)

  #add checks for 0 genes in the region
  if(nrow(genes) == 0){
    if(return_empty){
      if(verbose){
        message("No genes found within the specified region, an empty data frame will be returned")
      }
      return(genes) 
    }else{
      stop(paste0("The specified region does not overlap with any genes in the selected projection.", "\n", 
                "  If you want the empty data frame back, run this function with `return_empty = TRUE`"))
    }
  }

  # chromosome names follow the build's convention, as in gene_to_region():
  # no prefix for grch37 ("19"), "chr" prefix for hg38 ("chr19")

  if(verbose){
    message(paste0(nrow(genes), " gene(s) returned for ", nrow(region), " region(s)")) 
  }

  return(genes)
}
