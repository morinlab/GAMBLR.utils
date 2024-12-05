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
#' @import dplyr
#' @export
#'
#' @examples
#' #get a region, used for input in Example 1
#' genes_of_interest = gene_to_region(gene_symbol = c("MYC", "BCL2"),
#'                                    projection = "grch37",
#'                                    return_as = "df")
#' 
#' #Example 1 - Use a data frame with 1 region:
#' in_my_region = region_to_gene(region = genes_of_interest)
#'                         
#' #Example 2 - Use a region specified as a vector of characters:
#' in_my_region = region_to_gene(region = "chr8:127735433-127742951", 
#'                               projection = "hg38", 
#'                               gene_format = "ensembl", 
#'                               verbose = TRUE)
#'                               
#' #Example 3 - Use a region with no genes in it, and return the empty data frame anyway:
#' in_my_region = region_to_gene(region = "chr14:105862865-105863058",
#'                               projection = "hg38", 
#'                               verbose = TRUE, 
#'                               return_empty = TRUE)
#'
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

  if(!is.data.frame(region)){
    region_table = as.data.frame(region)
  }else if(is.character(region)){
    split_chunks = unlist(strsplit(region, ":"))
    split_chunks = unlist(strsplit(split_chunks, "-"))
    chromosome = split_chunks[1]
    start = split_chunks[2]
    end = split_chunks[3]
    region = cbind(chromosome, start, end) %>%
      as.data.frame()

    region$chromosome = as.character(region$chromosome)
    region$start = as.double(region$start)
    region$end = as.double(region$end)
  }

  #intersect regions
  intersect = cool_overlaps(
    region,
    gene_list,
    columns1 = c("chromosome", "start", "end"),
    columns2 = c("chromosome", "start", "end")
  )
  colnames(intersect)[2] = "region_start"
  colnames(intersect)[3] = "region_end"
  colnames(intersect) <- gsub("\\.y", "", colnames(intersect))
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

  #paste chr in chromosome column, if not there
  if(!grepl("chr", genes$chromosome[1])){
    genes = mutate(genes, chromosome = paste0("chr", chromosome))}

  if(verbose){
    message(paste0(nrow(genes), " gene(s) returned for ", nrow(region), " region(s)")) 
  }

  return(genes)
}
