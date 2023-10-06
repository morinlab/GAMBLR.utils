#' @title Gene To Region.
#'
#' @description Return coordinates for a given gene or a set of genes.
#'
#' @details This function takes one or multiple gene names, either as hugo symbols or Ensembl IDs
#' and returns the coordinates for the selected genes in multiple formats (`return_as`).
#' The possible return formats are; bed, data frame and in "region" format (chr:start-end).
#' For returning genes residing in specific regions, see [GAMBLR.utils::region_to_gene].
#'
#' @param gene_symbol A vector of one or more gene symbols.
#' @param ensembl_id A vector of one or more Ensembl IDs.
#' @param genome_build Reference genome build. Possible values are "grch37" (default) or "hg38".
#' @param return_as Specify the type of return. Default is "region" (chr:start-end), other acceptable arguments are "bed" and "df".
#' @param sort_regions A boolean parameter (TRUE is the default) indicating whether regions should be sorted by chomosome and start location.
#'
#' @return A data frame, or a string with region(s) for the provided gene(s).
#'
#' @import dplyr
#' @export
#'
#' @examples
#' bcl2_region = gene_to_region(gene_symbol = "BCL2",
#'                              genome_build = "grch37")
#'
#' bcl2_region = gene_to_region(ensembl_id = "ENSG00000171791",
#'                              genome_build = "grch37")
#'
gene_to_region = function(gene_symbol,
                          ensembl_id,
                          genome_build = "grch37",
                          return_as = "region",
                          sort_regions = TRUE){
  
  stopifnot('`genome_build` parameter must be "grch37" or "hg38"' = genome_build %in% c("grch37", "hg38"))
  stopifnot('`return_as` parameter must be "region", "bed" or "df"' = return_as %in% c("region", "bed", "df"))
  stopifnot('One and only one of the `gene_symbol` and `ensembl_id` parameters must be given to this function' = sum(missing(gene_symbol), missing(ensembl_id)) == 1)
  
  #set mart based on selected genome projection
  if(genome_build == "grch37"){
    gene_coordinates = GAMBLR.data::grch37_gene_coordinates
    chr_select = paste0(c(c(1:22),"X","Y"))
  }else{
    gene_coordinates = GAMBLR.data::hg38_gene_coordinates
    chr_select = paste0("chr", c(c(1:22),"X","Y"))
  }
  
  #filter on gene_symbol/ensembl_id
  if(!missing(gene_symbol) && missing(ensembl_id)){
    gene_coordinates = dplyr::filter(gene_coordinates, hugo_symbol %in% gene_symbol)
    genes_not_vailable = gene_symbol[! gene_symbol %in% gene_coordinates$hugo_symbol]
  }
  if(missing(gene_symbol) && !missing(ensembl_id)){
    gene_coordinates = dplyr::filter(gene_coordinates, ensembl_gene_id %in% ensembl_id)
    genes_not_vailable = ensembl_id[! ensembl_id %in% gene_coordinates$ensembl_gene_id]
  }
  
  #print list of genes that have no region info available
  if(length(genes_not_vailable) > 0){
    paste(genes_not_vailable, collapse = ", ") %>%
      gettextf("Some input gene(s) have no region info available. They are:\n%s.", .) %>%
      message
  }
  
  region = dplyr::select(gene_coordinates, chromosome, start, end, gene_name, hugo_symbol, ensembl_gene_id) %>%
    as.data.frame() %>% 
    dplyr::filter(chromosome %in% chr_select)
  
  if(sort_regions){
    if(genome_build == "grch37"){
      chrm_num = region$chromosome
    }else{
      chrm_num = sub("chr", "", region$chromosome)
    }
    chrm_num = factor(chrm_num, levels = c(1:22, "X", "Y"), ordered = TRUE)
    region = dplyr::arrange(region, chrm_num, start)
  }else{
    #make the output gene order the same as the input
    if(!missing(gene_symbol) && missing(ensembl_id)){
      region = arrange(region, match(hugo_symbol, gene_symbol))
    }
    if(missing(gene_symbol) && !missing(ensembl_id)){
      region = arrange(region, match(hugo_symbol, ensembl_id))
    }
  }
  
  region[region == ""] = NA
  region = distinct(region, .keep_all = TRUE)
  
  if(return_as == "bed"){
    #return one-row data frame with first 4 standard BED columns. TODO: Ideally also include strand if we have access to it in the initial data frame
    region = dplyr::select(region, chromosome, start, end, hugo_symbol)
    
  }else if(return_as == "df"){
    region = region
    
  }else{
    #default: return in chr:start-end format
    if(!missing(gene_symbol) && missing(ensembl_id)){
      ids = dplyr::pull(region, hugo_symbol)
    }
    if(missing(gene_symbol) && !missing(ensembl_id)){
      ids = dplyr::pull(region, ensembl_gene_id)
    }
    region = setNames(
      paste0(region$chromosome, ":", region$start, "-", region$end, recycle0 = TRUE),
      ids
    )
  }
  
  if(return_as %in% c("bed", "df")){
    if(!missing(gene_symbol)){
      message(paste0(nrow(region[!is.na(region$chromosome),]), " region(s) returned for ", length(gene_symbol), " gene(s)"))
    }
    
    if(!missing(ensembl_id)){
      message(paste0(nrow(region[!is.na(region$chromosome),]), " region(s) returned for ", length(ensembl_id), " gene(s)"))
    }
  }else{
    if(!missing(gene_symbol)){
      message(paste0(length(region[!is.na(region)]), " region(s) returned for ", length(gene_symbol), " gene(s)"))
    }
    
    if(!missing(ensembl_id)){
      message(paste0(length(region[!is.na(region)]), " region(s) returned for ", length(ensembl_id), " gene(s)"))
    }
  }
  return(region)
}
