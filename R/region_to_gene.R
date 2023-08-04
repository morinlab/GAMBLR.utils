#' @title Region To Gene.
#'
#' @description Return genes residing in defined region(s).
#'
#' @details This function takes a region as a vector of characters, or a data frame with of regions (e.g output from `gene_to_region(return_as="df")`).
#' and returns the genes residing withing the specified region. For the other way around (i.e gene to regions, refer to [GAMBLR::gene_to_region]).
#'
#' @param region Regions to intersect genes with, this can be either a data frame with regions sorted in the following columns; chromosome, start, end. Or it can be a character vector in "region" format, i.e chromosome:start-end.
#' @param gene_format Parameter for specifying the format of returned genes, default is "hugo", other acceptable inputs are "ensembl".
#' @param genome_build Reference genome build.
#'
#' @return A data frame with gene(s) that are residing in the specified region(s).
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr
#' @export
#'
#' @examples
#' myc_region = gene_to_region(gene_symbol = "MYC",
#'                             genome_build = "grch37",
#'                             return_as = "df")
#'
#' region = region_to_gene(region = myc_region,
#'                         gene_format = "hugo",
#'                         genome_build = "grch37")
#'
region_to_gene = function(region,
                          gene_format = "hugo",
                          genome_build = "grch37"){

  #set mart based on selected genome projection
  if(genome_build == "grch37"){
    gene_list = GAMBLR.data::grch37_gene_coordinates
  }else if(genome_build == "hg38"){
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

  if(is.data.frame(region)){
    region_table = as.data.table(region)
  }else if(is.character(region)){
    split_chunks = unlist(strsplit(region, ":"))
    split_chunks = unlist(strsplit(split_chunks, "-"))
    chromosome = split_chunks[1]
    start = split_chunks[2]
    end = split_chunks[3]
    region = cbind(chromosome, start, end) %>%
      as.data.frame()

    region_table = as.data.table(region)

    region_table$chromosome = as.character(region_table$chromosome)
    region_table$start = as.double(region_table$start)
    region_table$end = as.double(region_table$end)
  }

  #transform regions to data tables
  gene_table = as.data.table(gene_list)

  #set keys
  data.table::setkey(region_table, chromosome, start, end)
  data.table::setkey(gene_table, chromosome, start, end)

  #intersect regions
  intersect = data.table::foverlaps(region_table, gene_table, nomatch = 0)
  colnames(intersect)[7] = "region_start"
  colnames(intersect)[8] = "region_end"

  #transform object to data frame
  inter_df = as.data.frame(intersect)

  #organize columns to match the expected format
  if(gene_format == "hugo"){
    genes = dplyr::select(inter_df, chromosome, start, end, hugo_symbol, region_start, region_end)
  }else if(gene_format == "ensembl"){
    genes = dplyr::select(inter_df, chromosome, start, end, ensembl_gene_id, region_start, region_end)
  }

  #paste chr in chromosome column, if not there
  if(!str_detect(genes$chromosome[1], "chr")){
    genes = mutate(genes, chromosome = paste0("chr", chromosome))}

  genes = as.data.frame(genes) %>%
    dplyr::arrange(chromosome, start) %>%
    distinct(.keep_all = TRUE)

  message(paste0(nrow(genes), " gene(s) returned for ", nrow(region), " region(s)"))

  return(genes)
}
