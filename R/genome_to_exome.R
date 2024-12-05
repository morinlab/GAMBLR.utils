#' @title Genome To Exome.
#'
#' @description Subset maf file to only features that would be available in the WEX data.
#'
#' @details To subset an incoming MAF data frame to only show features that would be available in WEX data this function was developed.
#' Pass the incoming MAF (genome) to the `maf` parameter as the only required parameter to run this function. Other parameters such as `custom_bed`,
#' `genome_build`, `padding`, and `chr_prefixed` are also available for greater control of how this function operates.
#' Refer to parameter descriptions for more information on how to use the available parameters.
#'
#' @param maf Incoming maf object. Must be maf-like data frame. Required parameter. Minimum columns that should be present are Chromosome, Start_Position, and End_Position.
#' @param custom_bed Optional argument specifying a path to custom bed file for covered regions. Must be bed-like and contain chrom, start, and end position information in the first 3 columns. Other columns are disregarded if provided.
#' @param genome_build String indicating genome build of the maf file. Default is grch37, but can accept modifications of both grch37- and hg38-based builds.
#' @param padding Numeric value that will be used to pad probes in WEX data from both ends. Default is 100. After padding, overlapping features are squished together.
#' @param chr_prefixed Is the data chr-prefixed or not? Default is FALSE.
#'
#' @return A data frame of a maf-like object with the same columns as in input, but where rows are only kept for features that would be present as if the sample is WEX.
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#'
#' #get all ssm in the MYC aSHM region
#' myc_ashm_maf = get_ssm_by_region(region = "8:128748352-128749427")
#'
#' #get mutations with 100 bp padding (default)
#' maf = genome_to_exome(maf = myc_ashm_maf)
#'
#' #get mutations covered in WEX with no padding
#' maf = genome_to_exome(maf = myc_ashm_maf,
#'                 padding = 0)
#'
genome_to_exome = function(maf,
                           custom_bed,
                           genome_build = "grch37",
                           padding = 100,
                           chr_prefixed = FALSE){

  if(missing(custom_bed)){
      # first check that the genome build provided is supported
      if(! genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37", "hg38", "GRCh38", "grch38")){
        stop("The genome build specified is not currently supported. Please refer to genome build in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38.")
      }else if(genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
        this_genome_coordinates = GAMBLR.data::target_regions_grch37 # if the genome build is a flavour of hg19, get its exome space
      }else if(genome_build %in% c("hg38", "GRCh38", "grch38")){
        this_genome_coordinates = GAMBLR.data::target_regions_hg38 # exome space for the variations of hg38
      }
  }else{
      this_genome_coordinates = read_tsv(custom_bed)
      this_genome_coordinates = this_genome_coordinates[,1:3]
      colnames(this_genome_coordinates)[1:3] = c("chrom", "start", "end")
  }

  # pad the ends of the probes with the padding length
  this_genome_coordinates = this_genome_coordinates %>%
    dplyr::mutate(chrom = as.character(chrom), start = start - padding, end = end + padding)

  # collapse regions if the padding results in 2 features that are overlapping
  this_genome_coordinates = this_genome_coordinates %>%
    dplyr::arrange(chrom, start) %>%
    group_by(chrom) %>%
    dplyr::mutate(indx = c(0, cumsum(as.numeric(lead(start)) > cummax(as.numeric(end)))[-n()])) %>%
    group_by(chrom, indx) %>%
    dplyr::summarise(start = first(start), end = last(end)) %>%
    dplyr::select(-indx) %>%
    ungroup %>%
    dplyr::mutate_if(is.factor, as.character)

  # handle the chr prefixes
  if(chr_prefixed & ! grepl("chr", this_genome_coordinates$chrom[1])){
    this_genome_coordinates = this_genome_coordinates %>%
      dplyr::mutate(chrom = paste0("chr", chrom))
  }else if(!chr_prefixed & grepl("chr", this_genome_coordinates$chrom[1])){
    this_genome_coordinates = this_genome_coordinates %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom, ignore.case = TRUE))
  }

  # subset to only features covered in exome with provided padding
  features_in_exome = cool_overlaps(
    maf,
    this_genome_coordinates,
    columns2 = c("chrom", "start", "end")
    ) %>%
    as.data.frame() %>%
    dplyr::select(colnames(maf)) # make sure columns and their order is consitent with the input maf
  return(features_in_exome)
}
