#' @title Annotate SVs.
#'
#' @description Annotate a data frame of SV breakpoints represented in an extended BEDPE format.
#'
#' @details Specify a data frame with SVs (preferably the output from [GAMBLR::get_manta_sv]) to the `sv_df` parameter and get back the same data frame with SV annotations.
#'
#' @param sv_data A data frame of SVs. This should be the output of get_manta_sv. If you aren't using the database backend you can supply your own data frame in the format show below.
#' Most of this data is directly from the bedpe files that are obtained by converting the Manta outputs from VCF.
#' Only the columns for both chromosomes, coordinates and strand plus SOMATIC_SCORE and tumour_sample_id are absolutely required.
#'  CHROM_A  START_A    END_A CHROM_B  START_B    END_B NAME SOMATIC_SCORE STRAND_A STRAND_B TYPE FILTER VAF_tumour VAF_normal DP_tumour DP_normal tumour_sample_id normal_sample_id pair_status
#'   1  1556541  1556547       1  1556664  1556670    .            40        -        -  BND   PASS      0.145          0        55        73  00-14595_tumorA  00-14595_normal     matched
#' @param partner_bed Optional bed-format data frame to use for annotating oncogene partners (e.g. enhancers). required columns are: chrom,start,end,gene
#' @param with_chr_prefix Optionally request that chromosome names are returned with a chr prefix. Default is FALSE.
#' @param collapse_redundant Remove reciprocal events and only return one per event. Default is FALSE.
#' @param return_as Stated format for returned output, default is "bedpe". Other accepted output formats are "bed" and "bedpe_entrez" (to keep entrez_ids for compatibility with portal.R and cBioPortal).
#' @param blacklist A vector of regions to be removed from annotations. Default coordinates are in respect to hg19.
#' @param genome_build Reference genome build parameter, default is grch37.
#'
#' @return A data frame with annotated SVs (gene symbol and entrez ID).
#' 
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import tidyr dplyr stringr
#' @export
#'
#' @examples
#' sv_df = get_manta_sv(verbose = FALSE)
#' annotated_entrez = annotate_sv(sv_data = sv_df,
#'                                with_chr_prefix = FALSE,
#'                                collapse_redundant = FALSE,
#'                                return_as = "bedpe_entrez", 
#'                                genome_build = "grch37")
#'
annotate_sv = function(sv_data,
                       partner_bed,
                       with_chr_prefix = FALSE,
                       collapse_redundant = FALSE,
                       return_as = "bedpe",
                       blacklist = c(60565248, 30303126, 187728894, 101357565, 101359747, 161734970, 69400840, 65217851, 187728889, 187728888,187728892, 187728893,188305164),
                       genome_build = "grch37"){
  
  # remove duplicate rows in sv_data, if any
  sv_data <- filter(sv_data, !duplicated(sv_data))
  
  bedpe1 = sv_data %>%
    dplyr::select("CHROM_A", "START_A", "END_A", "tumour_sample_id", ends_with("SCORE"), "STRAND_A")
  
  bedpe2 = sv_data %>%
    dplyr::select("CHROM_B", "START_B", "END_B", "tumour_sample_id", ends_with("SCORE"), "STRAND_B")
  
  colnames(bedpe1) = c("chrom", "start", "end", "tumour_sample_id", "score", "strand1")
  colnames(bedpe2) = c("chrom", "start", "end", "tumour_sample_id", "score", "strand2")
  suppressWarnings({
    if(any(grepl("chr", bedpe1$chrom))){
      bedpe1 = dplyr::mutate(bedpe1, chrom = str_replace(chrom, "chr", ""))
      bedpe2 = dplyr::mutate(bedpe2, chrom = str_replace(chrom, "chr", ""))
    }
  })
  if(missing(partner_bed)){
    if(genome_build == "hg38"){
      ig_regions = GAMBLR.data::hg38_partners
    }else{
      ig_regions = GAMBLR.data::grch37_partners
    }
  }else{
    ig_regions = partner_bed
    if(!"entrez" %in% colnames(ig_regions)){
      ig_regions$entrez = 0
    }
  }
  if(genome_build == "hg38"){
    oncogene_regions = GAMBLR.data::hg38_oncogene
  }else{
    oncogene_regions = GAMBLR.data::grch37_oncogene
  }
  y = data.table::as.data.table(oncogene_regions)
  data.table::setkey(y, chrom, start, end)
  
  #use foverlaps to get oncogene SVs
  a = data.table::as.data.table(bedpe1)
  a.onco = data.table::foverlaps(a, y, type = "any", mult = "first") #oncogene-annotated bedpe for the first breakpoints
  
  b = data.table::as.data.table(bedpe2)
  b.onco = data.table::foverlaps(b, y, type = "any", mult = "first") #oncogene-annotated bedpe for the first breakpoints
  
  #insist oncogene breakpoints are anchored in an IG or superenhancer region (currently just IG or BCL6)
  #other end of breakpoint
  a.onco.break = a.onco[which(!is.na(a.onco$start)), c("chrom", "i.start", "i.end", "tumour_sample_id", "gene", "entrez", "score", "strand1")]
  b.onco.break = b.onco[which(!is.na(b.onco$start)), c("chrom", "i.start", "i.end", "tumour_sample_id", "gene", "entrez", "score", "strand2")]
  
  a.partner = b[which(!is.na(a.onco$start)),]
  b.partner = a[which(!is.na(b.onco$start)),]
  
  y = data.table::as.data.table(ig_regions)
  data.table::setkey(y, chrom, start, end)
  
  a.ig = data.table::foverlaps(a.partner, y, type = "any", mult = "first")
  b.ig = data.table::foverlaps(b.partner, y, type = "any", mult = "first")
  
  a.ig = a.ig[,c("chrom", "i.start", "i.end", "strand2", "gene")]
  b.ig = b.ig[,c("chrom", "i.start", "i.end", "strand1", "gene")]
  
  a.annotated.both = cbind(a.onco.break, a.ig)
  colnames(a.annotated.both) = c("chrom1", "start1", "end1", "tumour_sample_id", "gene", "entrez", "score", "strand1", "chrom2", "start2", "end2", "strand2", "partner")
  
  b.annotated.both = cbind(b.onco.break, b.ig)
  colnames(b.annotated.both) = c("chrom2", "start2", "end2", "tumour_sample_id", "gene", "entrez", "score", "strand2", "chrom1", "start1", "end1", "strand1", "partner")
  
  all.annotated = rbind(a.annotated.both, b.annotated.both)
  
  all.annotated$fusion = dplyr::pull(tidyr::unite(all.annotated, fusion, partner, gene, sep = "-"), fusion)
  all.annotated = dplyr::filter(all.annotated, !fusion %in% c("BCL6-BCL6", "CIITA-CIITA", "FOXP1-FOXP1"))
  
  all.annotated = dplyr::filter(all.annotated, !start1 %in% blacklist)
  all.annotated = dplyr::filter(all.annotated, !start2 %in% blacklist)
  
  
  ### filter out any duplicate row in `all.annotated`
  
  # remove duplicated rows (considering all columns; may happen when, e.g., both are oncogenes)
  all.annotated = unique(all.annotated)
  
  # give a unique marker to each non-duplicate row (don't consider gene, partner and fusion columns)
  all.annotated <- select(all.annotated, !c(entrez, gene, partner, fusion)) %>% 
    tidyr::unite(all_cols_marker, everything(), remove = TRUE) %>% 
    cbind(all.annotated, .)
  
  # get duplicate row markers
  dup_marker <- all.annotated$all_cols_marker [ duplicated(all.annotated$all_cols_marker) ] %>% 
    unique
  
  # remove duplicate rows that contain NA for partner column
  for(dup_marker_i in dup_marker){
    indx_dup_rows = which( all.annotated$all_cols_marker == dup_marker_i )
    partner_is_na = all.annotated$partner[indx_dup_rows] %>% 
      is.na
    if( any(partner_is_na) & any(!partner_is_na) ){
      all.annotated = -indx_dup_rows[partner_is_na] %>% 
        slice(all.annotated, .)
    }
  }
  
  # update duplicate row markers
  dup_marker <- all.annotated$all_cols_marker [ duplicated(all.annotated$all_cols_marker) ] %>% 
    unique
  
  # remove duplicate rows by prioritizing "gene" (oncogene) rather than "partner" 
  # these are the genes that have priority to be oncogene. firsts have higher priority. all others have the lowest priority.
  priority_to_be_oncogene <- c("MYC", "BCL6")
  for(dup_marker_i in dup_marker){
    indx_dup_rows = which( all.annotated$all_cols_marker == dup_marker_i )
    for(gene_i in priority_to_be_oncogene){
      is_oncogene <- all.annotated[indx_dup_rows,]$gene == gene_i
      if(any(is_oncogene))
        break
    }
    if( any(is_oncogene) & any(!is_oncogene) ){
      all.annotated = -indx_dup_rows[is_oncogene] %>% 
        slice(all.annotated, .)
    }
  }
  
  # check if there are still duplicate rows
  still_dup_rows <- all.annotated$all_cols_marker [ duplicated(all.annotated$all_cols_marker) ] %>% 
    unique %>% 
    {all.annotated$all_cols_marker %in% .} %>% 
    which
  if( length(still_dup_rows) > 0){
    paste(still_dup_rows, collapse = ", ") %>% 
      gettextf("Warning: There are still duplicate rows. They are: %s.", .) %>% 
      message  
  }
  
  # remove the row markers
  all.annotated = select(all.annotated, -all_cols_marker)
  
  
  if(return_as == "bedpe"){
    all.annotated$name = "."
    all.annotated = dplyr::select(all.annotated, chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, tumour_sample_id, gene, partner, fusion)
  }else if(return_as == "bedpe_entrez"){ #necessary for setting up cBioPOrtal instances (setup_study from portal.R)
    all.annotated$name = "."
    all.annotated = dplyr::select(all.annotated, chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, tumour_sample_id, gene, entrez, partner, fusion)
  }else if(return_as == "bed"){
    
    #lose the linkage but add a name that somewhat retains it
    if(!grepl("chr", all.annotated$chrom1)){
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom1 = paste0("chr", chrom1))
      
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom2 = paste0("chr", chrom2))
    }
    bed1 = dplyr::mutate(all.annotated, name = paste(tumour_sample_id, fusion, sep = "_")) %>%
      dplyr::select(chrom1, start1, end1, name, score, strand1)
    
    bed2 = dplyr::mutate(all.annotated, name = paste(tumour_sample_id, fusion, sep = "_")) %>%
      dplyr::select(chrom2, start2, end2, name, score, strand2)
    
    colnames(bed1) = c("chrom", "start", "end", "name", "score", "strand")
    colnames(bed2) = c("chrom", "start", "end", "name", "score", "strand")
    return(dplyr::arrange(rbind(bed1, bed2), name))
  }else{
    if(collapse_redundant){
      all.annotated = dplyr::distinct(all.annotated, tumour_sample_id, fusion, .keep_all = TRUE)
    }
  }
  if(with_chr_prefix){
    
    #add the prefix if necessary
    if(!grepl("chr", all.annotated$chrom1)){
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom1 = paste0("chr", chrom1))
      
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom2 = paste0("chr", chrom2))
    }
  }
  return(all.annotated)
}
