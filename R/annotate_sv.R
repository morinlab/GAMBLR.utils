#' @title Annotate SVs.
#'
#' @description Annotate a data frame of SV breakpoints represented in an extended BEDPE format.
#'
#' @details Specify a data frame with SVs (preferably the output from `GAMBLR.results::get_manta_sv` 
#' if the user has GSC-restricted access, or `GAMBLR.data::sample_data$grch37$bedpe` or 
#' `GAMBLR.data::sample_data$hg38$bedpe` otherwise) to the `sv_df` parameter and get back the same 
#' data frame with SV annotations. Alternatively, you can also provide your own data frame to the 
#' `sv_data` argument, which should contain the following columns: CHROM_A, START_A, END_A, CHROM_B, 
#' START_B, END_B, NAME, SOMATIC_SCORE, STRAND_A, STRAND_B, TYPE, FILTER, VAF_tumour, VAF_normal, 
#' DP_tumour, DP_normal, tumour_sample_id, normal_sample_id, pair_status.
#' 
#' Before printing the SV-annotated data frame, the function uses the `priority_to_be_oncogene` 
#' parameter to filter rows. If there is any pair of rows where (1) the values of their columns &mdash; 
#' other than "gene", "partner" and "fusion" &mdash; are the same, (2) the "gene" of a row is the 
#' "partner" of the other and vice versa, and (3) at least one of the genes is listed in 
#' `priority_to_be_oncogene`, the function keeps only that row where the "gene" column contains the 
#' gene with the highest priority to be oncogene (see how priority is given in the description of 
#' this parameter). If neither of the rows contains a gene listed in `priority_to_be_oncogene` in 
#' their "gene" column, both rows are kept and a warning message followed by the row values is printed.
#'
#' @param sv_data A data frame of SVs. This should be the output of get_manta_sv. If you aren't using the database backend you can supply your own data frame (see format in **Details** section).
#' Most of this data is directly from the bedpe files that are obtained by converting the Manta outputs from VCF.
#' Only the columns for both chromosomes, coordinates and strand plus SOMATIC_SCORE and tumour_sample_id are absolutely required.
#' @param partner_bed Optional bed-format data frame to use for annotating oncogene partners (e.g. enhancers). required columns are: chrom,start,end,gene
#' @param with_chr_prefix Optionally request that chromosome names are returned with a chr prefix. Default is FALSE.
#' @param collapse_redundant Remove reciprocal events and only return one per event. Default is FALSE.
#' @param return_as Stated format for returned output, default is "bedpe". Other accepted output formats are "bed" and "bedpe_entrez" (to keep entrez_ids for compatibility with portal.R and cBioPortal).
#' @param blacklist A vector of regions to be removed from annotations. Default coordinates are in respect to hg19.
#' @param genome_build Reference genome build parameter, default is grch37.
#' @param priority_to_be_oncogene Vector of gene names (default is `c("MYC", "BCL6")`) used to filter rows based on genes that have the highest priority to be considered oncogenes. Genes to the left (*i.e.* first elements of this vector) have higher priority; non-listed genes have the lowest priority. See **Details** section for more information.
#' @param oncogene_bed Optionally, provide bed for regions to be considered oncogenes.
#' @param include_shm_partners Set to TRUE if you want to allow aSHM regions to be considered valid oncogene partners when annotating SVs
#' 
#' @return A data frame with annotated SVs (gene symbol and entrez ID).
#' 
#' @import tidyr dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' 
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
                       blacklist = c(60565248, 30303126, 187728894, 101357565, 101359747, 161734970, 69400840, 65217851, 187728889, 187728888,187728892, 187728893,188305164, 72551424, 72551425, 72551554, 72551555, 72551558, 72551559, 72551562, 189255425, 189255426, 189255461, 189255462),
                       genome_build = "grch37",
                       priority_to_be_oncogene = c("MYC", "BCL6"),
                       oncogene_bed,
                       include_shm_partners = FALSE){
  if(!"SCORE" %in% colnames(sv_data)){
    stop("input is missing required column 'SCORE'")
  }
  # remove duplicate rows in sv_data, if any
  sv_data = unique(sv_data)
  
  bedpe1 = sv_data %>%
    dplyr::select(
        "CHROM_A", "START_A", "END_A", "tumour_sample_id",
        ends_with("SCORE"), "STRAND_A",
        if ("ID" %in% colnames(sv_data)) "ID" else "manta_name"
    )
  
  bedpe2 = sv_data %>%
    dplyr::select(
        "CHROM_B", "START_B", "END_B", "tumour_sample_id",
        ends_with("SCORE"), "STRAND_B",
        if ("ID" %in% colnames(sv_data)) "ID" else "manta_name"
    )
  
  colnames(bedpe1) = c("chrom", "start", "end", "tumour_sample_id", "score", "strand1", "ID")
  colnames(bedpe2) = c("chrom", "start", "end", "tumour_sample_id", "score", "strand2", "ID")
  suppressWarnings({
    if(any(grepl("chr", bedpe1$chrom))){
      bedpe1 = dplyr::mutate(bedpe1, chrom = gsub("chr", "", chrom))
      bedpe2 = dplyr::mutate(bedpe2, chrom = gsub("chr", "", chrom))
    }
  })
  if(missing(partner_bed)){
    if(genome_build == "hg38"){
      ig_regions = GAMBLR.data::hg38_partners
      if(include_shm_partners){
        extra_regions = GAMBLR.data::somatic_hypermutation_locations_GRCh38_v_latest %>% 
          filter(!gene %in% ig_regions$gene) %>%
          dplyr::select(chr_name,hg38_start,hg38_end,gene) %>%
          mutate(chr_name=str_remove(chr_name,"chr"))
        extra_regions = extra_regions %>% 
          mutate(start=ifelse(hg38_start<hg38_end,hg38_start,hg38_end)) %>%
          mutate(end=ifelse(hg38_start<hg38_end,hg38_end,hg38_start)) %>%
          mutate(chrom = chr_name) %>%
          dplyr::select(chrom,start,end,gene) %>%
          mutate(start=start-5000,end=end+5000)
        
        extra_regions$entrez = 0
        ig_regions = bind_rows(extra_regions,ig_regions)
      }
    }else{
      ig_regions = GAMBLR.data::grch37_partners 
      
      if(include_shm_partners){
       extra_regions = GAMBLR.data::somatic_hypermutation_locations_GRCh37_v_latest %>% 
        filter(!gene %in% ig_regions$gene) %>%
        dplyr::select(chr_name,hg19_start,hg19_end,gene) %>%
         mutate(chr_name=str_remove(chr_name,"chr"))
       extra_regions = extra_regions %>% 
         mutate(start=ifelse(hg19_start<hg19_end,hg19_start,hg19_end)) %>%
         mutate(end=ifelse(hg19_start<hg19_end,hg19_end,hg19_start)) %>%
         mutate(chrom = chr_name) %>%
         dplyr::select(chrom,start,end,gene) %>%
         mutate(start=start-5000,end=end+5000)
       
       extra_regions$entrez = 0
       ig_regions = bind_rows(extra_regions,ig_regions) %>% arrange(chrom)
      }
    }
  }else{
    ig_regions = partner_bed
    if(!"entrez" %in% colnames(ig_regions)){
      ig_regions$entrez = 0
    }
  }
  if(missing(oncogene_bed)){
    if(genome_build == "hg38"){
      oncogene_regions = GAMBLR.data::hg38_oncogene
    }else{
      oncogene_regions = GAMBLR.data::grch37_oncogene
    }
  }else{
    oncogene_regions = oncogene_bed
  }

  y = as.data.frame(oncogene_regions)
  
  #use overlaps to get oncogene SVs
  a = as.data.frame(bedpe1)
  a.onco = cool_overlaps(a, y, columns1 = c("chrom", "start", "end"), columns2 = c("chrom", "start", "end"), nomatch = TRUE) #oncogene-annotated bedpe for the first breakpoints
  b = as.data.frame(bedpe2)
  b.onco = cool_overlaps(b, y, columns1 = c("chrom", "start", "end"), columns2 = c("chrom", "start", "end"), nomatch = TRUE) #oncogene-annotated bedpe for the first breakpoints

  #insist oncogene breakpoints are anchored in an IG or superenhancer region (currently just IG or BCL6)
  #other end of breakpoint
  a.onco.break = a.onco[which(!is.na(a.onco$start.y)), c("chrom", "start", "end", "tumour_sample_id", "gene", "entrez", "score", "strand1")]
  b.onco.break = b.onco[which(!is.na(b.onco$start.y)), c("chrom", "start", "end", "tumour_sample_id", "gene", "entrez", "score", "strand2")]
  
  a.partner = b[which(!is.na(a.onco$start.y)),]
  b.partner = a[which(!is.na(b.onco$start.y)),]
  
  y = as.data.frame(ig_regions)
  
  a.ig = cool_overlaps(a.partner, y, type = "any", columns1 = c("chrom", "start", "end"), columns2 = c("chrom", "start", "end"), nomatch = TRUE)
  b.ig = cool_overlaps(b.partner, y, type = "any", columns1 = c("chrom", "start", "end"), columns2 = c("chrom", "start", "end"), nomatch = TRUE)

  a.ig <- a.ig %>%
    mutate(row_id = row_number()) %>% # make sure the order is consistent
    group_by(chrom, start, end, tumour_sample_id, score, ID) %>%
    slice_head() %>%
    ungroup %>%
    arrange(row_id) %>%
    as.data.frame() %>%
    select(-row_id)

  b.ig <- b.ig %>%
    mutate(row_id = row_number()) %>% # make sure the order is consistent
    group_by(chrom, start, end, tumour_sample_id, score, ID) %>%
    slice_head() %>%
    ungroup %>%
    arrange(row_id) %>%
    as.data.frame() %>%
    select(-row_id)

  a.ig = a.ig[,c("chrom", "start", "end", "strand2", "gene")]
  b.ig = b.ig[,c("chrom", "start", "end", "strand1", "gene")]

  a.annotated.both = bind_cols(a.onco.break, a.ig)
  colnames(a.annotated.both) = c("chrom1", "start1", "end1", "tumour_sample_id", "gene", "entrez", "score", "strand1", "chrom2", "start2", "end2", "strand2", "partner")

  b.annotated.both = cbind(b.onco.break, b.ig)
  colnames(b.annotated.both) = c("chrom2", "start2", "end2", "tumour_sample_id", "gene", "entrez", "score", "strand2", "chrom1", "start1", "end1", "strand1", "partner")

  all.annotated = rbind(a.annotated.both, b.annotated.both) %>%
    as.data.frame %>%
    `rownames<-`(NULL)

  all.annotated$fusion = dplyr::pull(tidyr::unite(all.annotated, fusion, partner, gene, sep = "-"), fusion)
  all.annotated = dplyr::filter(all.annotated, !fusion %in% c("BCL6-BCL6", "CIITA-CIITA", "FOXP1-FOXP1","PAX5-PAX5"))
  
  all.annotated = dplyr::filter(all.annotated, !start1 %in% blacklist)
  all.annotated = dplyr::filter(all.annotated, !start2 %in% blacklist)
  
  ### filter out any duplicate row in `all.annotated`
  
  # remove duplicated rows (considering all columns; may happen when, e.g., both are oncogenes)
  all.annotated = unique(all.annotated)

  # give a unique marker to each non-duplicate row (don't consider gene, partner and fusion columns)
  all.annotated = select(all.annotated, !c(entrez, gene, partner, fusion)) %>% 
    tidyr::unite(all_cols_marker, everything(), remove = TRUE) %>% 
    cbind(all.annotated, .)
  
  # get duplicate row markers
  dup_marker = all.annotated$all_cols_marker [ duplicated(all.annotated$all_cols_marker) ] %>% 
    unique
  
  # add column that stores row numbers, it will be used later to rearrange the rows
  all.annotated = mutate(all.annotated, row_num = row_number())
  
  # from all rows, separate those with any duplicate marker
  all.annotated = (all.annotated$all_cols_marker %in% dup_marker) %>% 
    ifelse("duplicate", "not_duplicate") %>% 
    factor(levels = c("duplicate", "not_duplicate")) %>% 
    split(all.annotated, .)
  
  # fix duplicated rows (rows with the same marker)
  if( nrow(all.annotated$duplicate) > 0 ){
    # separate duplicate rows by markers
    all.annotated$duplicate = split(all.annotated$duplicate, dup_marker)
    
    # check whether duplications are in pairs
    stopifnot( "There are duplicate rows that are not in pairs." = all( lapply(all.annotated$duplicate, nrow) == 2 ) )
    
    # if both genes are oncogenes and partners, remove row according to gene prioritization of being oncogene.
    # `priority_to_be_oncogene` contains the genes with high priority to be oncogene. genes to the left have higher priority. all non-listed genes have the lowest priority.
    # if both genes have the lowest priority of being oncogene, no rows are removed, but reported.
    all.annotated$duplicate = lapply(all.annotated$duplicate, function(dup_pair){
      for(gene_i in priority_to_be_oncogene){
        is_priority_gene = dup_pair$gene == gene_i
        if(any(is_priority_gene))
          break
      }
      switch(
        as.character( sum(is_priority_gene) ),
        "0" = {
          # both genes have the lowest priority to be oncogene
          # if they are partner of each other, report this error
          if( all( !is.na(dup_pair$partner) ) ){
            message("Warning: There are two rows with the same values, but inverted ones for the `gene` and `partner` columns (i.e. both genes as oncogene and partner of each other). However, neither of them is listed as a gene with high priority to be oncogene (included in the `priority_to_be_oncogene` object). See these rows below:")
            select(dup_pair, -all_cols_marker) %>% 
              print
            message("Despite that, both rows were kept.\n")
          }
          dup_pair
        },
        "1" = {
          # one gene with a higher priority
          # if both genes are partners of each other, remove the row with the lower priority gene as oncogene
          if( all( !is.na(dup_pair$partner) ) ){
            filter(dup_pair, is_priority_gene)
          }else{
            dup_pair
          }
        },
        {
          message("Error: `as.character(sum(is_priority_gene))` when `gene_i` is the data frame")
          select(dup_pair, -all_cols_marker) %>% 
            print
          k = as.character(sum(is_priority_gene)) %>% 
            gettextf("\r is \"%s\", but it should be either \"0\" or \"1\"." , .)
          stop(k)
        }
      )
    }) %>% 
      bind_rows
  }
  
  # merge the not_duplicate and duplicate data frames
  all.annotated = bind_rows(all.annotated)
  
  # rearrange rows
  all.annotated = arrange(all.annotated, row_num)
  
  # remove temporary columns
  all.annotated = select(all.annotated, -all_cols_marker, -row_num)
  
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
    k = !grepl("chr", all.annotated$chrom1)
    if( all(k) ){
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom1 = paste0("chr", chrom1))
      
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom2 = paste0("chr", chrom2))
    }
  }
  all.annotated <- all.annotated %>%
    group_by(chrom1, start1, end1, chrom2, start2, end2) %>%
    # Keep rows where partner is not NA if both NA and value exist
    filter(
        !(is.na(partner) & any(!is.na(partner)))
    ) %>%
    ungroup() %>%
    as.data.frame
  return(all.annotated)
}
