#' @title Reannotate SSMs in a MAF format
#'
#' @description Pick another annotation for SSMs annotated by VEP using the additional isoforms
#'
#' @details To Do
#'
#'
#' @return A data frame in maf_data format
#'
#' @import tidyr dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#'\dontrun{
#' csnk2b_region_maf = get_ssm_by_regions(regions_list="chr6:31665606-31674128",
#'                                      these_samples_metadata = get_gambl_metadata(),
#'                                       projection="hg38",
#'                                       basic_columns=F,streamlined=F,this_seq_type = "genome")
#' csnk2b_region_maf_cap = get_ssm_by_regions(regions_list="chr6:31665606-31674128",
#'                                      these_samples_metadata = get_gambl_metadata(),
#'                                       projection="hg38",
#'                                       basic_columns=F,streamlined=F,this_seq_type = "capture")
#'}
reannotate_ssm = function(maf_data,
                          preferred_transcripts){
  #example:
  preferred_transcripts = c("ENST00000375865")
  preferred_transcripts = c("ENST00000375880")

  rename_variant_classification = function(effect){
    if(effect == "missense_variant"){
      return("Missense_Mutation")
    }
    if(effect == "stop_gained"){
      return("Nonsense_Mutation")
    }
    if(effect == "synonymous_variant"){
      return("Silent")
    }
    if(effect == "frameshift_variant"){
      return("Frame_Shift_Del")
    }
    if(grepl("splice",effect)){
      return("Splice_Site")
    }
    if(effect == "upstream_gene_variant"){
      return("5'Flank")
    }
    if(effect == "downstream_gene_variant"){
      return("3'Flank")
    }
    if(effect == "intron_variant"){
      return("Intron")
    }
    if(effect == "3_prime_UTR_variant"){
      return("3'UTR")
    }
    if(effect == "5_prime_UTR_variant"){
      return("5'UTR")
    }
    print(paste("Unrecognized:",effect))
  }
  shorten_hgvsp = function(effect){
      # Three-letter to one-letter amino acid code mapping
      aa_map <- c(
        Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
        Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
        Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
        Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
        Ter = "*"
      )

      # Pattern: "p.AAA###BBB" where AAA and BBB are 3-letter codes
      pattern <- "^p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2}|Ter).*"

      if (grepl(pattern, effect)) {
        parsed <- sub(pattern, "\\1;\\2;\\3", effect)
        parts <- strsplit(parsed, ";")[[1]]

        from <- aa_map[[parts[1]]]
        pos <- parts[2]
        to <- aa_map[[parts[3]]]

        if (is.null(from) || is.null(to)) {
          stop("Unrecognized amino acid abbreviation.")
        }

        return(paste0("p.", from, pos, to))
      } else if(effect == "" || effect == "p.="){
        return(effect)
      }else{
        print(effect)
        stop("Input does not match expected format (e.g., 'p.Gln182Ter').")
      }
  }

  pick_mutation = function(effects,transcript){
    muts = unlist(str_split(effects,";"))
    muts = muts[grepl(transcript,muts)]
    chunks = unlist(str_split(muts,","))
    Hugo_Symbol = chunks[1]
    Variant_Classification = rename_variant_classification(chunks[2])
    HGVSp_Short = shorten_hgvsp(chunks[3])
    #return corresponding Hugo_Symbol, Variant_Classification and HGVSp_Short
    #print(effects)
    #print(length(Variant_Classification))
    #print(length(HGVSp_Short))
    #print(paste(Hugo_Symbol,Variant_Classification,HGVSp_Short))
    return(data.frame(Hugo_Symbol=Hugo_Symbol,Variant_Classification=Variant_Classification,HGVSp_Short=HGVSp_Short))
  }
  xx = lapply(maf_data$all_effects,pick_mutation,preferred_transcripts)
  maf_cols_new = do.call("bind_rows",xx)
  maf_data$Hugo_Symbol = maf_cols_new$Hugo_Symbol
  maf_data$Variant_Classification = maf_cols_new$Variant_Classification
  maf_data$HGVSp_Short = maf_cols_new$HGVSp_Short
  CSNK2B_aliases = c("CSNK2B-LY6G5B-1181","XXbac-BPG32J3.22")
  maf_data = mutate(maf_data,Hugo_Symbol=ifelse(Hugo_Symbol %in% CSNK2B_aliases,"CSNK2B",Hugo_Symbol))
  return(maf_data)
}
