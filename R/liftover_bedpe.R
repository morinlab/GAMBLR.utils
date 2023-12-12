#' @title Liftover Bedpe.
#'
#' @description Use liftOver to convert a bedpe file between the two main genome builds (grch37/hg38).
#'
#' @details The user can specify a path to the bedpe file that needs to be lifted with `bedpe_file`,
#' or, the suer can specify the bedpe data in a data frame with `bedpe_df`.
#' The other required parameter is `target_build`, this parameter decides the final projection of the lifted bedpe file.
#'
#' @param bedpe_file Either specify the path to a bedpe file.
#' @param bedpe_df Or specify the bedpe data in a data frame.
#' @param target_build Specify which build the data should be lifted to (must be one of hg19, grch37, hg38, grch38).
#' @param col_names If not provided, the column names will be imposed.
#' @param col_types Specify column types if column names are also defined with `col_names`.
#' @param standard_bed Boolean parameter for defining the type of bed file that is provided with `bedpe_file`. Deafult is FALSE.
#' @param verbose Set to TRUE for verbose output. Default is FALSE.
#'
#' @return Data frame containing original bedpe data with new coordinates.
#'
#' @rawNamespace import(S4Vectors, except = c("merge", "second", "first", "union", "intersect", "setdiff", "setequal", "rename", "expand"))
#' @import dplyr tidyr readr rtracklayer GAMBLR.data
#' @export
#'
#' @examples
#' hg19_sv = get_manta_sv(verbose = FALSE)
#' hg19_sv = head(hg19_sv, 100)
#'
#' hg38_sv = liftover_bedpe(bedpe_df = hg19_sv,
#'                          target_build = "hg38")
#'
liftover_bedpe = function(bedpe_file,
                          bedpe_df,
                          target_build = "grch37",
                          col_names,
                          col_types,
                          standard_bed = FALSE,
                          verbose = FALSE){

  if(!missing(bedpe_file)){
    if(missing(col_names)){
      message("imposing column names")
      original_bedpe = suppressMessages(read_tsv(bedpe_file, comment = "##", col_types="cddcddccccccccccccccccc"))
    }else{
      message(paste("using column names", col_names, sep = ": "))
      original_bedpe = suppressMessages(read_tsv(bedpe_file, col_names = col_names, col_types = col_types))
    }
  }else{
    original_bedpe = bedpe_df
  }
  if(!standard_bed){
    colnames(original_bedpe)[1] = "CHROM_A"
    original_bedpe = as.data.frame(original_bedpe)
    original_bedpe = original_bedpe %>%
      dplyr::mutate(CHROM_A = ifelse(!grepl("chr", CHROM_A), paste0("chr", CHROM_A), CHROM_A),
                    CHROM_B = ifelse(!grepl("chr", CHROM_B), paste0("chr", CHROM_B), CHROM_B))
    #convert to strings manually to avoid caused by scientific notation in rare cases when R coerces to strings
    #Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec, :
    #scan() expected 'an integer', got '4.7e+07'

    original_bedpe = original_bedpe %>% mutate(START_A = format(START_A,scientific=F),
                                               START_B = format(START_B,scientific=F),
                                               END_A = format(END_A,scientific=F),
                                               END_B = format(END_B,scientific=F))

    if(verbose){
      print(head(original_bedpe))
    }

    char_vec = original_bedpe %>%
      tidyr::unite(united, sep = "\t") %>%
        dplyr::pull(united)

    bedpe_obj = rtracklayer::import(text = char_vec, format = "bedpe")
    if(length(colnames(original_bedpe)) > 22){
      this_patient = colnames(original_bedpe)[23]
      this_normal = colnames(original_bedpe)[22]
    }else{
      this_patient = original_bedpe$tumour_sample_id
      this_normal = original_bedpe$normal_sample_id
    }
  }else{
    colnames(original_bedpe)[1] = "chrom"
    if(!grepl("chr", original_bedpe$chrom)){
      original_bedpe = mutate(original_bedpe, chrom = paste0("chr", chrom))
    }
    char_vec = original_bedpe %>%
      tidyr::unite(united, sep = "\t") %>%
      dplyr::pull(united)

    bedpe_obj = rtracklayer::import(text = char_vec, format = "bed")
  }
  if(target_build == "grch37" | target_build == "hg19"){
    chain = rtracklayer::import.chain(system.file("extdata", "hg38ToHg19.over.chain", package = "GAMBLR.data"))
  }else if(target_build == "grch38" | target_build == "hg38"){
    chain = rtracklayer::import.chain(system.file("extdata", "hg19ToHg38.over.chain", package = "GAMBLR.data"))
  }
  if(!standard_bed){
    colnames(original_bedpe)[1] = "CHROM_A"
    original_columns = colnames(original_bedpe)

    first_sv_lifted = rtracklayer::liftOver(bedpe_obj@first, chain)
    second_sv_lifted = rtracklayer::liftOver(bedpe_obj@second, chain)
    no_problem = !((elementNROWS(first_sv_lifted) != 1) | (elementNROWS(second_sv_lifted) != 1))

    first_ok = subset(first_sv_lifted, no_problem)
    second_ok = subset(second_sv_lifted, no_problem)

    first_ok_df = data.frame(first_ok@unlistData) %>%
      dplyr::select(seqnames, start, end, strand) %>%
      dplyr::mutate(start = start - 1, end = as.numeric(end)) %>%
      `names<-`(c("CHROM_A", "START_A", "END_A", "STRAND_A"))

    second_ok_df = data.frame(second_ok@unlistData) %>%
      dplyr::select(seqnames, start, end, strand) %>%
      dplyr::mutate(start = start - 1, end = as.numeric(end)) %>%
      `names<-`(c("CHROM_B", "START_B", "END_B", "STRAND_B"))

    ok_bedpe = original_bedpe[no_problem,]

    kept_cols = ok_bedpe %>%
      dplyr::select(-c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "STRAND_A", "STRAND_B"))

    fully_lifted = bind_cols(first_ok_df, second_ok_df, kept_cols) %>%
      dplyr::select(all_of(original_columns)) %>%
      dplyr::mutate_if(is.factor, as.character)

    return(fully_lifted)

  }else{
    lifted = rtracklayer::liftOver(bedpe_obj, chain)
    no_problem = !((elementNROWS(lifted) != 1))
    first_ok = subset(lifted, no_problem)

    output = rtracklayer::export(first_ok, format = "bed") %>%
      suppressMessages(read_tsv(col_names = c("chrom", "start", "end", "score", "strand", "nothing", "s1", "e1", "junk", "more", "stuff", "nada"))) %>%
      dplyr::select("chrom", "start", "end")

    return(output)
  }
}
