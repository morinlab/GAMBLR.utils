#' @title Maf To Custom Track.
#'
#' @description Convert a maf-formatted data frame into a bed custom track file for UCSC.
#'
#' @details This function takes an incoming MAF and converts it to a UCSC Genome Browser ready BED (or bigbed/biglolly) file.
#' Optional parameters available for further customization of the returned file. For more information, refer to the parameter descriptions and function examples.
#'
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function.
#' @param these_samples_metadata Optional argument, a metadata table subset to the samples of interest. If not provided, the function will return metadata for all available samples.
#' @param seq_type The seq type you want back, default is "genome".
#' @param output_file Name for your new bed file that can be uploaded as a custom track to UCSC.
#' @param as_bigbed Boolean parameter controlling the format of the returned file. Default is FALSE.
#' @param colour_column Set the colouring properties of the returned bed file. Per default, this function will assign colour based on "lymphgen".
#' @param as_biglolly Boolean parameter controlling the format of the returned file. Default is FALSE (i.e a BED file will be returned).
#' @param track_name Track name. Default is "GAMBL mutations"
#' @param track_description Track description. Default is "mutations from GAMBL"
#' @param verbose Default is FALSE.
#' @param padding_size Optional parameter specifying the padding size in the returned file, default is 0.
#'
#' @return Nothing.
#'
#' @import tidyr dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' my_maf <- get_coding_ssm()
#' maf_to_custom_track(maf_data = my_maf, output_file = "../mutations.bed")
#'
maf_to_custom_track = function(maf_data,
                               these_samples_metadata,
                               seq_type = "genome",
                               output_file,
                               as_bigbed = FALSE,
                               colour_column = "lymphgen",
                               as_biglolly = FALSE,
                               track_name = "GAMBL mutations",
                               track_description = "mutations from GAMBL",
                               verbose = FALSE,
                               padding_size = 0){

  #reduce to a bed-like format
  maf_data = dplyr::select(maf_data, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)
  colnames(maf_data) = c("chrom", "start", "end", "sample_id")
  maf_data = mutate(maf_data,end = end + padding_size)
  if(!any(grepl("chr", maf_data[,1]))){
    #add chr
    maf_data[,1] = unlist(lapply(maf_data[,1], function(x){paste0("chr", x)}))
  }
  lymphgen_cols = GAMBLR.helpers::get_gambl_colours(colour_column,verbose=verbose)

  colour_df = data.frame(group = names(lymphgen_cols), colour = lymphgen_cols)

  rgb_df = data.frame(t(col2rgb(lymphgen_cols))) %>%
    mutate(group = names(lymphgen_cols),hex=unname(lymphgen_cols)) %>%
    unite(col = "rgb", red, green, blue, sep = ",")
  if(verbose){
    print(rgb_df)
  }
  if(missing(these_samples_metadata)){
    meta = GAMBLR.helpers::handle_metadata(this_seq_type = seq_type) %>% dplyr::select(sample_id,all_of(colour_column))
  }else{
    meta = these_samples_metadata %>% dplyr::select(sample_id,all_of(colour_column))
  }
  colnames(meta)[2]="group"


  samples_coloured = left_join(meta, rgb_df)
  if(verbose){
    print(samples_coloured)
  }

  maf_bed = maf_data %>%
    mutate(score = 0, strand = "+", start1 = start-1,start=start1, end1 = end)
  if(verbose){
    print(head(maf_bed))
  }
  maf_coloured = left_join(maf_bed, samples_coloured, by = "sample_id") %>%
    dplyr::select(-group) %>%
    mutate(rgb=ifelse(is.na(rgb),"0,0,0",rgb))
  maf_summary = group_by(maf_coloured,hex) %>% tally()
  if(verbose){
    print(maf_summary)
    print(head(maf_coloured))
  }
  maf_coloured = dplyr::select(maf_coloured,-hex)
  if(as_bigbed | as_biglolly){

    if(grepl(pattern = ".bb$",x = output_file)){
      #temp file will be .bed
      temp_bed = gsub(".bb$",".bed",output_file)

    }else{
      stop("please provide an output file name ending in .bb to create a bigBed file")
    }

    maf_coloured = mutate(maf_coloured,sample_id="redacted") %>%
      arrange(chrom,start)
    if(as_biglolly){
      #currently the same code is run either way but this may change so I've separated this until we settle on format
      #TO DO: collapse based on hot spot definition and update column 4 (score) based on recurrence
      #needs to have size column
      maf_score_options = factor(maf_coloured$rgb)
      maf_coloured$score = as.numeric(maf_score_options)

      #determine frequency of each event per group to assign the size
      maf_coloured = group_by(maf_coloured,start,rgb) %>% mutate(size=n())

      #maf_coloured = mutate(maf_coloured,size=10)

      write.table(maf_coloured, file = temp_bed, quote = F, sep = "\t", row.names = F, col.names = F)
      #conversion:
      autosql_file = "/Users/rmorin/git/LLMPP/resources/reference/ucsc/bigLollyExample3.as"

      bigbedtobed = "/Users/rmorin/miniconda3/envs/ucsc/bin/bedToBigBed"
      bigbed_conversion = paste0(bigbedtobed," -as=",autosql_file," -type=bed9+1 ",temp_bed," /Users/rmorin/git/LLMPP/resources/reference/ucsc/hg19.chrom.sizes ",output_file)
      print(bigbed_conversion)
      system(bigbed_conversion)
    }else{
      write.table(maf_coloured, file = temp_bed, quote = F, sep = "\t", row.names = F, col.names = F)
      #conversion:
      bigbedtobed = "/Users/rmorin/miniconda3/envs/ucsc/bin/bedToBigBed"
      bigbed_conversion = paste(bigbedtobed,"-type=bed9",temp_bed,"/Users/rmorin/git/LLMPP/resources/reference/ucsc/hg19.chrom.sizes",output_file)

      system(bigbed_conversion)
    }
  }else{
    header_ucsc = paste0('track name="',track_name,'" description="', track_description, '" visibility=2 itemRgb="On"\n')
    cat(header_ucsc,file = output_file)
    write.table(maf_coloured, file = output_file, quote = F, sep = "\t", row.names = F, col.names = F, append = TRUE)
  }
}
