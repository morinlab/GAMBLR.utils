#' @title Maf To Custom Track.
#'
#' @description Convert a maf-formatted data frame into a bed custom track file for UCSC.
#'
#' @details This function takes an incoming MAF and converts it to a UCSC Genome Browser ready BED (or bigbed/biglolly) file.
#' Optional parameters available for further customization of the returned file. For more information, refer to the parameter descriptions and function examples.
#'
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function.
#' @param these_sample_ids A vector of sample IDs to subset the samples of interest from the input `maf_data`. If NULL (the default), all samples in `maf_data` are kept.
#' @param these_samples_metadata A metadata table to subset the samples of interest from the input `maf_data`. If NULL (the default), all samples in `maf_data` are kept.
#' @param this_seq_type The seq type you want back, default is "genome".
#' @param output_file Name for your new bed file that can be uploaded as a custom track to UCSC.
#' @param as_bigbed Boolean parameter controlling the format of the returned file. Default is FALSE.
#' @param colour_column Set the colouring properties of the returned bed file. Per default, this function will assign colour based on "lymphgen".
#' @param as_biglolly Boolean parameter controlling the format of the returned file. Default is FALSE (i.e a BED file will be returned).
#' @param track_name Track name. Default is "GAMBL mutations"
#' @param track_description Track description. Default is "mutations from GAMBL"
#' @param verbose Default is FALSE.
#' @param padding_size Optional parameter specifying the padding size in the returned file, default is 0.
#' @param projection Specify which genome build to use. Possible values are "grch37" (default) or "hg38". This parameter has 
#'   effect only when `as_bigbed` or `as_biglolly` is TRUE.
#' @param bedToBigBed_path Path to your local `bedToBigBed` UCSC tool or the string 
#'   `"config"` (default). If set to `"config"`, `GAMBLR.helpers::check_config_value` 
#'   is called internally and the `bedToBigBed` path is obtained from the `config.yml` 
#'   file saved in the current working directory. This parameter is ignored if both
#'   `as_bigbed` and `as_biglolly` is set to `FALSE`. 
#'
#' @return Nothing.
#'
#' @import tidyr dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' region_myc <- dplyr::filter(grch37_lymphoma_genes_bed, hgnc_symbol == "MYC")
#' my_maf <- get_ssm_by_regions(regions_bed = region_myc, streamlined = FALSE)
#' maf_to_custom_track(maf_data = my_maf, output_file = "../mutations.bed")
#'
maf_to_custom_track = function(maf_data,
                               these_sample_ids = NULL,
                               these_samples_metadata = NULL,
                               this_seq_type = "genome",
                               output_file,
                               as_bigbed = FALSE,
                               colour_column = "lymphgen",
                               as_biglolly = FALSE,
                               track_name = "GAMBL mutations",
                               track_description = "mutations from GAMBL",
                               verbose = FALSE,
                               padding_size = 0,
                               projection = "grch37",
                               bedToBigBed_path = "config"){
  
  # check some provided parameter 
  if(bedToBigBed_path == "config"){
    bedToBigBed_path = tryCatch(
      GAMBLR.helpers::check_config_value(config::get("dependencies")$bedToBigBed),
      error=function(e){
        k = paste0("You set bedToBigBed_path parameter to \"config\". However...\n", e)
        stop(k, call. = FALSE)
      }
    )
  }else{
    stopifnot("`bedToBigBed_path` points to a non-existent file." = 
                file.exists(bedToBigBed_path))
  }
  
  # get metadata with the dedicated helper function
  these_samples_metadata = id_ease(
    these_samples_metadata = these_samples_metadata,
    these_sample_ids = these_sample_ids,
    this_seq_type = this_seq_type,
    verbose = FALSE
  )
  
  # subset input maf according to metadata samples
  maf_data = filter(maf_data, Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)
  
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
  meta = dplyr::select(these_samples_metadata, sample_id, all_of(colour_column))
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
      temp_bed = tempfile(pattern = "bed_")
      
    }else{
      stop("please provide an output file name ending in .bb to create a bigBed file")
    }
    
    maf_coloured = mutate(maf_coloured,sample_id="redacted") %>%
      arrange(chrom,start)
    
    # create temp file chrom.sizes
    if(projection == "grch37"){
      chr_arms <- GAMBLR.data::chromosome_arms_grch37 %>% 
        mutate(chromosome = paste0("chr", chromosome))
    }else if(projection == "hg38"){
      chr_arms <- GAMBLR.data::chromosome_arms_hg38
    }else{
      stop("projection parameter must be \"grch37\" or \"hg38\".")
    }
    chr_sizes <- chr_arms %>%
      dplyr::filter(arm == "q") %>%
      dplyr::select(chromosome, end) %>%
      rename(size = "end")
    temp_chr_sizes <- tempfile(pattern = "chrom.sizes_")
    write.table(chr_sizes, file = temp_chr_sizes, quote = FALSE, row.names = FALSE, 
                col.names = FALSE, sep = "\t")
    
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
      
      bigbed_conversion = paste(bedToBigBed_path, "-as=", autosql_file, "-tab -type=bed9+1", temp_bed,
                                temp_chr_sizes, output_file)
      print(bigbed_conversion)
      system(bigbed_conversion)
    }else{
      write.table(maf_coloured, file = temp_bed, quote = F, sep = "\t", row.names = F, col.names = F)
      #conversion:
      bigbed_conversion = paste(bedToBigBed_path, "-tab -type=bed9", temp_bed, temp_chr_sizes, output_file)
      
      system(bigbed_conversion)
    }
    unlink(c(temp_bed, temp_chr_sizes))
  }else{
    header_ucsc = paste0('track name="',track_name,'" description="', track_description, '" visibility=2 itemRgb="On"\n')
    cat(header_ucsc,file = output_file)
    write.table(maf_coloured, file = output_file, quote = F, sep = "\t", row.names = F, col.names = F, append = TRUE)
  }
}
