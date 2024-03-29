#' @title MAF to ProteinPaint table
#' 
#' @description This function takes a MAF-like data frame and convert it to the right format 
#' for visualization with ProteinPaint.
#' 
#' @details For visualization with ProteinPaint, the output data frame must be saved to a 
#' file and uploaded to https://proteinpaint.stjude.org/. 
#' 
#' For a valid ProteinPaint table for visualization, the MAF table must contain the required 
#' columns to be converted to the ProteinPaint format. They are: Hugo_Symbol, RefSeq, 
#' Chromosome, Start_Position, HGVSp_Short, Variant_Classification. Optional columns used 
#' from the MAF table are: Mutation_Status, t_alt_count, t_depth, n_alt_count, n_depth, 
#' Tumor_Sample_Barcode, Reference_Allele. In the case of using functions `get_ssm_by_samples`, 
#' `get_coding_ssm`, or `get_ssm_by_patients` to get your MAF table and required/desired 
#' columns are missing, consider using parameter `basic_columns = FALSE` to ensure that your 
#' MAF has all needed columns. Other columns are taken from the metadata provided by the
#' `these_samples_metadata` parameter or generated internally.
#' 
#' @param maf_data A data frame in MAF format.
#' @param these_sample_ids A vector of sample IDs to subset the samples of interest from the 
#' input `maf_data`. If NULL (the default), samples are not filtered by this parameter.
#' @param these_samples_metadata A metadata table, with sample IDs in a column, to subset 
#' the samples of interest from the input `maf_data`. Also, important information from the 
#' metadata is added to the output table. If NULL (the default), `ssm_to_proteinpaint` 
#' internally creates the metadata table with samples according to the parameters 
#' `these_sample_ids` and `this_seq_type`.
#' @param this_seq_type The seq type to get results for. One of "genome" (default) or 
#' "capture".
#' @param sample_type A single sting with the name of the optional column (from `maf_data` 
#' or the metadata) that you want to use to distinguish multiple samples from the same 
#' patient. If NULL, the ProteinPaint visualization is split only by pathologies. If not 
#' NULL, the visualization will be split by each combination between pathologies and the 
#' column specified by `sample_type`. The default is "time_point".
#' @param coding_only Boolean parameter. Set to TRUE to restrict to only coding mutations.
#' The default is FALSE.
#' @param these_genes A vector of strings with name of genes that you want results for. 
#' If NULL, all genes of the input `maf_data` are kept. The default is all lymphoma genes. 
#' @param debug_flag Boolean parameter. Set to TRUE for returning rows from the 
#' incoming MAF that do not contain any values in the required columns. Commonly used for 
#' checking purposes only. Setting this to TRUE, does not produce an output compatible with 
#' Protein paint. The default is FALSE.
#'
#' @return A data frame.
#' 
#' @import tidyr
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' 
#' # define same somples
#' my_samples = c("DOHH-2", "OCI-Ly10", "OCI-Ly3", "SU-DHL-10", "SU-DHL-4")
#' 
#' # get maf data frame
#' my_maf = get_ssm_by_samples(these_sample_ids = my_samples, 
#'                             this_seq_type = "genome", basic_columns = FALSE)
#' 
#' # convert maf to ProteinPaint format
#' pp_df = ssm_to_proteinpaint(maf_data = my_maf, this_seq_type = "genome")
#' 
#' # convert only coding mutations to ProteinPaint format
#' pp_df = ssm_to_proteinpaint(maf_data = my_maf, this_seq_type = "genome",
#'                             coding_only = TRUE)
#' 
ssm_to_proteinpaint = function(maf_data,
                               these_sample_ids = NULL,
                               these_samples_metadata = NULL,
                               this_seq_type = "genome",
                               sample_type = "time_point",
                               coding_only = FALSE,
                               these_genes = GAMBLR.data::lymphoma_genes_comprehensive$Gene,
                               debug_flag = FALSE){
  
  # check parameters
  stopifnot( "`this_seq_type` must be one of \"genome\" or \"capture\"." = 
               this_seq_type %in% c("genome", "capture") & length(this_seq_type) == 1 )
  
  if(!is.null(sample_type)){
    stopifnot( "`sample_type` must be either NULL or a single string." = 
                 is.character(sample_type) & length(sample_type) == 1 )
  }
  
  # check for required columns in maf_data 
  maf_req_cols = c("Hugo_Symbol", "RefSeq", "Chromosome", "Start_Position", "HGVSp_Short", 
                    "Variant_Classification")
  maf_req_cols_not_present = maf_req_cols %>% 
    "["(! . %in% names(maf_data))
  if( length(maf_req_cols_not_present) > 0 ){
    k = paste(maf_req_cols_not_present, collapse = ", ") %>% 
      gettextf("Required columns missing in the input MAF data frame: %s.", .)
    stop(k)
  }
  
  # get metadata with the dedicated helper function
  these_samples_metadata = id_ease(
    these_samples_metadata = these_samples_metadata,
    these_sample_ids = these_sample_ids,
    this_seq_type = this_seq_type,
    verbose = FALSE
  )
  
  # filter maf according to the samples in the metadata
  maf_data = dplyr::filter(maf_data, Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)
  
  # keep only those ssms in lymphoma genes
  if( ! is.null(these_genes) ){
    maf_data = filter(maf_data, Hugo_Symbol %in% these_genes)
  }
  
  # add metadata columns to maf_data
  maf_data = left_join(maf_data, these_samples_metadata, by = join_by(Tumor_Sample_Barcode == sample_id))
  
  # check for optional columns in maf_data (columns from the metadata included)
  maf_opt_cols = c("Mutation_Status", "t_alt_count", "t_depth", "n_alt_count", "n_depth", 
                   "Tumor_Sample_Barcode", "Reference_Allele", "pathology", "patient_id",
                   sample_type)
  maf_opt_cols_not_present = maf_opt_cols %>% 
    "["(! . %in% names(maf_data))
  if( length(maf_opt_cols_not_present) > 0 ){
    k = paste(maf_opt_cols_not_present, collapse = ", ") %>% 
      gettextf("Warning: Optional columns missing in the input MAF data frame or in the metadata and will be ignored: %s.", .)
    message(k)
  }
  
  # adequate maf_data to ProteinPaint visualization format
  maf_data = rename(maf_data, gene = Hugo_Symbol, chromosome = Chromosome, 
                    start = Start_Position, class = Variant_Classification,
                    refseq = RefSeq, aachange = HGVSp_Short, origin = Mutation_Status,
                    disease = pathology, mutant_in_tumor = t_alt_count, 
                    total_in_tumor = t_depth, mutant_in_normal = n_alt_count, 
                    total_in_normal = n_depth, patient = patient_id, 
                    sample = Tumor_Sample_Barcode, REF = Reference_Allele,
                    sampletype = any_of(sample_type)) %>% 
    mutate(class = recode(class,
                          Missense_Mutation = "missense",
                          In_Frame_Del = "proteinDel",
                          In_Frame_Ins = "proteinIns",
                          Nonsense_Mutation = "nonsense",
                          Nonstop_Mutation = "nonsense",
                          Silent = "silent",
                          Splice_Region = "splice_region",
                          Splice_Site = "splice",
                          Frame_Shift_Del = "frameshift",
                          Frame_Shift_Ins = "frameshift",
                          `5'UTR` = "utr_5",
                          `3'UTR` = "utr_3",
                          Intron = "intron",
                          `3'Flank` = "noncoding",
                          `5'Flank` = "noncoding",
                          IGR = "noncoding"))
  
  # if only coding mutations are desired
  if(coding_only){
    noncoding_muts = c("utr_5", "utr_3", "intron", "noncoding")
    maf_data = dplyr::filter(maf_data, !(class %in% noncoding_muts))
  }
  
  # remove version names in refseq column
  maf_data$refseq = maf_data$refseq %>% 
    gsub("\\.[0-9]+,", ",", .) %>% 
    gsub("\\.[0-9]+$", "", .)
  
  # if alt allele is heterogeneous alternative, show both allele separated 
  # by comma (like in VCF files)
  is_het_alt = ! maf_data$REF == maf_data$Tumor_Seq_Allele1
  maf_data$ALT = is_het_alt %>% 
    { paste( maf_data$Tumor_Seq_Allele1[.], maf_data$Tumor_Seq_Allele2[.], sep = "," ) } %>% 
    replace(maf_data$Tumor_Seq_Allele2, is_het_alt, .)
  
  # add VAF (variant allele fraction)
  maf_data = mutate(maf_data, VAF = mutant_in_tumor/total_in_tumor)
  
  # keep only columns that are meaningful for proteinpaint
  if( !is.null(sample_type) ){
    sample_type = "sampletype"
  }
  maf_data = dplyr::select(maf_data, gene, refseq, chromosome, start, aachange, class, 
                           disease, origin, patient, sample, any_of(sample_type), mutant_in_tumor, 
                           total_in_tumor, mutant_in_normal, total_in_normal, REF, ALT, VAF)
  
  # if aachange is empty, change it to NA. this is to avoid filtering out desirable rows
  maf_data = (maf_data$aachange == "") %>% 
    { mutate(maf_data, aachange = ifelse(., "NA", aachange)) }
  
  # to filter out rows that don't contain required column values and output warning message
  required_cols = dplyr::select(maf_data, gene, refseq, chromosome, start, aachange, class) %>% 
    mutate(i = row_number()) %>% 
    tidyr::drop_na()
  required_cols = apply(required_cols != "", 1, all) %>% 
    dplyr::filter(required_cols, .)
  
  # to filter out rows with unsupported values in column class
  suported_class = c("missense", "proteinDel", "proteinIns", "nonsense", "silent", 
                     "splice_region", "splice", "frameshift", "utr_5", "utr_3", "intron", 
                     "noncoding")
  required_cols = dplyr::filter(required_cols, class %in% suported_class)
  
  # if only removed rows should be returned (for checking purpose)
  if(debug_flag){  
    removed_rows = 1:nrow(maf_data) %>% 
      "["(! . %in% required_cols$i) %>% 
      slice(maf_data, .)
    return(removed_rows)
  }
  
  # print a message stating how many mutations were lost
  num_removed_rows = (nrow(maf_data) - nrow(required_cols)) %>% 
    format(big.mark="'")
  k = format( nrow(maf_data), big.mark="'" ) %>% 
    gettextf("Warning: %s rows out of %s were removed from the output table because there were missing or unsupported values in required columns. Run `ssm_to_proteinpaint` again with `debug_flag = TRUE` to see these rows.",
             num_removed_rows, .)
  message(k)
  maf_data = slice(maf_data, required_cols$i)
  
  # return final output table
  maf_data = as.data.frame(maf_data)
  maf_data
}
