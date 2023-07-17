#' @title Estimate Purity.
#'
#' @description Annotate a MAF with segmented absolute copy number data and added additional columns (VAF, Ploidy and Final_purity).
#'
#' @details This function takes a sample ID with the `this_sample_id` parameter and calls [GAMBLR::assign_cn_to_ssm] to get CN information.
#' The user can also use an already loaded maf file with `maf_df`. In addition, a path to the maf/seq file of interest can also be passed to this function with
#' `in_maf` and `in_seg`. To visualize VAF and purity distributions, set the `show_plots` to TRUE (default is FALSE).
#' For more information on how to run this function with the parameters at hand, refer to the parameter descriptions and function examples.
#'
#' @param in_maf Path to a local maf file.
#' @param maf_df Optional. Instead of using the path to a maf file, use a local dataframe as the maf file.
#' @param in_seg Path to a local corresponding seg file for the same sample ID as the input maf.
#' @param this_sample_id Specify the sample_id or any other string you want embedded in the file name.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, battenberg, etc.
#' @param show_plots Optional. Show two faceted plots that display the VAF and purity distributions for each copy number state in the sample. Default is FALSE.
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral. Default is FALSE.
#' @param coding_only Optional. set to TRUE to restrict to only coding variants. Default is FALSE.
#' @param genes Genes of interest.
#'
#' @return A list containing a data frame (MAF-like format) with the segmented absolute copy number data and three extra columns:
#' VAF is the variant allele frequency calculated from the t_ref_count and t_alt_count
#' Ploidy is the number of copies of an allele in the tumour cell
#' Final_purity is the finalized purity estimation per mutation after considering different copy number states and LOH events.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' #load a maf
#' this_maf = get_ssm_by_sample(this_sample_id = "HTMCP-01-06-00422-01A-01D",
#'                              this_seq_type = "genome")
#'
#' #estimate purity based on an already loaded maf object
#' estimate_purity(maf_df = this_maf,
#'                 show_plots = TRUE)
#'
#' #estimate purity based sole on a smaple ID + added seg data.
#' estimate_purity(this_sample_id = "HTMCP-01-06-00422-01A-01D",
#'                 show_plots = TRUE, coding_only = TRUE)
#'
estimate_purity = function(in_maf,
                           maf_df,
                           in_seg,
                           this_sample_id,
                           this_seq_type = "genome",
                           seg_file_source = "battenberg",
                           show_plots = FALSE,
                           assume_diploid = FALSE,
                           coding_only = FALSE,
                           genes){

  # Merge the CN info to the corresponding MAF file, uses GAMBLR function
  if(missing(in_maf) & missing(in_seg) & missing(maf_df)){
    CN_new = assign_cn_to_ssm(
      this_sample_id = this_sample_id,
      coding_only = coding_only,
      assume_diploid = assume_diploid,
      genes = genes,
      seg_file_source = seg_file_source,
      this_seq_type = this_seq_type)$maf
  }else if(!missing(in_seg)){
    CN_new = assign_cn_to_ssm(
      this_sample_id = this_sample_id,
      maf_file = in_maf,
      maf_df = maf_df,
      seg_file = in_seg,
      seg_file_source = seg_file_source,
      coding_only = coding_only,
      genes = genes,
      this_seq_type = this_seq_type)$maf
  }else{
    # If no seg file was provided, assume_diploid parameter is automatically set to true
    if(missing(in_seg)){
      CN_new = assign_cn_to_ssm(
        this_sample_id = this_sample_id,
        maf_file = in_maf,
        maf_df = maf_df,
        assume_diploid = TRUE,
        coding_only = coding_only,
        genes = genes,
        this_seq_type = this_seq_type)$maf
    }
  }
  # Change any homozygous deletions (CN = 0) to 1 for calculation purposes
  CN_new$CN[CN_new$CN<1] = 1
  # Select only the relevant columns, add new columns for VAF, Ploidy, and Purity
  CN_new = CN_new %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, t_ref_count, t_alt_count, CN) %>%
    dplyr::mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
    dplyr::mutate(Ploidy = ifelse(CN %in% 1:2, 1, 0)) %>%
    dplyr::mutate(Purity = "") %>%
    tidyr::drop_na(CN) %>%
    tidyr::unite("Chrom_pos", Chromosome:Start_Position, sep = ":", remove = FALSE)

  # Create an empty list
  indiv_CN = list()

  # For each copy number state, separate mutations in each CN state into separate datatables
    ## "for" loop, and this function to get different copy number states: unique(CN_new$CN)
  # Duplicate rows based on the copy number value (CN of 2 = 2 rows, etc):
    ##CN_maf[rep(seq(nrow(CN_maf)), CN_maf$CN),]
  # Fill out the ploidy column for each duplicated column:
    ##rep(seq(i),nrow(CN_max))
  # Add a column for Purity and Calculate for CNs 3 or larger:
    ##mutate(CN_max_dup, Purity = (CN*VAF)/Ploidy)

  for(i in unique(CN_new$CN)){
    CN_max = CN_new %>%
      dplyr::filter(CN == i)
    if(i > 2){
      CN_max_dup = CN_max[rep(seq(nrow(CN_max)), CN_max$CN),]
      CN_max_dup$Ploidy = rep(seq(i),nrow(CN_max))
      CN_max_dup = dplyr::mutate(CN_max_dup, Purity = (CN*VAF)/Ploidy)
      indiv_CN[[as.character(i)]] = CN_max_dup
    }else(indiv_CN[[as.character(i)]] = CN_max)
  }

  # Merge all copy state tables together into one table
  merged_CN = do.call("rbind", indiv_CN)

  # Filter for copy nnumber states 1 or 2
  # Calculate purity, and if the number is larger 1 (100%) use the VAF calue instead
  merged_CN_neut = merged_CN %>%
    dplyr::filter(merged_CN$CN < 3) %>%
    dplyr::mutate(Purity = (CN*VAF)/Ploidy) #%>%
    #dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity))

    # Calculate a temporary purity based on the mean of these purity values
    merged_CN_neut <- merged_CN_neut %>%
      drop_na(Purity)

    mean_neut_purity = mean(merged_CN_neut$Purity)

  # For CN of 3 or larger:
    ## Filter for copy number states 1 or 2
    ## Calculate purity, and if the number is larger 1 (100%) use the VAF calue instead
    ## Group by chromosonal position,
    ## For each group, choose only the purity value that is closest to the mean_neut_purity (the temporary purity calculated from the copy neutral values)
  if (is.na(mean_neut_purity)){
    merged_CN_gain = merged_CN %>%
      dplyr::filter(merged_CN$CN > 2) %>%
      dplyr::mutate(Purity = (VAF*2)) %>%
      #dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity)) %>%
      group_by(Chrom_pos) %>%
      slice_head() %>%
      dplyr::mutate(Assumed_CN = 2) %>%
      dplyr::mutate(Assumed_ploidy = 1)
  }else{
    merged_CN_gain = merged_CN %>%
      dplyr::filter(merged_CN$CN > 2) %>%
      dplyr::mutate(Purity = (CN*VAF)/Ploidy) %>%
      #dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity)) %>%
      group_by(Chrom_pos) %>%
      slice(which.min(abs(Purity - mean_neut_purity)))
  }

  # Bind both datatables back together (the first contains CNs 1 and/or 2, the second contains CNs 3 and/or higher)
  CN_final = bind_rows(merged_CN_neut, merged_CN_gain)

  # Estimate the mean of all purity values from all available copy number states
  sample_purity_estimation = mean(CN_final$Purity)

  # If the final sample purity is above 1, make it equal to 1
  if(sample_purity_estimation > 1){
    sample_purity_estimation = 1
  }

  # Calculate CCF (cancer cell fraction) 2 ways:
   ## With the maximum purity estimation for the mutations in the same (largest value will be 1 for the mutation with the highest purity estimation)
   ## With the purity estimation of the sample in total. This will give values over 1 though which is problematic, as the maximum value should be 1 (100%) for clonal events
  CN_final = CN_final %>%
    dplyr::mutate(CCF_mutation = Purity/max(Purity)) %>%
    dplyr::mutate(CCF_sample = Purity/sample_purity_estimation)

  output = list()
  if(show_plots){
    # Figure 1 : VAF distribution
    # Creates facet wraps showing the VAF distribution of CN-annotated mutations for each available copy number state
    VAF_plot = CN_final %>%
      ggplot(aes(x = VAF)) +
      geom_histogram() +
      facet_wrap(~CN)

    # Figure 2: Final purity distribution
    Purity_plot = CN_final %>%
      ggplot(aes(x = Purity)) +
      geom_histogram() +
      facet_wrap(~CN)

     output[["VAF_plot"]] = VAF_plot
     output[["Purity_plot"]] = Purity_plot
  }
  output[["sample_purity_estimation"]] = sample_purity_estimation
  output[["CN_final"]] = CN_final
  return(output)
}
