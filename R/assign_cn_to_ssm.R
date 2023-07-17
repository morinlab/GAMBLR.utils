#' @title Assign CN to SSM.
#'
#' @description Annotate mutations with their copy number information.
#'
#' @details This function takes a sample ID with the `this_sample_id` parameter and annotates mutations with copy number information.
#' A variety of parameters are at hand for a customized workflow. For example, the user can specify if only coding mutations are of interest.
#' To do so, set `coding_only = TRUE`. It is also possible to point the function to already loaded maf/seq files, or a path to these files.
#' See parameters; `maf_file`, `maf_path`, `seq_file` and `seg_path` for more information on how to use these parameters.
#' This function can also take a vector with genes of interest (`genes`) that the returned data frame will be restricted to.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR::get_cn_segments], [GAMBLR::get_cn_states], [GAMBLR::get_sample_cn_segments]
#'
#' @param this_sample_id Sample ID of the sample you want to annotate.
#' @param coding_only Optional. set to TRUE to restrict to only coding variants.
#' @param from_flatfile Optional. Instead of the database, load the data from a local MAF and seg file.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is FALSE.
#' @param tool_name name of tool to be used, default is "battenberg".
#' @param maf_file Path to maf file.
#' @param maf_df Optional. Use a maf dataframe instead of a path.
#' @param seg_file path to seq file.
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, Battenberg, etc.
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral.
#' @param genes Genes of interest.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE
#' @param this_seq_type Specified seq type for returned data.
#' @param projection specified genome projection that returned data is in reference to.
#'
#' @return A list containing a data frame (MAF-like format) with two extra columns:
#' log.ratio is the log ratio from the seg file (NA when no overlap was found)
#' as well as the segmented copy number data with the same copy number information
#' CN is the rounded absolute copy number estimate of the region based on log.ratio (NA when no overlap was found)
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr readr RMariaDB DBI ssh
#' @export
#'
#' @examples
#' cn_list = assign_cn_to_ssm(this_sample_id = "HTMCP-01-06-00422-01A-01D",
#'                            coding_only = TRUE)
#'
assign_cn_to_ssm = function(this_sample_id,
                            coding_only = FALSE,
                            from_flatfile = TRUE,
                            use_augmented_maf = TRUE,
                            tool_name = "battenberg",
                            maf_file,
                            maf_df,
                            seg_file,
                            seg_file_source = "battenberg",
                            assume_diploid = FALSE,
                            genes,
                            include_silent = FALSE,
                            this_seq_type = "genome",
                            projection = "grch37"){

  seq_type = this_seq_type

  remote_session = check_remote_configuration(auto_connect = TRUE)
  database_name = check_config_value(config::get("database_name"))
  project_base = check_config_value(config::get("project_base"))
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  if(!missing(maf_file)){
    maf_sample = fread_maf(maf_file) %>%
      dplyr::mutate(Chromosome = gsub("chr", "", Chromosome))
  }
  else if(!missing(maf_df)){
    maf_sample = maf_df %>%
      dplyr::mutate(Chromosome = gsub("chr", "", Chromosome))
  }
  else if(from_flatfile){
    #get the genome_build and other wildcards for this sample
    wildcards = get_sample_wildcards(this_sample_id,seq_type)
    genome_build = wildcards$genome_build
    unix_group = wildcards$unix_group
    seq_type = wildcards$seq_type
    tumour_sample_id = wildcards$tumour_sample_id
    normal_sample_id = wildcards$normal_sample_id
    pairing_status = wildcards$pairing_status
    maf_sample = get_ssm_by_sample(this_sample_id = this_sample_id, this_seq_type = this_seq_type, augmented = use_augmented_maf)

  }else{
    #get all the segments for a sample and filter the small ones then assign CN value from the segment to all SSMs in that region
    con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    maf_table = check_config_value(config::get("results_tables")$ssm)
    maf_sample = dplyr::tbl(con, maf_table) %>%
      dplyr::filter(Tumor_Sample_Barcode == this_sample_id) %>%
      as.data.frame()
  }
  if(coding_only){
    maf_sample = dplyr::filter(maf_sample, Variant_Classification %in% coding_class)
  }

  if(!missing(genes)){
    maf_sample = dplyr::filter(maf_sample, Hugo_Symbol %in% genes)
  }

  if(!missing(seg_file)){
    seg_sample = suppressMessages(read_tsv(seg_file)) %>%
      dplyr::mutate(size = end - start) %>%
      dplyr::filter(size > 100)

    colnames(seg_sample)[c(1:4)] = c("ID", "chrom", "start", "end")
    seg_sample = seg_sample %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
      dplyr::rename(Chromosome = chrom, Start_Position = start, End_Position = end) %>%
      data.table::as.data.table()

    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
  }else if(assume_diploid == TRUE){
    if(missing(seg_file)){
      print("WARNING: A seg file was not provided! Annotating all mutation calls as copy neutral")
    }

    a = data.table::as.data.table(maf_sample)
    a_diploid = dplyr::mutate(a, CN = 2)
    return(list(maf = a_diploid))

  }else if(from_flatfile){
    message(paste("trying to find output from:", tool_name))
    project_base = check_config_value(config::get("project_base",config="default"))
    local_project_base = check_config_value(config::get("project_base"))

    results_path_template = check_config_value(config::get("results_flatfiles")$cnv$battenberg)
    results_path = paste0(project_base, results_path_template)
    local_results_path = paste0(local_project_base, results_path_template)

    ## NEED TO FIX THIS TO contain tumour/normal ID from metadata and pairing status
    battenberg_file = glue::glue(results_path)
    local_battenberg_file = glue::glue(local_results_path)


    message(paste("looking for flatfile:", battenberg_file))
    if(remote_session){
      print(local_battenberg_file)
      dirN = dirname(local_battenberg_file)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_battenberg_file)){

        ssh::scp_download(ssh_session,battenberg_file,dirN)
      }
      battenberg_file = local_battenberg_file
    }

    #check for missingness
    if(!file.exists(battenberg_file)){
      print(paste("missing: ", battenberg_file))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    seg_sample = suppressMessages(read_tsv(battenberg_file)) %>%
      as.data.table() %>%
      dplyr::mutate(size = end - start) %>%
      dplyr::filter(size > 100) %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
      dplyr::rename(Chromosome = chrom, Start_Position = start, End_Position = end)

    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
  }else{
    seg_sample = get_sample_cn_segments(this_sample_id = this_sample_id, this_seq_type = this_seq_type) %>%
      dplyr::mutate(size = end - start) %>%
      dplyr::filter(size > 100) %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
      dplyr::rename(Chromosome = chrom, Start_Position = start, End_Position = end) %>%
      data.table::as.data.table()

    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
    a.seg = data.table::foverlaps(a, seg_sample, type = "any")
    a$log.ratio = a.seg$log.ratio
    a$LOH = factor(a.seg$LOH_flag)
    a = dplyr::mutate(a, CN = round(2*2^log.ratio))
    seg_sample = dplyr::mutate(seg_sample, CN = round(2*2^log.ratio))
    seg_sample$LOH_flag = factor(seg_sample$LOH_flag)
  }
  if(seg_file_source == "ichorCNA"){
      message("defaulting to ichorCNA format")
      seg_sample = dplyr::rename(seg_sample, c("log.ratio" = "median", "CN" = "copy.number"))
      a.seg = data.table::foverlaps(a, seg_sample, type = "any")
      a$log.ratio = a.seg$log.ratio
      a$LOH = factor(a.seg$LOH_flag)
      a$CN = a.seg$CN

    }else{
      a.seg = data.table::foverlaps(a, seg_sample, type = "any")
      a$log.ratio = a.seg$log.ratio
      a$LOH = factor(a.seg$LOH_flag)
      a = dplyr::mutate(a, CN = round(2*2^log.ratio))
      seg_sample = dplyr::mutate(seg_sample, CN = round(2*2^log.ratio))
      seg_sample$LOH_flag = factor(seg_sample$LOH_flag)
  }
  if(!missing(seg_sample)){
    return(list(maf = a, seg = seg_sample))
  }
}
