#' @title Write Sample Set Hash
#'
#' @description Update or create a file to track unique identifiers for sample sets in GAMBL
#'
#' @details Run this function with `update = TRUE` (default) to use an existing sample table.
#' If this table does not exist, perhaps you need to pull from the master branch.
#' If this function is run with the default for `update`, the user must also provide the new sample sets with the `new_sample_sets_df`.
#'
#' @param update Leave as TRUE for default functionality (i.e. updating the existing table). If the table doesn't exist you probably need to pull from Master.
#' @param new_sample_sets_df Data frame of all existing and new sample sets. Required when running in default update mode.
#'
#' @import dplyr readr
#' @export
#'
write_sample_set_hash = function(update = TRUE,
                                 new_sample_sets_df){

  sample_sets_file = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("sample_sets")$default))
  md5_file = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("sample_sets")$hashes))

  if(update){
    # load the existing file and update it using the contents of sample_sets_df as well as checking for consistency for existing sample sets
    if(missing(new_sample_sets_df)){
      stop("You must provide a data frame containing all the sample sets to update the digests")
    }
    original_digests = suppressMessages(read_tsv(md5_file))

    set_names = dplyr::select(new_sample_sets_df,-sample_id) %>% colnames()
    #only compare for sample sets that we have in the current file
    md5_values= c()
    for(set_name in set_names){

      this_md5 = get_samples_md5_hash(sample_set_name=set_name,sample_sets_df = new_sample_sets_df)
      md5_values=c(md5_values,this_md5)

    }
    all_md5 = data.frame(sample_set = set_names,new_md5_digest=md5_values)
    oldnew = right_join(original_digests,all_md5,by="sample_set")
    #check the rows where md5_digest is not NA (i.e. rows that were there before)
    to_check = dplyr::filter(oldnew,!is.na(md5_digest))
    if(any(to_check$md5_digest != to_check$new_md5_digest)){
      problems = dplyr::filter(to_check,md5_digest!=new_md5_digest)
      print(problems)
      stop("some md5 digests do not match. Have these sample sets changed???")

    }

  }else{
    # just create a file that records the md5 digests for existing sample sets
    sample_sets = suppressMessages(read_tsv(sample_sets_file))
    set_names = select(sample_sets,-sample_id) %>% colnames()
    message(paste("Will log digests for",length(set_names),"sample sets"))
    md5_values= c()
    for(set_name in set_names){
      this_md5 = get_samples_md5_hash(sample_set_name=set_name)
      md5_values=c(md5_values,this_md5)
    }
    all_md5 = data.frame(sample_set = set_names,md5_digest=md5_values)
    write_tsv(all_md5,file=md5_file)
  }
}
