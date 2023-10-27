#' @title SV To Custom Track.
#'
#' @description Make a UCSC-ready custom track file from SV data.
#'
#' @details This function takes an incoming SV data frame and outputs a bed file, ready for visualization on the UCSC Genome Browser.
#' Specify the output file with `output_file`, indicate if the incoming SVs are annotated with `is_annotated` (default is TRUE).
#' Lastly, the user can also specify if the incoming SV data frame should be subset to specific mutation types (e.g deletions, duplications, insertions, etc.).
#' This is specified with the `sv_name` parameter. Default is to include all SV subtypes.
#'
#' @param sv_bedpe A bedpe formatted data frame of SVs.
#' @param output_file A bed file with UCSC custom header.
#' @param is_annotated Set to TRUE if input SV bedpe is annotated, default is TRUE.
#' @param sv_name SV name. Default is set to "all" = include all subtypes of SVs.
#'
#' @return Nothing.
#'
#' @import dplyr tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' #custom track with annotations
#' all_sv = GAMBLR.data::sample_data$grch37$bedpe
#' annotated_sv = annotate_sv(sv_data = all_sv)
#' sv_to_custom_track(annotated_sv,
#'                    output_file = "GAMBL_sv_custom_track_annotated.bed",
#'                    is_annotated = TRUE)
#'
#' #custom track (no anotatated SVs)
#' sv_to_custom_track(all_sv,
#'                    output_file = "GAMBL_sv_custom_track_annotated.bed",
#'                    is_annotated = FALSE)
#' }
#'
sv_to_custom_track = function(sv_bedpe,
                              output_file,
                              is_annotated = TRUE,
                              sv_name = "all"){

  if(is_annotated){
    #reduce to a bed-like format
    sv_data1 = mutate(annotated_sv, annotation = paste0(chrom1, ":", start1, "_", fusion)) %>%
      dplyr::select(chrom2, start2, end2, tumour_sample_id, annotation, fusion)

    sv_data2 = mutate(annotated_sv, annotation = paste0(chrom2, ":", start2, "_", fusion)) %>%
      dplyr::select(chrom1, start1, end1, tumour_sample_id, annotation, fusion)

    print(head(sv_data1))
    print(head(sv_data2))
    colnames(sv_data1) = c("chrom", "start", "end", "sample_id", "annotation", "fusion")
    colnames(sv_data2) = c("chrom", "start", "end", "sample_id", "annotation", "fusion")
    sv_data = bind_rows(sv_data1, sv_data2)
    sv_data = mutate(sv_data, end = end + 10)
  }else{
    sv_data_1 = mutate(sv_bedpe, annotation = paste0( CHROM_B, ":", START_B)) %>%
      dplyr::select(CHROM_A, START_A, END_A, tumour_sample_id, annotation)

    sv_data_2 = mutate(sv_bedpe, annotation = paste0( CHROM_A, ":", START_A)) %>%
      dplyr::select(CHROM_B, START_B, END_B, tumour_sample_id, annotation)

    colnames(sv_data_1) = c("chrom", "start", "end", "sample_id", "annotation")
    colnames(sv_data_2) = c("chrom", "start", "end", "sample_id", "annotation")
    sv_data = bind_rows(sv_data_1, sv_data_2)

  }
  if(!any(grepl("chr", sv_data[,1]))){
    #add chr
    sv_data[,1] = unlist(lapply(sv_data[,1], function(x){paste0("chr", x)}))
  }

  coo_cols = GAMBLR.helpers::get_gambl_colours("COO")
  path_cols = GAMBLR.helpers::get_gambl_colours("pathology")
  all_cols = c(coo_cols, path_cols)
  colour_df = data.frame(coo = names(all_cols), colour = all_cols)

  rgb_df = data.frame(t(col2rgb(all_cols))) %>%
    mutate(consensus_coo_dhitsig = names(all_cols)) %>%
    unite(col = "rgb", red, green, blue, sep = ",")

  meta = GAMBLR.helpers::handle_metadata() %>%
    dplyr::select(sample_id, "consensus_coo_dhitsig", pathology) %>%
    mutate(consensus_coo_dhitsig = if_else(consensus_coo_dhitsig == "NA", pathology, consensus_coo_dhitsig))

  samples_coloured = left_join(meta, rgb_df)
  sv_bed_coloured = left_join(sv_data, samples_coloured) %>%
    arrange(pathology)

  write_bed = function(coloured_svs, sv_name, output_file_base){
    data_bed = coloured_svs %>%
      mutate(details = paste0(sample_id, "_", annotation)) %>%
      mutate(score = 0, strand = "+", end = end + 1, start1 = start, end1 = end) %>%
      dplyr::select(chrom, start, end, details, score, strand, start1, end1, rgb) %>%
      dplyr::filter(!is.na(rgb)) %>%
      unique()

    header_content = paste0('track name="GAMBL SVs ', sv_name, '" description="SV breakpoints ', sv_name, '" visibility=2 itemRgb="On"\n')
    cat(header_content, file = output_file)
    message(paste("writing to", output_file))
    tabular = write.table(data_bed, file = output_file, quote = F, sep = "\t", row.names = F, col.names = F, append = TRUE)
  }
  write_bed(sv_bed_coloured, sv_name = sv_name)
}
