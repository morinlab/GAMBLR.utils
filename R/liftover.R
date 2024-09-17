#' @title Liftover.
#'
#' @description Use liftOver to convert a coordinate-based file between the two main genome builds (grch37/hg38).
#'
#' @details The user can specify the data to be converted to another genome build in a data frame with `data_df`.
#' For the bedpe data, the user can optionally also specify a path to the bedpe file that needs to be lifted with `bedpe_file`
#' (remained here for backwards compatibility).
#' 
#' The other required parameter is `target_build`, this parameter decides the final projection of the lifted data.
#' 
#' The data type can be specified with the `mode` argument, which accepts one of the bedpe (default, maf, or bed files).
#' Other formats are not supported yet.
#'
#' @param data_df Input data to be lifted to another projection in the format of data frame.
#' @param mode Specify the type of data to be converted. Available options are "bedpe" (default), "maf", and "bed".
#' @param target_build Specify which build the data should be lifted to (must be one of hg19, grch37, hg38, grch38).
#' @param bedpe_file Optionally, specify the path to a bedpe file if it is not already loaded to the environment. Only used in the bedpe mode.
#' @param col_names If not provided, the column names will be imposed. Only used in the bedpe mode.
#' @param col_types Specify column types if column names are also defined with `col_names`. Only used in the bedpe mode.
#' @param standard_bed Boolean parameter for defining the type of bed file that is provided with `bedpe_file`. Deafult is FALSE. Only used in the bedpe mode.
#' @param verbose Set to TRUE for verbose output. Default is FALSE.
#'
#' @return Data frame containing original bedpe data with new coordinates.
#'
#' @rawNamespace import(S4Vectors, except = c("merge", "second", "first", "union", "intersect", "setdiff", "setequal", "rename", "expand", "head", "tail", "end", "start"))
#' @rawNamespace import(GenomicRanges, except = c("start", "end", "merge", "shift", "union", "intersect", "setdiff", "reduce", "trim"))
#' @rawNamespace import(rtracklayer, except = c("start", "end"))
#' @import dplyr tidyr readr
#' @export
#'
#' @examples
#' bed <- grch37_ashm_regions
#'
#' test <- liftover(
#'     data_df = bed,
#'     mode = "bed",
#'     target_build = "hg38"
#' )
#'
#' maf <- get_coding_ssm(
#'      these_sample_ids = "DOHH-2"
#' )
#'
#' test <- liftover(
#'     data_df = maf,
#'     mode = "maf",
#'     target_build = "hg38"
#' )
#'
#' bedpe <- get_manta_sv(
#'      these_sample_ids = "DOHH-2"
#' )
#'
#' test <- liftover(
#'     data_df = bedpe,
#'     mode = "bedpe",
#'     target_build = "hg38"
#' )
#'

liftover = function(
        data_df,
        mode = "bedpe",
        target_build = "grch37",
        bedpe_file,
        col_names,
        col_types,
        standard_bed = FALSE,
        verbose = FALSE
    ){

    if(target_build == "grch37" | target_build == "hg19") {
        chain <- rtracklayer::import.chain(
            system.file(
                "extdata",
                "hg38ToHg19.over.chain",
                package = "GAMBLR.data"
            )
        )
    } else if(target_build == "grch38" | target_build == "hg38") {
        chain <- rtracklayer::import.chain(
            system.file(
                "extdata",
                "hg19ToHg38.over.chain",
                package = "GAMBLR.data"
            )
        )
    }

    # Original implementation of this function made to work with bedpe
    if(mode == "bedpe") {
        if(!missing(bedpe_file)) {
            if(missing(col_names)) {
                message("imposing column names")
                original_bedpe <- suppressMessages(
                    read_tsv(
                        bedpe_file,
                        comment = "##",
                        col_types = "cddcddccccccccccccccccc"
                    )
                )
            } else {
                message(
                    paste(
                        "using column names",
                        col_names,
                        sep = ": "
                    )
                )
                original_bedpe <- suppressMessages(
                    read_tsv(
                        bedpe_file,
                        col_names = col_names,
                        col_types = col_types
                    )
                )
            }
        } else {
            original_bedpe = data_df
        }


        if(!standard_bed){
            colnames(original_bedpe)[1] = "CHROM_A"
            original_bedpe = as.data.frame(original_bedpe)
            original_bedpe = original_bedpe %>%
                dplyr::mutate(
                    CHROM_A = ifelse(
                        !grepl("chr", CHROM_A),
                        paste0("chr", CHROM_A),
                        CHROM_A
                    ),
                    CHROM_B = ifelse(
                        !grepl("chr", CHROM_B),
                        paste0("chr", CHROM_B),
                        CHROM_B
                    )
                )

            original_bedpe = original_bedpe %>%
                dplyr::mutate(
                    START_A = format(START_A, scientific = FALSE),
                    START_B = format(START_B, scientific = FALSE),
                    END_A = format(END_A, scientific = FALSE),
                    END_B = format(END_B, scientific = FALSE)
                )

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
                original_bedpe = mutate(
                    original_bedpe,
                    chrom = paste0("chr", chrom)
                )
            }
            char_vec = original_bedpe %>%
                tidyr::unite(united, sep = "\t") %>%
                dplyr::pull(united)

            bedpe_obj = rtracklayer::import(text = char_vec, format = "bed")
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
                dplyr::select(
                    -c(
                        "CHROM_A", "START_A", "END_A",
                        "CHROM_B", "START_B", "END_B",
                        "STRAND_A", "STRAND_B"
                    )
                )

            fully_lifted = bind_cols(first_ok_df, second_ok_df, kept_cols) %>%
                dplyr::select(all_of(original_columns)) %>%
                dplyr::mutate_if(is.factor, as.character)

            return(fully_lifted)

        } else {
            lifted = rtracklayer::liftOver(bedpe_obj, chain)
            no_problem = !((elementNROWS(lifted) != 1))
            first_ok = subset(lifted, no_problem)

            output = rtracklayer::export(first_ok, format = "bed") %>%
                suppressMessages(
                    read_tsv(
                        col_names = c(
                            "chrom", "start", "end", "score", "strand",
                            "nothing", "s1", "e1", "junk", "more",
                            "stuff", "nada"
                        )
                    )
                ) %>%
                dplyr::select("chrom", "start", "end")

            return(output)
        }
    } else if(mode %in% c("maf", "bed")) {
        # handle different styles of column names (chrom, chr, chr_hg19, etc.)
        data_df <- data_df %>%
            rename_at(
                vars(matches("chr")), ~ "Chromosome"
            ) %>%
            rename_at(
                vars(matches("start")), ~ "Start_Position"
            ) %>%
            rename_at(
                vars(matches("end")), ~ "End_Position"
            )

        # Always make sure liftover input is chr-prefixed and avoid scenarios
        # of creating erroneous chrchr1 instances
        data_df$Chromosome <- gsub("chr", "", data_df$Chromosome)
        data_df$Chromosome <- paste0("chr", data_df$Chromosome)

        # liftover wants certain specific column names
        over <- data_df %>%
            dplyr::rename(
                "chromosome" = "Chromosome",
                "start" = "Start_Position",
                "end" = "End_Position"
            )
        over <- GenomicRanges::makeGRangesFromDataFrame(
            over,
            keep.extra.columns = TRUE
        )
        lifted <- rtracklayer::liftOver(over, chain)
        lifted <- data.frame(
            lifted@unlistData
        ) %>%
            dplyr::select(-width) %>%
            dplyr::rename(
                "Chromosome" = "seqnames",
                "Start_Position" = "start",
                "End_Position" = "end",
                "Strand" = "strand"
            ) %>%
            rename_at(
                vars(matches("strand")), ~ "Strand"
            ) %>%
            dplyr::select(all_of(colnames(data_df))) %>%
            as.data.frame()

        # If this is maf we also need to reset the NCBI_Build column value
        # to keep things sane
        if(mode == "maf" & target_build %in% c("grch37", "hg19")) {
            lifted$NCBI_Build <- "GRCh37"
        } else if(mode == "maf" & target_build %in% c("grch38", "hg38")) {
            lifted$NCBI_Build <- "GRCh38"
        }

        return(lifted)
    } else {
        stop(
            "The specified mode is not supported yet."
        )
    }
}
