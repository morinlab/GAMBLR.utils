#' @title Annotate intrachromosomal SV.
#'
#' @description Annotate the intrachromosomal events in a data frame of SV breakpoints.
#'
#' @details Specify a data frame with SVs (preferably the output from get_manta_sv() to the `sv_data` parameter and get back the same data frame with SV annotations.
#'
#' @param sv_data A data frame of SVs. This should be the output of get_manta_sv.
#' The required columns are columns for both chromosomes, coordinates and strand plus SOMATIC_SCORE, tumour_sample_id, and manta_name (can be NAs for output from other tools).
#' @param return_as Stated format for returned output, default is "bedpe". Other accepted output formats are "bed" and "bedpe_entrez" (to keep entrez_ids for compatibility with portal.R and cBioPortal).
#' @param projection Reference genome build parameter. Can be one of grch37 or hg38. The default is grch37.
#' @param types A vector of sv types to be annotated. The default is set to "DEL", "DUP", "INS", and "INV".
#' @param max_size An integer specifying the maximum distance between start of first breakpoint and end of second breakpoint. The default is 100000 (100000 bp).
#' @param keep_genes Optional argument allowing ot only return svs annotated to the list of genes provided in this argument. The default is "all" (return all svs annotated to all genes).
#' @param verbose Optionally control whether to print any status info to console. Useful for debugging. Default is FALSE (no messages).
#'
#'
#' @return A data frame with annotated SVs (gene symbol and entrez ID).
#'
#' @import tidyr dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' sv_df = get_manta_sv(verbose = FALSE)
#' annotated = annotate_intrachromosomal_sv(sv_data = sv_df)
#'}
annotate_intrachromosomal_sv = function(
    sv_data,
    return_as = "bedpe",
    projection = "grch37",
    types = c("DEL", "DUP", "INS", "INV", NA),
    max_size = 100000,
    keep_genes = "all",
    verbose = FALSE
){
    # Ensure input is supplied
    if(missing(sv_data)) {
        stop(
            "You did not provide input data frame with SVs to annotate."
        )
    }

    #get the genes in each region
    intra_sv <- sv_data %>%
        dplyr::mutate(
            type = gsub("Manta|:.*", "", manta_name)
        ) %>%
        dplyr::mutate(size = END_B - START_A) %>%
        dplyr::filter(
            CHROM_A == CHROM_B,
            type %in% types,
            size <= max_size
        )
    
    if(nrow(intra_sv) == 0){
        stop(
            paste(
                "There are no intrachromosomal SVs in the provided data",
                "that satisfy specified type and max_size arguments."
            )
        )
    }

    get_affected_genes = function(region, projection){
        if(verbose) {
            genes <- region_to_gene(
                region = region,
                projection = projection
            )
        } else {
            genes <- suppressMessages(
                region_to_gene(
                    region = region,
                    projection = projection
                )
            )
        }

        if("all" %in% keep_genes){
            genes <- genes %>%
            filter(hugo_symbol != "NA")
        }else{
            genes <- genes %>%
            filter(hugo_symbol %in% keep_genes)
        }
        genes <- unique(genes$hugo_symbol)

        return(paste(genes,collapse = ","))
    }

    intra_sv$genes <- NA

    for(i in c(1:nrow(intra_sv))){
        region <- paste0(
            intra_sv[i,"CHROM_A"],
            ":",
            intra_sv[i,"START_A"],
            "-",
            intra_sv[i,"END_B"]
        )
        if(verbose) {
            print(region)
        }
        g <- tryCatch(
            {
                get_affected_genes(region = region, projection = projection)
            },
                error = function(cond) {
                    if(verbose) {
                        message("Here's the original error message:")
                        message(cond)
                    }
                    # Choose a return value in case of error
                    return(NA)
                }
            )
        intra_sv[i, "genes"] <- g
        if(verbose) {
            print(g)
        }
    }
    #TODO: implement return_as feature for eventual planned cBioPortal functionality
    return(intra_sv)
}
