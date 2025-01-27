#' Copy Number (CN) Core Processor
#' 
#' @description 
#' The primary purpose is to extract statistics from `pretty_CN_heatmap` with customized modifications, such as filtering and normalization of copy number (CN) state data and calculating Copy Number Alteration Ratios (CNARs).
#' 
#' @param cn_state_matrix A matrix of copy number (CN) states exported from get_cn_states function, where rows represent samples and columns represent genomic bins (segments).
#' @param remove_chrom A character vector of chromosomes to exclude from analysis (e.g., "chrX"). Default is `"chrX"`.
#' @param round_first Logical. If `TRUE`, CN values are rounded to integers before further analysis. Default is `FALSE`.
#' @param drop_bin_if_sd_below Numeric. Bins (segments) with a standard deviation below this threshold are excluded. Default is `0`.
#' @param nrm_by_avg_ploidy Logical. If `TRUE`, CN values are normalized by the average ploidy across all samples. Default is `FALSE`.
#' @param simplyfy_nrm Logical. If `TRUE`, normalized CN values are rounded to integers for simplicity. Default is `TRUE`.
#' @param auto_calc_cutoff Logical. If `TRUE`, cutoff values for CN state classification are calculated automatically using the method specified in `auto_calc_cutoff_method`. Default is `FALSE`.
#' @param auto_calc_cutoff_method Character. Specifies the method for automatic cutoff calculation. Options are `"median"`, `"1st_quartile"`, `"3rd_quartile"`, or `"mean"`. Default is `"median"`.
#' @param cutoff_gain_custom Numeric. Custom cutoff value for classifying gains in CN states, used when `auto_calc_cutoff` is `FALSE`. Default is `0.05`.
#' @param cutoff_loss_custom Numeric. Custom cutoff value for classifying losses in CN states, used when `auto_calc_cutoff` is `FALSE`. Default is `0.005`.
#' 
#' @return 
#' An object of S4 class `all` containing:
#' \itemize{
#'   \item \code{mat}: A `mat` object with filtered CN state matrix and optionally normalized matrix.
#'   \item \code{data}: A `data` object containing filtered bin data, CNAR statistics for bins and lengths.
#' }
#' 
#' @import dplyr tidyr
#' @export 
#'
#' @examples
#' metadata = get_gambl_metadata() %>% filter(pathology == "DLBCL")
#' cn_state_matrix = get_cn_states(n_bins_split=5000,
#'                                 missing_data_as_diploid = TRUE, # if FALSE, missing data will be NA
#'                                 seg_data = all_segments,
#'                                 these_samples_metadata = metadata)
#' result = cn_states(cn_state_matrix, 
#'                    metadata, 
#'                    remove_chrom = c("chrX"), 
#'                    auto_calc_cutoff=FALSE,
#'                    cutoff_gain_custom = 0.05,
#'                    cutoff_loss_custom = 0.005)

######################################################################################
###################                 MAIN FUNCTION                  ###################
######################################################################################
cn_states = function(cn_state_matrix,
                     remove_chrom = c("chrX"),
                     round_first=FALSE,
                     drop_bin_if_sd_below=0,
                     nrm_by_avg_ploidy=FALSE,
                     simplyfy_nrm=TRUE,
                     auto_calc_cutoff=FALSE,
                     auto_calc_cutoff_method="median",
                     cutoff_gain_custom = 0.05,
                     cutoff_loss_custom = 0.005) {

    # Check segments
    bin_data = data.frame(segments = colnames(cn_state_matrix)) %>%
        tidyr::separate(segments, into = c("chrom", "range"), sep = ":", remove = FALSE) %>%
        tidyr::separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
        mutate(genome_length = end-start+1) %>%
        select(segments, chrom, start, end, genome_length)

    # (Optional) Remove specific chromosomes
    if(!missing(remove_chrom)) {
        bin_data = bin_data %>% filter(! chrom%in%remove_chrom )
        nchrom_removed = ncol(cn_state_matrix) - nrow(bin_data)
        bin_data_to_retain = bin_data %>% pull(segments)
        cn_state_matrix = cn_state_matrix[colnames(cn_state_matrix) %in% bin_data_to_retain]
        cat("\n segments from", remove_chrom, paste0("(n =", nchrom_removed, ")"), " were removed\n")
    }

    # (Optional) Round CN values
    if (round_first==TRUE) { cn_state_matrix = round(cn_state_matrix) }

    # Filtering segmented regions (=bin_data) by SD
    bin_cn_sd = apply(cn_state_matrix, 2, sd)
    bin_data = cbind(bin_data, bin_cn_sd)
    nbin_pre = nrow(bin_data)
    bin_data = bin_data[bin_cn_sd > drop_bin_if_sd_below,]
    nbin_post = nrow(bin_data)
    cn_state_matrix = cn_state_matrix[colnames(cn_state_matrix) %in% bin_data$segments]
    cn_state_matrix_flt = cn_state_matrix
    if(nbin_pre - nbin_post > 0) {
        cat("\n", nbin_pre - nbin_post, "segment(s) was(were) excluded due to zero variability (SD=0). \n")
        } else { cat("\n---- All segmetns have variability (SD > 0). ----\n") }


    # Let's define useful parameters
    params = cn_state_matrix %>%
        list(nbin_data = ncol(.),
            mask_diploid = . == 2,
            mask_nondiploid = . != 2,
            mask_gain = . > 2,
            mask_loss = . < 2)

    # (Optional) Normalize CN values by average ploidy for the comparison among samples
    avg_ploidy = apply(cn_state_matrix, 1, mean)
    if(nrm_by_avg_ploidy==TRUE) {
        cn_state_matrix = cn_state_matrix - avg_ploidy + 2 
        if(simplyfy_nrm==TRUE) { cn_state_matrix = round(cn_state_matrix_nrm) }
        cn_state_matrix_flt_nrm = cn_state_matrix
        cat("\n CN values were normalized by average ploidy.\n")
    } else {
        cn_state_matrix_flt_nrm = NULL
    }

    # Calculate Copy Number Alteration Ratio (CNAR)
    cnar_bins = calc_cnar(cn_state_matrix, bin_data, params, mode="bins", auto_calc_cutoff, cutoff_gain_custom, cutoff_loss_custom)
    cnar_length = calc_cnar(cn_state_matrix, bin_data, params, mode="length", auto_calc_cutoff, cutoff_gain_custom, cutoff_loss_custom)

    # Returns
    result = return_as_S4(
        cn_state_matrix_flt, 
        cn_state_matrix_flt_nrm, 
        bin_data, 
        cnar_bins, 
        cnar_length, 
        nrm_by_avg_ploidy
    )

    return(result)
}


######################################################################################
###################               HELPER FUNCTIONS                 ###################
######################################################################################
calc_cutoff = function(values, method) {
    if (method == "median") {
        return(median(values))
    } else if (method == "1st_quartile") {
            return(quantile(values, 0.25))
    } else if (method == "3rd_quartile") {
        return(quantile(values, 0.75))
    } else if (method == "mean") {
        return(mean(values))
    } else if (method == "custom") {
        return(mean(values))
    } else {
        stop("Invalid method. Use 'median', '1st_quartile', '3rd_quartile', or 'mean'.")
    }
}

calc_cnar = function(cn_state_matrix, 
                     bin_data, params, 
                     mode, 
                     auto_calc_cutoff,
                     cutoff_gain_custom, 
                     cutoff_loss_custom) {
    mode = match.arg(mode, choices = c("bins", "length"))
    cat("\nCalculating Copy Number Alteration Ratio (CNAR) with mode =",mode,"\n")

    if (mode=="length") {
        genome_length_checked = bin_data$genome_length[match(colnames(cn_state_matrix), bin_data$segments)]
        cn_state_matrix = sweep(cn_state_matrix, 2, genome_length_checked, FUN = "*")
    }

    stats = cn_state_matrix %>%
        reframe(!!sym(paste0("diploid_", mode)) := rowSums(. * params$mask_diploid),
                !!sym(paste0("gain_", mode)) := rowSums(. * params$mask_gain),
                !!sym(paste0("loss_", mode)) := rowSums(.* params$mask_loss),
                !!sym(paste0("nondiploid_", mode)) := rowSums(. * params$mask_nondiploid),
                !!sym(paste0("total_", mode)) := rowSums(.)) %>%
        mutate(!!sym(paste0("diploid_ratio_", mode)) := !!sym(paste0("diploid_", mode)) / !!sym(paste0("total_", mode)),
               !!sym(paste0("gain_ratio_", mode)) := !!sym(paste0("gain_", mode)) / !!sym(paste0("total_", mode)),
               !!sym(paste0("loss_ratio_", mode)) := !!sym(paste0("loss_", mode)) / !!sym(paste0("total_", mode)),
               !!sym(paste0("nondiploid_ratio_", mode)) := !!sym(paste0("nondiploid_", mode)) / !!sym(paste0("total_", mode)))
    rownames(stats) = rownames(cn_state_matrix)

    # Let's define Copy Number states
    if (auto_calc_cutoff==TRUE) {
        cutoff_gain = stats %>% pull(!!sym(paste0("gain_ratio_", mode))) %>% calc_cutoff(., method=auto_calc_cutoff_method)
        cutoff_loss = stats %>% pull(!!sym(paste0("loss_ratio_", mode))) %>% calc_cutoff(., method=auto_calc_cutoff_method)
        cat("\nCutoff values were calculated using the", auto_calc_cutoff_method, "values\n")
        cat("\ncutoff_gain:", cutoff_gain, ", cutoff_loss:", cutoff_loss, "\n")
    } else {
        cutoff_gain = cutoff_gain_custom
        cutoff_loss = cutoff_loss_custom
        cat("\nCustom cutoff values were used\n")
        cat("\ncutoff_gain:", cutoff_gain, ", cutoff_loss:", cutoff_loss, "\n")
    }

    stats = stats %>%
        mutate(!!sym(paste0("cn_state_", mode)) := case_when(
            !!sym(paste0("gain_ratio_", mode)) >= cutoff_gain & !!sym(paste0("loss_ratio_", mode)) >= cutoff_loss ~ "gain_and_loss",
            !!sym(paste0("gain_ratio_", mode)) >= cutoff_gain & !!sym(paste0("loss_ratio_", mode)) < cutoff_loss ~ "gain",
            !!sym(paste0("gain_ratio_", mode)) < cutoff_gain & !!sym(paste0("loss_ratio_", mode)) >= cutoff_loss ~ "loss",
            !!sym(paste0("gain_ratio_", mode)) < cutoff_gain & !!sym(paste0("loss_ratio_", mode)) < cutoff_loss ~ "non-altered",
            TRUE ~ NA))

    cat("\nCopy Number States\n")
    print( stats %>% count(!!sym(paste0("cn_state_", mode))) )

    return(stats)
}

# Define S4 Class
repr_mat = representation(flt = "data.frame", flt_nrm = "ANY")
setClass("mat", representation = repr_mat)

repr_data = representation(bin_data = "data.frame", cnar_bins = "data.frame", cnar_length = "data.frame")
setClass("data", representation = repr_data)

repr_all = representation(mat = "mat", data = "data")
setClass("all", representation = repr_all)

return_as_S4 = function(cn_state_matrix_flt, 
                        cn_state_matrix_flt_nrm,  
                        bin_data, 
                        cnar_bins, 
                        cnar_length,
                        nrm_by_avg_ploidy=FALSE) {
    mat = new("mat",
        flt = cn_state_matrix_flt,
        flt_nrm = if (nrm_by_avg_ploidy==TRUE) cn_state_matrix_flt_nrm else NULL
    )

    data = new("data",
        bin_data = bin_data,
        cnar_bins = cnar_bins,
        cnar_length = cnar_length
    )

    new("all", mat = mat, data = data)
}
