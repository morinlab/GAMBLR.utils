#' @title Fishers Exact Test (CNV).
#'
#' @description Using GISTIC2.0 outputs, perform Fisher's exact test to compare
#' CNV frequencies between 2 groups.
#'
#' @details This function was developed to compare (Fisher's exact test) CNV
#' frequencies between two groups. To do so, set the path to the GISTIC2.0
#' all_lesions file with `gistic_lesions`, together with a metadata table with
#' the sample IDs of interest (`metadata`). The last remaining required
#' parameter is `comparison`, this parameter takes the name of the column
#' annotating the groups of interest, e.g pathology, cohort, etc.
#' For more information on how to run this function, refer to the function
#' examples and parameter descriptions.
#'
#' @param gistic_lesions Path to the GISTIC2.0 all_lesions output file.
#' @param metadata Data frame containing sample ids and column with annotated
#'      data for the 2 groups of interest. All other columns are ignored.
#'      The function will exit with error if asked to compare more than two
#'      groups, so please ensure to provide correct metadata.
#' @param comparison Specify column annotating groups of interest.
#' @param fdr.method FDR method to adjust p values. Uses p.adjust function, and
#'      therefore accepts its method for FDR ("holm", "hochberg", "hommel",
#'      "bonferroni", "BH", "BY", "fdr", "none"). By default, this function
#'      uses "fdr".
#' @param fdr.cutoff Specify FDR significance cut-off. By default, this
#'      function uses 0.1.
#' @param confidence_level Numeric value between 0 and 1 representing the
#'      confidence level for the plotted confidence intervals.
#' @param text_size Size of the text on the forest plot of differentially
#'      enriched CNV. Default text-size is 7.
#' @param blacklisted_regions Optionally, specify any descriptors (value from
#'      column `Descriptor` of GISTIC2.0 all_lesions output file) to filter out
#'      before any comparisons are done. It is possible to specify a list of
#'      multiple descriptors, for example, c("3p12.3", "12p13.2"). Default
#'      is NULL.
#' @param point_size Number specifying size of the red point on the resulting
#'      forestplot.
#'
#' @return list
#'
#' @importFrom stats qnorm
#' @import dplyr ggpubr readr tidyr ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' # basic usage
#' FtestCNV(gistic_lesions = "path_to_GISTIC2.0_output/all_lesions.conf_{confidence_level}.txt",
#'          metadata = derived_data,
#'          comparison = "pathology")
#'
#' # advanced usage
#' FtestCNV(gistic_lesions = "path_to_GISTIC2.0_output/all_lesions.conf_{confidence_level}.txt",
#'          metadata = derived_data,
#'          comparison = "pathology",
#'          fdr.method = "bonferroni",
#'          fdr.cutoff = 0.05,
#'          blacklisted_regions = c("3p12.3", "12p13.2"))
#' }
#'
FtestCNV = function(
    gistic_lesions,
    metadata,
    comparison,
    fdr.method = "fdr",
    fdr.cutoff = 0.1,
    confidence_level = 0.95,
    text_size = 7,
    blacklisted_regions = NULL,
    point_size = 3
){

    # get groups of interest for comparison
    GROUPS.TO.COMPARE = unique(metadata[,comparison]) %>% unlist

    # quick check that only 2 groups are provided forr comparison
    if(length(GROUPS.TO.COMPARE) > 2){
        message(
            paste0(
                "The current implementation of function only accepts 2 groups ",
                "for comparison. You provided 3 groups."
            )
        )
        message("Please modify metadata accordingly to compare only 2 groups.")
        message(
            "Groups you provided are: ", paste(
                c(GROUPS.TO.COMPARE),
                collapse = ","
            )
        )
        return(NULL)
    }
    # read lesions from gistic utput to collect event/case
    lesions = suppressMessages(read_tsv(gistic_lesions, col_names = TRUE)) %>%
        dplyr::filter(!grepl("CN", `Unique Name`)) %>%
        dplyr::select(
            -tail(names(.), 1),
            -`Residual q values after removing segments shared with higher peaks`,
            -`Broad or Focal`,
            -`Amplitude Threshold`,
            -`q values`
        ) %>%
        dplyr::filter (! Descriptor %in% blacklisted_regions)

    # prepare metadata
    # save names of first few columns of gistic lesions file
    columns = colnames(lesions)[1:5]
    # subset lesions file to only lesions metadata and samples of interest
    columns = c(columns, pull(metadata, Tumor_Sample_Barcode))
    lesions = subset(lesions, select = intersect(columns, colnames(lesions)))

    # subset ids for each category for future comparisons
    GROUP1 = dplyr::pull(
        metadata %>%
            dplyr::filter(
                base::get(comparison) == GROUPS.TO.COMPARE[1]
            ),
        Tumor_Sample_Barcode
    )

    GROUP2 = dplyr::pull(
        metadata %>%
            dplyr::filter(
                base::get(comparison) == GROUPS.TO.COMPARE[2]
            ),
        Tumor_Sample_Barcode
    )

    # subset lesions file for each group and count number of CNV
    GROUP1.CNV = lesions[,colnames(lesions) %in% GROUP1] %>%
        dplyr::mutate(count = rowSums(.!=0))

    GROUP2.CNV = lesions[,colnames(lesions) %in% GROUP2] %>%
        dplyr::mutate(count = rowSums(.!=0))

    ########## PART 1
    # compare significance of differences between 2 groups

    # first, get the counts of mutated samples for each group
    GROUP1_vs_GROUP2 = as.data.frame(
        cbind(lesions[,1:5],
        GROUP1.CNV$count,
        GROUP2.CNV$count)
    )
    colnames(GROUP1_vs_GROUP2)[6:7] = c("Mutated_GROUP1", "Mutated_GROUP2")

    # second, get the counts of unmutated samples for each group
    GROUP1_vs_GROUP2 = GROUP1_vs_GROUP2 %>%
        dplyr::mutate(
            Not_Mutated_GROUP1 = length(GROUP1) - Mutated_GROUP1,
            .after = Mutated_GROUP1
        ) %>%
        dplyr::mutate(
            Not_Mutated_GROUP2 = length(GROUP2) - Mutated_GROUP2,
            .after = Mutated_GROUP2
        ) %>%
        dplyr::mutate(
            Region = ifelse(
                grepl("Amplification", `Unique Name`),
                paste0(Descriptor, "_Amp"),
                paste0(Descriptor, "_Del")
            ),
            .after = Descriptor
        )

    # row-by-row, calculate 2-tailed Fisher's exact test for each CNV
    GROUP1_vs_GROUP2 = GROUP1_vs_GROUP2 %>%
        dplyr::filter(
            (Mutated_GROUP1 / (Mutated_GROUP1 + Not_Mutated_GROUP1) > 0.05) |
            (Mutated_GROUP2 / (Mutated_GROUP2 + Not_Mutated_GROUP2) > 0.05)
        ) %>%
        dplyr::mutate(
            pVal = apply(
                .,
                1,
                function(x) {
                    tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T);
                    fisher.test(tbl, alternative = "two.sided")$p.value
                }
            )
        ) %>%
        dplyr::mutate(
            OddsRatio = apply(
                .,
                1,
                function(x) {
                    tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T);
                    log(fisher.test(tbl, alternative = "two.sided")$estimate)
                }
            )
        ) %>%
        dplyr::mutate(
            LowConfInt = apply(
                .,
                1,
                function(x) {
                    tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T);
                    log(fisher.test(tbl, alternative = "two.sided")$conf.int[1])
                }
            )
        ) %>%
        dplyr::mutate(
            HighConfInt = apply(
                .,
                1,
                function(x) {
                    tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T);
                    log(fisher.test(tbl, alternative = "two.sided")$conf.int[2])
                }
            )
        )

    # adjust FDR for multiple comparisons
    GROUP1_vs_GROUP2$FDR = p.adjust(GROUP1_vs_GROUP2$pVal, method = fdr.method)

    # filter for only those CNV that pass FDR cutoff
    GROUP1_vs_GROUP2.PASSED = GROUP1_vs_GROUP2[GROUP1_vs_GROUP2$FDR<fdr.cutoff,]
    GROUP1_vs_GROUP2.PASSED = GROUP1_vs_GROUP2.PASSED %>%
        dplyr::distinct(Region, .keep_all = TRUE)

    # rename columns to match with comparison groups and save resulting
    # df as part of the output
    DISTINCT = GROUP1_vs_GROUP2.PASSED %>%
        `colnames<-`(gsub("GROUP1", GROUPS.TO.COMPARE[1], colnames(.))) %>%
        `colnames<-`(gsub("GROUP2", GROUPS.TO.COMPARE[2], colnames(.)))

    message("Successfully completed step 1/3...")

    ########## PART 2
    # prepare data with CNV for more convenient downstream analyses, like survival

    # First, collect events for the group 1
    # this just extracts metadata for each CNV from lesions file
    REGIONS = as.data.frame(lesions[,1:5]) %>%
        dplyr::mutate(
            Region = ifelse(
                grepl("Amplification", `Unique Name`),
                paste0(Descriptor, "_Amp"),
                paste0(Descriptor, "_Del")
            ),
            .after = Descriptor
        )

    # Collect events for samples in group 1. Since same band can contain more
    # than 1 peak, keep peak limits for future unique differentiation
    survival_object = 0
    # here, i is a CNV feature, anf j is a sample_id
    for(i in 1:length(REGIONS$Region)) {
        for(j in 1:1:ncol(GROUP1.CNV[,-ncol(GROUP1.CNV)])) {
            row = c(
                colnames(GROUP1.CNV[j]),
                REGIONS$Region[i],
                REGIONS$`Wide Peak Limits`[i],
                GROUP1.CNV[i,j]
            )
            survival_object = rbind(survival_object, row)
        }
    }
    # tidy output for convenience
    survival_object = as.data.frame(survival_object)
    colnames(survival_object) = c(
        "sample_id", "Region", "Wide Peak Limits", "Event"
    )
    survival_object = survival_object[2:length(survival_object$sample),]
    survival_object = survival_object %>%
        dplyr::filter(!Event==0) %>%
        mutate(Event = ifelse(
            grepl("Amp", Region),
            "GAIN",
            "LOSS"
            )
        )

    # Create final object that will be concatenated with events for second group further down
    ALL.EVENTS = survival_object

    # Collect events for samples in group 2
    survival_object = 0

    # here, i is a CNV feature, anf j is a sample_id
    for(i in 1:length(REGIONS$Region)) {
        for(j in 1:ncol(GROUP2.CNV[,-ncol(GROUP2.CNV)])) {
            row = c(
                colnames(GROUP2.CNV[j]),
                REGIONS$Region[i],
                REGIONS$`Wide Peak Limits`[i],
                GROUP2.CNV[i,j]
            )
            survival_object = rbind(survival_object, row)
        }
    }
    # tidy output for convenience
    survival_object = as.data.frame(survival_object)
    colnames(survival_object) = c(
        "sample_id", "Region", "Wide Peak Limits", "Event"
    )
    survival_object = survival_object[2:length(survival_object$sample),]
    survival_object = survival_object %>% dplyr::filter(!Event==0) %>%
        mutate(Event = ifelse(
            grepl("Amp", Region),
            "GAIN",
            "LOSS"
            )
        )

    # Now, merge with the object prepared earlier and save it
    # as part of the output
    CNV.EVENTS = tidyr::unnest(rbind(ALL.EVENTS, survival_object))

    message("Successfully completed step 2/3...")

    ########## PART 3
    # visualize forest plot in a convenient way

    # calculate SE
    mergedPassed = GROUP1_vs_GROUP2.PASSED %>%
        dplyr::mutate(SE = (HighConfInt - LowConfInt) / 5.95) %>%
        dplyr::mutate(interval = qnorm(1 - (1 - confidence_level)/2)*SE)

    # order in decreasing order for better visualization
    mergedPassed = mergedPassed[
        order(mergedPassed$OddsRatio, decreasing = TRUE),
    ]

    # prepare table that will be plotted
    study_table = data.frame(
        name = mergedPassed[, "Region"],
        Events_GROUP1 = paste(
            mergedPassed$Mutated_GROUP1,
            "/",
            mergedPassed$Mutated_GROUP1+mergedPassed$Not_Mutated_GROUP1,
            sep = ""
        ),
        Events_GROUP2 = paste(
            mergedPassed$Mutated_GROUP2,
            "/",
            mergedPassed$Mutated_GROUP2+mergedPassed$Not_Mutated_GROUP2,
            sep = ""
        )
    )

    # rename columns to match with comparison groups
    study_table = study_table %>%
        `colnames<-`(gsub("GROUP1", GROUPS.TO.COMPARE[1], colnames(.))) %>%
        `colnames<-`(gsub("GROUP2", GROUPS.TO.COMPARE[2], colnames(.)))

    # actual plot
    res <- bind_cols(
        mergedPassed,
        study_table
    ) %>%
    mutate(Region = factor(Region, levels = rev(Region)))


    p_mid <- res %>%
        ggplot(aes(y = Region)) +
        geom_linerange(
            aes(
                xmin=OddsRatio-interval,
                xmax=OddsRatio+interval,
                linewidth = (Mutated_GROUP1 + Mutated_GROUP2)
            )
        ) +
        geom_point(
            aes(x = OddsRatio),
            shape = 15,
            size = point_size,
            color = "red"
        ) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(
            x = "Log(OddsRatio)",
            y = ""
        ) +
        theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(
                size = text_size,
                face = "bold",
                color = "black"
            ),
            axis.title = element_text(
                size = text_size,
                face = "bold"
            ),
            panel.background = element_blank(),
            panel.grid.major.x = element_line(color = "grey"),
            panel.border = element_rect(colour = "black", fill = NA)
        )


    p_left <- res %>%
        mutate(name = Region) %>%
        select(Region, name, starts_with("Events")) %>%
        pivot_longer(
            !Region,
            names_to = "column",
            values_to = "label"
        ) %>%
        mutate(column = factor(column, levels = unname(c("name", paste0("Events_", GROUPS.TO.COMPARE)))
    )) %>%
        ggplot(aes(y = Region)) +
        geom_text(
            aes(
                x = column,
                label = label,
                size = text_size
            )
        ) +
        scale_x_discrete(position = "top") +
        theme(
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(
                size = text_size,
                face = "bold",
                color = "black"
            )
        )


    p_right <- res %>%
        select(Region, OddsRatio, interval, FDR) %>%
        mutate(
            label = paste0(
                round(OddsRatio, 2),
                " [",
                round(OddsRatio - interval, 2),
                ", ",
                round(OddsRatio + interval, 2),
                "]"
            ),
            column = "Log(OddsRatio)"
        ) %>%
        ggplot(aes(y = Region)) +
        geom_text(
            aes(
                x = column,
                label = label,
                size = text_size
            )
        ) +
        scale_x_discrete(position = "top") +
        theme(
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(
                size = text_size,
                face = "bold",
                color = "black"
            )
        )

    GRAPH <- ggarrange(
        p_left,
        p_mid,
        p_right,
        widths = c(1,1,0.3),
        ncol = 3,
        align = "h"
    )

    message("Successfully completed step 3/3...")

    OUTPUT = list(
        "DISTINCT" = DISTINCT,
        "CNV.EVENTS" = CNV.EVENTS,
        "GRAPH" = GRAPH
    )

    return(OUTPUT)
    message("Done!")
}
