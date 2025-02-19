#' Prepare GISTIC input files
#'
#' @param these_samples_metadata
#' @param projection
#' @param output_dir
#' @param flavour
#' @param overwrite
#'
#' @return
#' @export
#'
#' @examples
#'
#' BL_meta = get_gambl_metadata() %>% dplyr::filter(pathology=="BL")
#' prepare_gistic_inputs(BL_meta,"hg38","/Users/rmorin/Desktop/GISTIC2/BL/","battenberg")
#'
#'
#' FL_meta = get_gambl_metadata() %>% dplyr::filter(pathology=="FL",seq_type!="mrna")
#' prepare_gistic_inputs(FL_meta,"grch37","/Users/rmorin/Desktop/GISTIC2/FL/","battenberg")
prepare_gistic_inputs <- function(these_samples_metadata,
                                  projection = "grch37",
                                  output_dir = "/Users/rmorin/Desktop/GISTIC2/",
                                  flavour = "combined",
                                  overwrite = FALSE) {
    #check that output directory exists
    if(!dir.exists(output_dir)){
      stop("you must provide an existing path")
    }
    markers <- paste0(output_dir, "markers")
    out_seg <- paste0(output_dir, "all.seg")
    if(!overwrite){
      if(file.exists(markers)|file.exists(out_seg)){
        print("output files already exist in this directory")
        stop("If you want to clobber these, re-run with overwrite = TRUE")
      }
    }
    # get the seg_data
    seg_for_gistic <- get_cn_segments(
        these_samples_metadata = these_samples_metadata,
        flavour = flavour,
        fill_missing_with = "avg_ploidy",
        projection = projection
    )
    #check how many samples we have data from
    u_id = length(unique(seg_for_gistic$ID))
    print(paste("got data from",u_id,"sample IDs"))
    seg_for_gistic <- seg_for_gistic %>%
        dplyr::select(ID:log.ratio) %>%
        mutate(LOH_flag = ifelse(LOH_flag == 99, 0, LOH_flag))
    print("checking for overlapping segments")
    seg_for_gistic <- check_overlap(seg_for_gistic)

    if ("overlap" %in% seg_for_gistic$overlap_status) {
        print("resolving overlapping segments")
        seg_for_gistic <- solve_overlap(seg_for_gistic)
    } else {
      seg_for_gistic <- seg_for_gistic %>%
            select(-overlap_status, -region_size)
    }

    # write to disk and make markers file
    write_tsv(seg_for_gistic, file = out_seg)
    temp_markers <- paste0(output_dir, "temp_markers")
    print("processing for GISTIC")
    cmd <- paste("sed '1d'", out_seg, "| cut -f2,3 > ", temp_markers)
    system(cmd)
    cmd <- paste("sed '1d'", out_seg, "| cut -f2,4 >> ", temp_markers)
    system(cmd)
    cmd <- paste("sort -V -k1,1 -k2,2nr", temp_markers, "|uniq | nl > ", markers)
    system(cmd)
    # sed '1d' {input.seg} | cut -f2,3 > {output.temp_markers}
    # sed '1d' {input.seg} | cut -f2,4 >> {output.temp_markers}
    # sort -V -k1,1 -k2,2nr {output.temp_markers} | uniq | nl > {output.markers}
    print("done!")
    print("Running gistic:")
    print("conda activate /projects/rmorin_scratch/conda_environments/30ee048690c736b4964ed302c88be45f_")
    print("gistic2 -b /projects/rmorin/projects/gambl-repos/gambl-rmorin/gistic_test")
    print("-seg ./inputs/all.seg -mk inputs/markers -refgene ref/grch37.refgene.mat -genegistic 1")
    print("-broad 1 -savegene 1 -conf 0.9 -v 30 -saveseg 0 -savedata 0")

}


check_overlap <- function(seg) {
    highest_end <- 0
    overlap <- c()
    region_sizes <- c()
    for (i in 1:nrow(seg)) {
        if (i > 1 && seg$ID[i] == seg$ID[i - 1] && seg$chrom[i] == seg$chrom[i - 1]) {
            if (seg$start[i] >= highest_end) {
                overlap[i] <- "NOToverlap"
                region_sizes[i] <- "FALSE"
            } else {
                overlap[i] <- "overlap"
                region_sizes[i] <- (seg$end[i] - seg$start[i])
            }
            if (seg$end[i] > highest_end) {
                highest_end <- seg$end[i]
            }
        } else {
            highest_end <- seg$end[i]
            overlap[i] <- "NA"
            region_sizes[i] <- (seg$end[i] - seg$start[i])
        }
    }
    seg <- seg %>% mutate(overlap_status = overlap, region_size = region_sizes)
    return(seg)
}

solve_overlap <- function(seg) {
    num_overlap <- which(seg$overlap_status == "overlap")
    num_pre_overlap_sorted <- (unique(sort(c(num_overlap - 1, num_overlap))))
    non_overlap <- seg[-num_pre_overlap_sorted, ]
    seg <- seg[num_pre_overlap_sorted, ]
    for (i in 1:nrow(seg)) {
        if (seg$overlap_status[i] == "overlap") {
            if (seg$end[i] < seg$end[i - 1]) {
                new_row1 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$start[i - 1],
                    end = seg$start[i],
                    LOH_flag = seg$LOH_flag[i - 1],
                    log.ratio = seg$log.ratio[i - 1],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                new_row2 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$end[i],
                    end = seg$end[i - 1],
                    LOH_flag = seg$LOH_flag[i - 1],
                    log.ratio = seg$log.ratio[i - 1],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                seg <- seg[-(i - 1), ]
                seg <- rbind(seg, new_row1, new_row2)
            } else if (seg$start[i - 1] == seg$start[i]) {
                new_row <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$end[i - 1],
                    end = seg$end[i],
                    LOH_flag = seg$LOH_flag[i],
                    log.ratio = seg$log.ratio[i],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                seg <- seg[-(i), ]
                seg <- rbind(seg, new_row)
            } else if (seg$region_size[i] < seg$region_size[i - 1]) {
                new_row1 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$start[i - 1],
                    end = seg$start[i],
                    LOH_flag = seg$LOH_flag[i - 1],
                    log.ratio = seg$log.ratio[i - 1],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                new_row2 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$start[i],
                    end = seg$end[i - 1],
                    LOH_flag = seg$LOH_flag[i],
                    log.ratio = seg$log.ratio[i],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                new_row3 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$end[i - 1],
                    end = seg$end[i],
                    LOH_flag = seg$LOH_flag[i],
                    log.ratio = seg$log.ratio[i],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                seg <- seg[-c(i, i - 1), ]
                seg <- rbind(seg, new_row1, new_row2, new_row3)
            } else if (seg$region_size[i] > seg$region_size[i - 1]) {
                new_row1 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$start[i - 1],
                    end = seg$start[i],
                    LOH_flag = seg$LOH_flag[i - 1],
                    log.ratio = seg$log.ratio[i - 1],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                new_row2 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$start[i],
                    end = seg$end[i - 1],
                    LOH_flag = seg$LOH_flag[i - 1],
                    log.ratio = seg$log.ratio[i - 1],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                new_row3 <- data.frame(
                    ID = seg$ID[i],
                    chrom = seg$chrom[i],
                    start = seg$end[i - 1],
                    end = seg$end[i],
                    LOH_flag = seg$LOH_flag[i],
                    log.ratio = seg$log.ratio[i],
                    overlap_status = "NOToverlap",
                    region_size = "FALSE"
                )
                seg <- seg[-c(i, i - 1), ]
                seg <- rbind(seg, new_row1, new_row2, new_row3)
            }
        }
        seg <- seg %>%
            arrange(ID, chrom, start, end)
    }
    seg <- seg %>% arrange(ID, chrom, start, end)
    seg <- check_overlap(seg)
    while ("overlap" %in% seg$overlap_status) {
        seg <- check_overlap(solve_overlap(seg))
    }
    seg <- rbind(non_overlap, seg) %>%
        arrange(ID, chrom, start, end) %>%
        select(-overlap_status, -region_size) %>%
        filter(!start == end)
    return(seg)
}
