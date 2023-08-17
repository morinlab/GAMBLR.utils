#' @title Annotate IGH Breakpoints.
#'
#' @description Add annotation to IGH breakpoints to infer mechanism based on location within IGH.
#'
#' @details Returns a modified bedpe with additional columns.
#' The most useful columns that are added are mechanism (one of CSR, AID, VDJ) and label (NA if unmatched, otherwise one of Emu, Smu, one of the J segments or switch regions).
#'
#' @param annotated_df Previously annotated data frame of SVs.
#' @param genome_build Version of reference build to be used, only grch37 currently accepted.
#'
#' @return A slightly modified bedpe with added columns.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr
#' @export
#'
#' @examples
#' manta_sv = get_manta_sv(verbose = FALSE)
#' all_annotated = annotate_sv(sv_data = manta_sv)
#' ig_annotated = annotate_igh_breakpoints(all_annotated)
#'
annotate_igh_breakpoints = function(annotated_df,
                                    genome_build = "grch37"){

  if(genome_build != "grch37"){
    message("Currently, only grch37 is supported")
    return()
  }
  sv_igh = dplyr::filter(annotated_df, partner == "IGH")

  switch_end = 106327893 #hg19 coordinate
  vdj_start = 106330635
  emu_start= 106327853
  emu_end = 106330900
  j_segs = data.frame(start = c(106329408, 106330010, 106330412, 106330789, emu_start, 106323800, 106054077, 106320301, 106304755, 106231780, 106206750, 106172700, 106107740, 106066513, 106089185, 106134348),
                      end = c(106329487, 106330175, 106330579, 106330991, emu_end, 106327600, 106058458, 106322397, 106312100, 106237884, 106209489, 106174976, 106112209, 106068295, 106092889, 106136399),
                      label = c("J6", "J5", "J4", "J3", "Emu", "Smu", "A2", "M", "D", "G3", "G1", "A1", "G2", "E", "G4", "GP"),
                      y = c(-0.5, -0.5, -0.5, -0.5, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))

  annoying_segments = c("J6", "J3P", "J5", "J4", "J3", "Smu", "M")
  j_segs = mutate(j_segs, start_wide = ifelse(!label %in% annoying_segments, start - 3500, start), end_wide = ifelse(!label %in% annoying_segments, end + 5000, end))
  sv_igh_anno = mutate(sv_igh, mechanism = case_when(chrom2 %in% c("14", "chr14") & start2 < switch_end ~ "CSR",
                                                     chrom2 %in% c("14", "chr14") & start2 < vdj_start ~ "AID",
                                                     chrom2 %in% c("14", "chr14") & start2 >= vdj_start ~ "VDJ",
                                                     chrom1 %in% c("14", "chr14") & start1 < switch_end ~ "CSR",
                                                     chrom1 %in% c("14", "chr14") & start1 < vdj_start ~ "AID",
                                                     chrom1 %in% c("14", "chr14") & start1 >= vdj_start ~ "VDJ",
                                                     TRUE ~ "OTHER"))
  sv_igh_anno = sv_igh_anno %>%
    mutate(IGH_start = ifelse(chrom2 %in% c("14", "chr14"), start2, start1)) %>%
    mutate(IGH_end = IGH_start + 1)

  #assign a higher resolution position by matching to nearest annotation
  sv_igh_anno_dt = as.data.table(sv_igh_anno)

  j_segs_dt = dplyr::filter(j_segs, label != "Emu") %>%
    as.data.table()

  #annotate Emu separately because it overlaps with some annoying J segments
  setkey(sv_igh_anno_dt, IGH_start, IGH_end)
  setkey(j_segs_dt, start_wide, end_wide)
  foverlaps(sv_igh_anno_dt, j_segs_dt, by.x = c("IGH_start", "IGH_end"), by.y = c("start_wide", "end_wide"))
  annotated = foverlaps(sv_igh_anno_dt,j_segs_dt, by.x = c("IGH_start", "IGH_end"), by.y = c("start_wide", "end_wide"))
  annotated[is.na(label) & IGH_start > emu_start & IGH_end < emu_end, "label"] = "Emu"
  return(annotated)
}
