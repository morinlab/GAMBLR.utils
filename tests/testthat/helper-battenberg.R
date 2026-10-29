# Write a synthetic Battenberg sample (subclones and, optionally, fitting
# files) to a temporary directory. `segments` is a data frame with chr,
# startpos, endpos, nMaj1_A, nMin1_A and optionally frac1_A, nMaj2_A, nMin2_A,
# frac2_A; missing columns get clonal defaults.
write_battenberg_sample = function(dir, sample, segments,
                                   rho_and_psi = NULL, cellularity = NULL,
                                   malformed_rho = FALSE, refit_rows = NULL){
  seg = segments
  if(is.null(seg$frac1_A)) seg$frac1_A = 1
  for(col in c("nMaj2_A", "nMin2_A", "frac2_A")){
    if(is.null(seg[[col]])) seg[[col]] = NA
  }
  tot = seg$nMaj1_A + seg$nMin1_A
  seg$BAF = ifelse(tot > 0, seg$nMaj1_A / pmax(tot, 1), 0.5)
  seg$pval = 1
  seg$LogR = log2(pmax(tot, 0.1) / 2)
  seg$ntot = tot
  seg = seg[, c("chr", "startpos", "endpos", "BAF", "pval", "LogR", "ntot",
                "nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]
  utils::write.table(seg, file.path(dir, paste0(sample, "_subclones.txt")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  rho_file = file.path(dir, paste0(sample, "_rho_and_psi.txt"))
  if(malformed_rho){
    writeLines(c("not a battenberg table", "1 2"), rho_file)
  }else if(!is.null(rho_and_psi)){
    # Battenberg writes the row names (ASCAT, FRAC_GENOME, REF_SEG) unlabelled
    utils::write.table(rho_and_psi, rho_file, sep = "\t", quote = FALSE,
                       col.names = NA)
    lines = readLines(rho_file); lines[1] = sub("^(\"\")?\t", "", lines[1])
    writeLines(lines, rho_file)
  }
  if(!is.null(cellularity)){
    utils::write.table(cellularity, file.path(dir, paste0(sample, "_cellularity_ploidy.txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if(!is.null(refit_rows)){
    utils::write.table(refit_rows, file.path(dir, paste0(sample, "_refit_suggestion.txt")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  file.path(dir, paste0(sample, "_subclones.txt"))
}

# One segment per autosome (5-50 Mb), all in the given state
genome_segments = function(major = 1, minor = 1){
  data.frame(chr = as.character(1:22), startpos = 5e6, endpos = 5e7,
             nMaj1_A = major, nMin1_A = minor)
}

fit_table = function(ascat_ploidy = 2.27, clonal_ploidy = 2.76, clonal_distance = 0.6,
                     ref_seg_ploidy = 2.91, ref_seg_distance = 0.59,
                     best = c("FRAC_GENOME", "REF_SEG")[1]){
  data.frame(rho = c(1, 0.95, 1), psi = c(1.55, 1.85, 1.98),
             ploidy = c(ascat_ploidy, clonal_ploidy, ref_seg_ploidy),
             distance = c(NA, clonal_distance, ref_seg_distance),
             is.best = c(NA, best == "FRAC_GENOME", best == "REF_SEG"),
             row.names = c("ASCAT", "FRAC_GENOME", "REF_SEG"))
}

summarize_one = function(dir, ...){
  summarize_battenberg(battenberg_dir = dir, loh_genes = NULL, ...)
}

# A diploid genome with one segment per chromosome arm (grch37 arm
# coordinates; acrocentric p arms left uncalled, as Battenberg leaves them)
arm_genome = function(major = 1, minor = 1){
  arms = as.data.frame(GAMBLR.data::chromosome_arms_grch37)
  arms$chromosome = gsub("^chr", "", as.character(arms$chromosome))
  arms = arms[arms$chromosome %in% as.character(1:22) &
                !(arms$chromosome %in% c("13", "14", "15", "21", "22") & arms$arm == "p"), ]
  data.frame(chr = arms$chromosome, startpos = arms$start, endpos = arms$end,
             nMaj1_A = major, nMin1_A = minor, frac1_A = 1,
             nMaj2_A = NA, nMin2_A = NA, frac2_A = NA)
}

# Set chrom:start-end to a state, splitting the segments it overlaps
set_region = function(seg, chrom, start, end, major, minor, ...){
  out = list()
  for(i in seq_len(nrow(seg))){
    r = seg[i, ]
    if(r$chr != chrom || r$endpos <= start || r$startpos >= end){
      out[[length(out) + 1]] = r; next
    }
    if(r$startpos < start){ left = r; left$endpos = start; out[[length(out) + 1]] = left }
    mid = r; mid$startpos = max(r$startpos, start); mid$endpos = min(r$endpos, end)
    mid$nMaj1_A = major; mid$nMin1_A = minor
    extra = list(...)
    for(col in names(extra)) mid[[col]] = extra[[col]]
    out[[length(out) + 1]] = mid
    if(r$endpos > end){ right = r; right$startpos = end; out[[length(out) + 1]] = right }
  }
  out = do.call(rbind, out)
  out[order(as.integer(out$chr), out$startpos), ]
}
