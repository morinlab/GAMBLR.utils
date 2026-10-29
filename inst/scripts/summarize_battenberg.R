#!/usr/bin/env Rscript
# Summarize Battenberg results with GAMBLR.utils::summarize_battenberg() and
# write two tables:
#   <out_prefix>.summary.tsv       one row per sample
#   <out_prefix>.loh_segments.tsv  one row per merged LOH / homozygous deletion event
#
# Usage:
#   Rscript summarize_battenberg.R <out_prefix>
#       all GAMBL tumour genomes with a matched normal, located with
#       GAMBLR.results::get_battenberg() (uses the GAMBLR config)
#   Rscript summarize_battenberg.R <out_prefix> <battenberg_dir> [genome_build]
#       every *_subclones.txt under a directory (sample = file prefix)
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("usage: Rscript summarize_battenberg.R <out_prefix> [battenberg_dir [genome_build]]")
}
out_prefix = args[1]
if(length(args) >= 2){
  genome_build = if(length(args) >= 3) args[3] else "grch37"
  bb = GAMBLR.utils::summarize_battenberg(battenberg_dir = args[2],
                                          genome_build = genome_build,
                                          verbose = TRUE)
}else{
  meta = GAMBLR.results::get_gambl_metadata(seq_type_filter = "genome")
  files = GAMBLR.results::get_battenberg(meta)
  write.table(files, paste0(out_prefix, ".files.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  bb = GAMBLR.utils::summarize_battenberg(files, verbose = TRUE)
}
write.table(bb$summary, paste0(out_prefix, ".summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(bb$loh_segments, paste0(out_prefix, ".loh_segments.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("summarized ", nrow(bb$summary), " samples")
