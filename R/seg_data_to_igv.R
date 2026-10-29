#' @title Write segmented data as an IGV seg file
#'
#' @description Write copy number segments in the .seg format that IGV
#' loads directly (tab-separated: ID, chrom, loc.start, loc.end, num.mark
#' and the value shown), e.g. after centring them with [center_segments].
#'
#' @details IGV colours segments by the last column, on a log2 scale (blue
#' below 0, red above), so `value` chooses what is shown:
#' \itemize{
#'   \item "log2_CN" (the default when there is a CN column): the copy
#'     number as `log2(CN / 2)`, so 2 copies show as neutral, 1 as -1, 4 as
#'     +1 (0 copies is floored at -4). After [center_segments] this is the
#'     centred copy number.
#'   \item "log.ratio" (the default otherwise): the caller's log ratio, a
#'     separate measurement relative to the sample's own average and
#'     compressed by tumour purity, not a transform of CN.
#'   \item "CN", or any other numeric column, as it is (for CN, IGV's colour
#'     scale then needs setting around 2).
#' }
#' Placeholder rows (`dummy_segment == 1`) are left out.
#'
#' With `these_samples_metadata`, a sample information file is also
#' written (next to the .seg file, as `<name>.sample_info.txt`, unless
#' `metadata_file` is given). Load it into IGV after the .seg file (File >
#' Load from File) to show the metadata as coloured attribute columns
#' beside the samples, and to sort, group or filter samples by them. Only
#' samples in both are written: segments of samples not in the metadata are
#' dropped (with a message), and metadata rows without segments are left
#' out. The ID column matches the .seg IDs.
#' Values that have a GAMBL colour (`GAMBLR.helpers::get_gambl_colours`,
#' e.g. pathology, seq_type) get it, in the file's `#colors` section.
#'
#' @param seg_data Segments with ID, chrom, start, end and the `value`
#' column (e.g. from `GAMBLR.results::get_cn_segments`).
#' @param file Path of the .seg file to write. When NULL, nothing is
#' written and the table is only returned.
#' @param value What IGV shows: "log2_CN" (default when there is a CN
#' column), "log.ratio" (default otherwise), "CN" or another numeric
#' column. See Details.
#' @param chr_prefix NULL (default) leaves chromosome names as they are;
#' TRUE adds a "chr" prefix and FALSE removes it, to match the IGV genome.
#' @param these_samples_metadata Optional metadata (with sample_id): only
#' its samples' segments are written, and it is written as an IGV sample
#' information file. See Details.
#' @param metadata_columns Metadata columns to write. Default sample_id,
#' seq_type, cohort, capture_space and pathology (those missing from the
#' metadata are skipped).
#' @param additional_columns Further columns to write after
#' `metadata_columns`, e.g. c("lymphgen", "COO_consensus"). Columns that
#' aren't in the metadata are taken from `seg_data` when it has them, one
#' value per sample (e.g. "module"; a sample with more than one value gets
#' them joined by "+").
#' @param metadata_file Path of the sample information file. Default: the
#' .seg file name with `.sample_info.txt` in place of `.seg`.
#' @param colours Add GAMBL colours for metadata values that have one.
#' Default TRUE.
#'
#' @return The seg table, invisibly; with `these_samples_metadata`, a list
#' with `seg` and `sample_info` (invisibly).
#'
#' @export
#'
#' @examples
#' # one segment per autosome, spanning each whole chromosome (a baseline
#' # needs segments covering at least half of the autosomal genome)
#' arms = GAMBLR.data::chromosome_arms_grch37
#' chrom_end = tapply(arms$end, arms$chromosome, max)[as.character(1:22)]
#' seg = data.frame(ID = "s1", chrom = as.character(1:22), start = 1,
#'                  end = as.numeric(chrom_end), CN = c(3, rep(4, 21)))
#' seg_data_to_igv(center_segments(seg), file.path(tempdir(), "centred.seg"))
#'
#' \dontrun{
#' # with a sample information file (centred.sample_info.txt) for IGV
#' seg_data_to_igv(center_segments(seg), "centred.seg",
#'                 these_samples_metadata = meta,
#'                 additional_columns = c("lymphgen", "pairing_status"))
#' }
seg_data_to_igv = function(seg_data, file = NULL, value = NULL, chr_prefix = NULL,
                           these_samples_metadata = NULL,
                           metadata_columns = c("sample_id", "seq_type", "cohort",
                                                "capture_space", "pathology"),
                           additional_columns = NULL,
                           metadata_file = NULL,
                           colours = TRUE){
  s = as.data.frame(seg_data)
  if(!all(c("ID", "chrom", "start", "end") %in% colnames(s))){
    stop("seg_data needs columns ID, chrom, start and end")
  }
  if(is.null(value)) value = if("CN" %in% colnames(s)) "log2_CN" else "log.ratio"
  if(value == "log2_CN"){
    # copy number (e.g. centred by center_segments) on IGV's log2 scale:
    # 2 copies is 0; 0 copies is floored at -4
    if(!"CN" %in% colnames(s)) stop("value = \"log2_CN\" needs a CN column")
    s$log2_CN = log2(pmax(as.numeric(s$CN), 2 * 2^-4) / 2)
  }else if(!value %in% colnames(s)){
    stop("seg_data has no ", value, " column")
  }
  if("dummy_segment" %in% colnames(s)){
    s = s[!s$dummy_segment %in% 1, , drop = FALSE]
  }
  # with metadata, only its samples are written
  if(!is.null(these_samples_metadata)){
    if(!"sample_id" %in% colnames(these_samples_metadata)){
      stop("these_samples_metadata needs a sample_id column")
    }
    keep = as.character(s$ID) %in% as.character(these_samples_metadata$sample_id)
    if(any(!keep)){
      message("dropping ", length(unique(s$ID[!keep])), " sample(s) not in these_samples_metadata")
    }
    s = s[keep, , drop = FALSE]
  }
  chrom = gsub("^chr", "", as.character(s$chrom))
  if(is.null(chr_prefix)){
    chrom = as.character(s$chrom)
  }else if(chr_prefix){
    chrom = paste0("chr", chrom)
  }
  out = data.frame(ID = as.character(s$ID), chrom = chrom,
                   loc.start = format(s$start, scientific = FALSE, trim = TRUE),
                   loc.end = format(s$end, scientific = FALSE, trim = TRUE),
                   num.mark = if("num.mark" %in% colnames(s)) s$num.mark else NA,
                   value = as.numeric(s[[value]]), stringsAsFactors = FALSE)
  colnames(out)[6] = if(value == "log.ratio") "seg.mean" else value
  out = out[order(out$ID, out$chrom, as.numeric(out$loc.start)), , drop = FALSE]
  rownames(out) = NULL
  if(!is.null(file)){
    utils::write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  if(is.null(these_samples_metadata)) return(invisible(out))
  info = .igv_sample_info(s, these_samples_metadata, metadata_columns, additional_columns)
  if(is.null(metadata_file) && !is.null(file)){
    metadata_file = paste0(sub("\\.seg$", "", file), ".sample_info.txt")
  }
  if(!is.null(metadata_file)){
    .write_igv_sample_info(info, metadata_file, colours)
  }
  invisible(list(seg = out, sample_info = info))
}

# One row of metadata per sample in the segments, in the requested columns
# (sample_id first, as the ID IGV matches to the .seg file)
.igv_sample_info = function(seg, these_samples_metadata, metadata_columns, additional_columns){
  meta = as.data.frame(these_samples_metadata)
  requested = setdiff(unique(c(metadata_columns, additional_columns)), "sample_id")
  # columns the metadata lacks can come from the segments (e.g. module),
  # one value per sample
  from_seg = setdiff(intersect(requested, colnames(seg)), c(colnames(meta), "ID"))
  missing_extra = setdiff(additional_columns, c(colnames(meta), colnames(seg)))
  if(length(missing_extra) > 0){
    stop("neither these_samples_metadata nor seg_data has column(s): ",
         paste(missing_extra, collapse = ", "))
  }
  skipped = setdiff(metadata_columns, c(colnames(meta), colnames(seg)))
  if(length(skipped) > 0){
    message("not in the metadata, so not written: ", paste(skipped, collapse = ", "))
  }
  ids = unique(as.character(seg$ID))
  meta = meta[meta$sample_id %in% ids, , drop = FALSE]
  # a sample_id with genome and capture rows: keep the row of the profile
  # in the segments when that's known, otherwise the first
  if("seg_seq_type" %in% colnames(seg) && "seq_type" %in% colnames(meta)){
    matched = paste(meta$sample_id, meta$seq_type) %in% paste(seg$ID, seg$seg_seq_type)
    meta = meta[matched | !meta$sample_id %in% meta$sample_id[matched], , drop = FALSE]
  }
  if(any(duplicated(meta$sample_id))){
    message("keeping the first metadata row of samples listed more than once")
    meta = meta[!duplicated(meta$sample_id), , drop = FALSE]
  }
  info = meta[, intersect(c("sample_id", requested), colnames(meta)), drop = FALSE]
  for(column in from_seg){
    # the sample's value; samples with more than one get them all, joined by "+"
    values = tapply(as.character(seg[[column]]), as.character(seg$ID), function(v){
      v = sort(unique(v[!is.na(v)]))
      if(length(v) == 0) NA_character_ else paste(v, collapse = "+")
    })
    info[[column]] = as.vector(values[as.character(info$sample_id)])
    mixed = sum(grepl("+", info[[column]], fixed = TRUE))
    if(mixed > 0){
      message(mixed, " sample(s) have more than one ", column, " value in seg_data, written joined by \"+\"")
    }
  }
  info = info[, c("sample_id", intersect(requested, colnames(info))), drop = FALSE]
  colnames(info)[1] = "ID"
  rownames(info) = NULL
  info
}

# IGV sample information file: a #sampleTable section, then a #colors
# section (attribute, value, R,G,B) for values with a GAMBL colour
.write_igv_sample_info = function(info, file, colours = TRUE){
  table = info
  table[] = lapply(table, function(v){ v = as.character(v); v[is.na(v)] = ""; v })
  lines = c("#sampleTable", paste(colnames(table), collapse = "\t"),
            apply(table, 1, paste, collapse = "\t"))
  if(colours){
    palette = GAMBLR.helpers::get_gambl_colours()
    col_lines = unlist(lapply(colnames(table)[-1], function(attribute){
      values = setdiff(unique(table[[attribute]]), "")
      values = values[values %in% names(palette)]
      if(length(values) == 0) return(NULL)
      rgb = grDevices::col2rgb(palette[values])
      paste(attribute, values, paste(rgb[1, ], rgb[2, ], rgb[3, ], sep = ","), sep = "\t")
    }))
    if(length(col_lines) > 0) lines = c(lines, "", "#colors", col_lines)
  }
  writeLines(lines, file)
  invisible(file)
}
