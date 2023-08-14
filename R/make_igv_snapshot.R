#' @title Make IGV Snapshot
#'
#' @description Load bams and generate an IGV screenshot for one or more regions.
#'
#' @details Specify the path to one or more bam files as a character vector to the `bams` parameter.
#' The user can also specify regions of interest with either the `region` parameter (chr:start-end),
#' or the user can directly supply the chromosome, start and end coordinates with the `chrom`, `start`, and `end` parameters.
#' For more information and examples, refer to the function examples and parameter descriptions.
#' IMPORTANT: you must be running IGV on the host that is running R and you need to have it listening on a port.
#' The simplest scenario is to run this command on a terminal (if using a Mac),
#' assuming you are using R on gphost10 and you have a ssh config that routes gp10 to that host
#'
#' ```
#' ssh -X gp10
#' ```
#'
#' then launch IGV (e.e. from a conda installation):
#'
#' ```
#' conda activate igv; igv &
#' ```
#'
#' @param these_sample_ids A vector of one or more sample_id (bams for these samples will be auto-loaded)
#' @param this_seq_type TO DO: automatically obtain this for the user from the metadata
#' @param genome_build String specifying the genome build for the bam files provided (TO DO: if it isn't already, it should be determined automatically if these_sample_ids was provided).
#' @param bams Character vector containing the full path to one or more bam files (specify if not providing these_sample_ids)
#' @param region Optionally specify the region as a single string (e.g. "chr1:1234-1235").
#' @param padding Optionally specify a positive value to broaden the region around the specified position. Default is 200.
#' @param chrom Optionally specify the region by specifying the chromosome, start and end (see below).
#' @param start Optionally specify the region by specifying the start.
#' @param end Optionally specify the region by specifying the end.
#' @param out_path Specify the output directory where the snapshot will be written.
#' @param igv_port Specify the port IGV is listening on. Default: 60506 (optional if using the default).
#' @param socket Provide the socket variable obtained by running this function with no arguments
#' @param gene Optionally provide a gene name that will be incorporated into the output file name
#' @param details Optionally provide any other text you want incorporated into the output file name
#' @param clobber Force existing file to be clobbered?
#' @param sort_by Specify whether and how to sort the reads (e.g. "base"; see IGV documentation)
#' @param colour_by Specify how IGV should colour the reads (see IGV documentation)
#' @param squish Force reads to be squished (see IGV documentation)
#' @param viewaspairs Set to TRUE if you want the reads to be shown as pairs rather than independently (see IGV documentation)
#'
#' @return Path to file (.png).
#'
#' @import SRAdb
#' @export
#'
#' @examples
#' \dontrun{
#' this_sv = annotated_sv %>%
#'  filter(gene=="ETV6")
#'
#' #you don't need to know the details for the bam file but you can supply it if you want
#' tumour_bam = GAMBLR.results::get_bams(sample = this_sv$tumour_sample_id)
#'
#' #run with no arguments to get the socket for a running IGV instance
#' socket = make_igv_snapshot()
#'
#' make_igv_snapshot(chrom = this_sv$chrom2,
#'                   start = this_sv$start2,
#'                   end = this_sv$end2,
#'                   this_sample_id = this_sv$tumour_sample_id,
#'                   out_path = "~/IGV_snapshots/")
#'
#' this_mutation = GAMBLR.data::sample_data$grch37$maf %>%
#'  head(1)
#'
#' make_igv_snapshot(socket = socket,
#'                   sample_ids = this_mutation$Tumor_Sample_Barcode,
#'                   this_seq_type = "capture",
#'                   colour_by = "READ_STRAND")
#'
#' #run on a bunch of variants using apply:
#' apply(to_snapshot,1,function(x){
#'  make_igv_snapshot(sample_ids = x["sample_id"],
#'                    seq_type_filter = "capture",
#'                    chrom = x["chr"],
#'                    start = x["start"],
#'                    end = x["end"],
#'                    details = paste0(x["ref"],"-",x["alt"]),
#'                    gene = x["Hugo_Symbol"],
#'                    socket = socket)})
#' }
#'
make_igv_snapshot = function(bams,
                             these_sample_ids,
                             this_seq_type="genome",
                             genome_build,
                             region,
                             padding = 200,
                             chrom,
                             start,
                             end,
                             out_path = "/tmp/",
                             igv_port = 60506,
                             socket,
                             gene="NA",
                             details="",
                             clobber=FALSE,
                             sort_by="base",
                             colour_by,
                             squish=FALSE,
                             viewaspairs=FALSE){
  if(missing(socket)){
    print("returning socket for future use")
    if(exists("currently_loaded_bam")){
      currently_loaded_bam = "" #avoid the function thinking iGV still has that loaded
    }
    sock = IGVsocket(port = igv_port)
    return(sock)
  }
  if(missing(these_sample_ids)){
    #don't load a bam. Assume the bam is already loaded
  }else{
    this_sample_id = paste0(these_sample_ids,collapse = ",")
  }
  if(missing(region) && missing(chrom)){
    stop("provide a region or coordinate as chrom, start, end")
  }

  if(missing(region)){
    if(missing(end)){
      end = as.numeric(start)+1
    }
    region = paste0(chrom, ":", as.numeric(start)-padding, "-", as.numeric(end) + padding)
  }
  filename = paste(region, gene, this_sample_id, details, "snapshot.png", sep = "--")
  if(!clobber){
    outfile = paste0(out_path,filename)
    #check if file exists already
    if(file.exists(outfile)){
      message(paste("file exists:",outfile,"skipping"))
      return()
    }
  }
  sock = socket
  if(missing(bams) & !missing(these_sample_ids)){
    meta = GAMBLR.helpers::handle_metadata(this_seq_type = this_seq_type) %>% dplyr::filter(sample_id %in% these_sample_ids)
    if(missing(genome_build)){
      genome_build = pull(meta,genome_build)
    }
    bam_path_pattern = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/{seq_type}_bams/{sample_id}.{genome_build}.bam"
    bams = mutate(meta,bam_path=glue::glue(bam_path_pattern)) %>% pull(bam_path)
    if(!length(bams)){

      message(paste("no bams found for",these_sample_ids))
      return()
    }
  }

  if(grepl("19",genome_build)||grepl("37",genome_build)){
    genome_build="hg19"
  }else{
    genome_build="hg38"
  }
  if(!missing(sample_ids)){
    if(exists("currently_loaded_bam") && currently_loaded_bam[1] == bams[1]){
        #skip loading
        print(paste("no need to load",bams))
      }
      else{
        IGVclear(sock)
        IGVgenome(sock, genome = genome_build)
        for(bam_file in bams){
          IGVload(sock, bam_file)
        }
        currently_loaded_bam <<- bams
    }

  }
  IGVgoto(sock, region)
  IGVsort(sock,sort_by)
  if(!missing(colour_by)){
    allowed = c("READ_STRAND","READ_GROUP","PAIR_ORIENTATION")
    if(colour_by %in% allowed){
      socketWrite(sock,paste("colorBy",colour_by, "\n"))
    }else{
      message(paste(colour_by,"must be one of",allowed))
    }
  }
  if(squish){
    socketWrite(sock,paste("squish", "\n"))
  }else{
    socketWrite(sock,paste("expand", "\n"))
  }
  if(viewaspairs){
    socketWrite(sock,paste("viewaspairs", "\n"))
  }
  filename = paste(this_sample_id, region, "snapshot.png", sep = "_")
  IGVsnapshot(sock, fname = filename, dirname = out_path)
  return(paste0(out_path, filename))
}
