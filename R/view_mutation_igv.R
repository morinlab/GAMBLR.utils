#' @title View a variant in IGV
#'
#' @description Load bam(s) and view the context around a mutation
#'
#' @details Load bam(s) and view the context around a mutation.
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
#' Then obtain a socket and run this function as per the example.
#'
#' @param this_mutation Specify the mutation of interest in MAF format.
#' @param this_seq_type Specify the seq type, default is genome.
#' @param igv_port Specify the port IGV is listening on. Default: 60506 (optional if using the default).
#' @param socket Provide the socket variable obtained by running this function with no arguments.
#' @param sort_by base, quality, sample or readGroup.
#' @param colour_by Specify how IGV should colour the reads (see IGV documentation).
#' @param squish Force reads to be squished (see IGV documentation). Default is FALSE.
#' @param viewaspairs Set to TRUE if you want the reads to be shown as pairs rather than independently (see IGV documentation), default is FALSE.
#'
#' @return Path to file (.png).
#'
#' @import SRAdb
#' @export
#'
#' @examples
#' \dontrun{
#' socket = make_igv_snapshot() #run with no arguments to get the socket for a running IGV instance
#' this_mutation = GAMBLR.data::sample_data$grch37$maf %>%
#'  head(1)
#' view_mutation_igv(this_mutation,
#'                   socket = socket,
#'                   this_seq_type = "capture",
#'                   colour_by = "READ_STRAND",
#'                   squish = TRUE,
#'                   viewaspairs = TRUE)
#' }
#'
view_mutation_igv = function(this_mutation,
                             this_seq_type = "genome",
                             igv_port = 60506,
                             socket,
                             sort_by = "base",
                             colour_by,
                             squish = FALSE,
                             viewaspairs = FALSE){
  if(missing(socket)){
    print("returning socket for future use")
    sock = IGVsocket(port = igv_port)
    if(exists("currently_loaded_bam")){
      currently_loaded_bam <<- "" #avoid the function mistakenly thinking IGV still has a bam loaded
    }
    return(sock)
  }
  this_sample_id = unique(this_mutation$Tumor_Sample_Barcode)
  if(length(this_sample_id)>1){
    stop("provide a MAF with only one sample_id")
  }

  start = pull(this_mutation,Start_Position)
  end = start
  chrom = pull(this_mutation,Chromosome)
  region = paste0(chrom,":",start-50,"-",end+50)
  sock = socket

    meta = GAMBLR.helpers::handle_metadata(this_seq_type = this_seq_type) %>%
      dplyr::filter(sample_id %in% this_sample_id)

      genome_build = pull(meta,genome_build)

    bam_path_pattern = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/{seq_type}_bams/{sample_id}.{genome_build}.bam"
    bams = mutate(meta,bam_path=glue::glue(bam_path_pattern)) %>% pull(bam_path)
    if(!length(bams)){
      message(paste("no bams found for",sample_ids))
      return()
    }


  if(grepl("19",genome_build)||grepl("37",genome_build)){
    genome_build="hg19"
  }else{
    genome_build="hg38"
  }
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
  }
  if(viewaspairs){
    socketWrite(sock,paste("viewaspairs", "\n"))
  }
}
