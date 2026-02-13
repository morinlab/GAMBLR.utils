#' @title Annotate mutations target motif
#'
#' @description Checks for the presence of mutations at a given motif
#'
#' @details In positions that reference allele has been mutated, it will capture (motif length - 1) before and (motif length + 1) alleles after the mutated position.
#' Then, it looks for the presence of motif in the captured sequence and check if the mutation has occurred in the indexed position, it will return 
#' SITE and if the the motif is present, but the mutation is not in the indexed position, it will return MOTIF. In other cases, it will return FALSE.
#' NA will be shown if the mutation is an indel mutation.
#'
#' @param maf MAF data frame (required columns: Reference_Allele, Chromosome, Start_Position, End_Position)
#' @param motif The motif sequence (default is WRCY)
#' @param index Position of the mutated allele in the motif
#' @param genome_build The genome build for the variants you are working with (default is grch37)
#' @param fastaPath Can be a path to a FASTA file
#'
#' @return A data frame with two extra columns (seq and motif).
#'
#' @rawNamespace import(IRanges, except = c("start", "end", "merge", "shift", "collapse", "union", "slice", "intersect", "setdiff", "desc", "reduce", "trim"))
#' @rawNamespace import(GenomicRanges, except = c("start", "end", "merge", "shift", "union", "intersect", "setdiff", "reduce", "trim"))
#' @import Rsamtools readr dplyr BSgenome
#' @export
#'
#' @examples
#' \dontrun{
#' my_maf <- GAMBLR.open::get_coding_ssm() %>% 
#'   dplyr::filter(Hugo_Symbol=="BCL2") %>%
#'   dplyr::arrange(Chromosome,Start_Position,Tumor_Sample_Barcode) %>%
#'   head()
#' 
#' annotated = annotate_ssm_context(maf = my_maf,
#'                                  motif = "WRCY",
#'                                  genome_build="grch37")
#' 
#' dplyr::select(annotated,1,5,6,11,13,16,seq,WRCY)
#'}
annotate_ssm_context <- function(maf,
                                       motif = "WRCY",
                                       index = 3,
                                       genome_build,
                                       fastaPath
){

    bsgenome_loaded = FALSE
    # If there is no fastaPath, it will read it from config key 
    # Based on the projection the fasta file which will be loaded is different
    
    # It checks for the presence of a local fastaPath
    if (missing(fastaPath)) {
        print("MISSING fastaPath, trying to load BSgenome package")
        #try BSgenome
      installed = BSgenome::installed.genomes()
      if(genome_build=="hg38"){
        bsgenome_name = "BSgenome.Hsapiens.UCSC.hg38"
      }else if(genome_build == "grch37"){
        bsgenome_name = "BSgenome.Hsapiens.UCSC.hg19"
      }else{
        stop(paste("unsupported genome:",genome_build))
      }
      if(bsgenome_name %in% installed){
          genome = BSgenome::getBSgenome(bsgenome_name)
          bsgenome_loaded = TRUE
      }else{
        print(installed)
        print(paste("Local Fasta file cannot be found and missing genome_build",
            bsgenome_name,"Supply a fastaPath for a local fasta file or install the missing BSGenome package and re-run"))
      }
    }else {
      
        stop(paste("fastaPath",
                    fastaPath,
                    "file cannot be found. Please provide a fastaPath to an existing fasta file"))
      
    }
    word <- motif
    splitWord <- strsplit(word,"")[[1]] # Split the word into its letters
    splitWordLen <- length(splitWord)

    if(!bsgenome_loaded){
      # Create a reference to an indexed fasta file.
      if (!file.exists(fastaPath)) {
        stop("Failed to find the fasta file and no compatible BSgenome found")
      }
      fasta = Rsamtools::FaFile(file = fastaPath)
    }
    if(bsgenome_loaded) {
      if(genome_build == "grch37") {
        maf = mutate(maf,original_chrom = Chromosome)
        maf = mutate(maf, Chromosome = paste0("chr",Chromosome))
      }
      sequences <- maf %>%
        dplyr::mutate(
          seq = as.character(
                Rsamtools::getSeq(
                    genome,
                    GenomicRanges::GRanges(
                      maf$Chromosome,
                      IRanges(
                        start = maf$Start_Position - (splitWordLen - 1),
                        end = maf$End_Position + (splitWordLen - 1)
                      )
                    )
                )
            )
        )
    } else{
      # This section provides the sequence
      # It will return one allele less than the length of motif before and after the indexed allele
      sequences <- maf %>%
        dplyr::mutate(
            seq = as.character(
                Rsamtools::getSeq(
                    fasta,
                    GenomicRanges::GRanges(
                        maf$Chromosome,
                        IRanges(
                            start = maf$Start_Position - (splitWordLen - 1),
                            end = maf$End_Position + (splitWordLen - 1)
                        )
                    )
                ))
        )
    }
    # This section provides motif and its reverse complement 
    compliment <- c(
        'A'= 'T',
        'T'= 'A',
        'C'= 'G',
        'G'= 'C',
        'U'= 'A',
        'R'= 'Y',
        'Y'= 'R',
        'S'= 'S',
        'W'= 'W',
        'K'= 'M',
        'M'= 'K',
        'B'= 'V',
        'D'= 'H',
        'H'= 'D',
        'V'= 'B',
        'N'= 'N'
    )
    IUPACCodeDict <- c(
        'A'= 'A',  # Adenine
        'C'= 'C',  # Cytosine
        'G'= 'G',  # Guanine
        'T'= 'T',  # Thymine
        'R'= '[AG]',  # A or G
        'Y'= '[CT]',  # C or T
        'S'= '[GC]',  # G or C
        'W'= '[AT]',  # A or T
        'K'= '[GT]',  # G or T
        'M'= '[AC]',  # A or C
        'B'= '[CGT]',  # C or G or T
        'D'= '[AGT]',  # A or G or T
        'H'= '[ACT]',  # A or C or T
        'V'= '[ACG]',  # A or C or G
        'N'= '[ACGT]'  # any base
    )
    forMotif <- character(splitWordLen) # forMotif, the same length as splitWord
    for (i in seq_along(splitWord)){
        # Convert incomplete nuc specification into their different nucleotides
        if (splitWord[i] %in% names(IUPACCodeDict)){
            forMotif[i] = IUPACCodeDict[splitWord[i]]
        }
    }
    # Collapsing all the letters of forward orientation and make it 
    # into a single string
    strForMotif <- paste(forMotif, collapse = "")
    RevCompMotif <- character(splitWordLen)
    for (i in seq_along(splitWord)){
        if (splitWord[i] %in% names(compliment)){
            # Provide the complement version of the motif
            RevCompMotif[i] = compliment[splitWord[i]]
        }
    }
    # Provide the reverse version of RevCompMotif
    vecRevComp <- rev(RevCompMotif)
    IUPACRevCompMotif <- character(splitWordLen)
    for (i in seq_along(vecRevComp)){
        if (vecRevComp[i] %in% names(IUPACCodeDict)){
            IUPACRevCompMotif[i] = IUPACCodeDict[vecRevComp[i]]
        } 
    }
    # Collapsing all the letters of reverse complement orientation and make it 
    # into a single string
    strRevComp <- paste(IUPACRevCompMotif, collapse = "")
    # This section checks for the presence of the motif in the sequence
    # Check if the desired mutation in the indexed postion is available in the sequence (If TRUE, it will return "SITE")
    finder <- sequences %>%
        dplyr::mutate(
            !!motif := case_when(
                Reference_Allele %in% unlist(
                    strsplit(
                        IUPACCodeDict[splitWord[[index]]], ""
                    )
                ) &
                    stringr::str_detect(
                        map_chr(
                            strsplit(
                                sequences$seq, ""
                            ), 
                            ~paste(.[((splitWordLen - index) + 1) : ((2 * splitWordLen)-index)], collapse = "")
                        ), strForMotif
                    )~ "SITE",
                Reference_Allele %in% unlist(
                    strsplit(
                        IUPACCodeDict[compliment[splitWord[[index]]]], ""
                    )
                ) &
                    stringr::str_detect(
                        map_chr(
                            strsplit(
                                sequences$seq, ""
                            ),
                            ~paste(.[(index) : (index + splitWordLen -1)], collapse = "")
                        ), strRevComp
                    ) ~ "SITE",
                # For indel mutations it will return "NA"
                Reference_Allele == "-" | Tumor_Seq_Allele2 == "-" ~ "NA",
                # If mutation is in the motif, but it is not in the indexed position, it will return "MOTIF"
                (stringr::str_detect(
                    sequences$seq, strForMotif
                ) |
                    stringr::str_detect(
                        sequences$seq, strRevComp
                    )
                ) ~ "MOTIF",
                # For the rest of mutations it will return "FALSE"
                TRUE ~ "FALSE"
            )
        )
    if("original_chrom" %in% colnames(finder)){
        finder = mutate(finder,Chromosome = original_chrom) %>% 
        dplyr::select(-original_chrom)
    }
    return(finder)
}