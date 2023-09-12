#' @title Annotate mutations target AID motif
#'
#' @description Checks for the presence of mutations at a given AID motif
#'
#' @details In positions that reference allele  "C" has been mutated, two nucleotides before and one nucleotide after it 
#' will be captured. In positions that reference allele  "G" has been mutated, one nucleotide before and two nucleotides after it 
#' will be captured. The other mutations will have "NO" for the sequence. 
#' Then, the forward orientation and reverse compliment of the motif will be generated and their presence in aquired reference
#' allele "C" and reference allele "G" sequence will be checked, respectively. ("NO" for the rows that have different reference
#' alleles) 
#'
#' @param maf MAF data frame (required columns: Reference_Allele, Chromosome, Start_Position, End_Position)
#' @param motif The motif sequence (default is WRCY)
#' @param projection The genome build projection for the variants you are working with (default is grch37)
#' @param fastaPath Can be a path to a FASTA file
#'
#' @return A data frame with two extra columns (seq and motif).
#'
#' @rawNamespace import(IRanges, except = c("start", "end", "merge", "shift", "collapse", "union", "slice", "intersect", "setdiff", "desc", "reduce"))
#' @rawNamespace import(GenomicRanges, except = c("start", "end", "merge", "shift", "union", "intersect", "setdiff", "reduce"))
#' @import Rsamtools readr dplyr
#' @export
#'
#' @examples
#' annotate_ssm_motif_context (maf = maf,
#'           motif = "WRCY",
#'           projection = "hg38"
#'          )


#This function check that if the motif pattern is in the sequence
annotate_ssm_motif_context <- function(maf,
                     motif = "WRCY",
                     projection = "grch37",
                     fastaPath
                     ){
    if (projection == "grch37") {
        maf$Chromosome <- gsub("chr", "", maf$Chromosome)
    } else {
        # If there is a mix of prefixed and non-prefixed options
        maf$Chromosome <- gsub("chr", "", maf$Chromosome) 
        maf$Chromosome <- paste0("chr", maf$Chromosome)
    }
    # If there is no fastaPath, it will read it from config key 
    # Based on the projection the fasta file which will be loaded is different
    if (missing(fastaPath)){
        base <- check_config_value(config::get("repo_base"))
        fastaPath <- paste0(
            base,
            "ref/lcr-modules-references-STABLE/genomes/",
            projection,
            "/genome_fasta/genome.fa"
        )
    }
    # It checks for the presence of a local fastaPath
    if (!file.exists(fastaPath)) {
        stop("Failed to find the fasta file")
    }
    # It checks that if the user considers "WRCY" as the motif
    if (!motif == "WRCY"){
        stop("It is only currently available for WRCY. ",
             "To use it for another motif, comment \"if\" code for the motif"
        )
    }
    # Create a reference to an indexed fasta file.
    fasta <- Rsamtools::FaFile(file = fastaPath)
    # This section provides the sequence
    # If the reference allele is C, it will return the one nucleotide after it
    # and 2 nucleotides before it
    sequences <- maf %>%
        dplyr::mutate(
            seq = ifelse(
                Reference_Allele == "C",
                    Rsamtools::getSeq(
                        fasta,
                        GenomicRanges::GRanges(
                            maf$Chromosome,
                            IRanges(
                                start = maf$Start_Position - 2,
                                end = maf$End_Position + 1
                            )
                        )
                    ),
                     # If the reference allele is G, it will return 
                     # 2 nucleotides after it and the one nucleotide before it
                     ifelse(
                         Reference_Allele == "G",
                         Rsamtools::getSeq(
                             fasta,
                             GenomicRanges::GRanges(
                                 maf$Chromosome,
                                 IRanges(
                                     start = maf$Start_Position - 1,
                                     end = maf$End_Position + 2
                                 )
                             )
                          ),
                           # In other cases it returns "NO" for that mutation
                           "NO"
                     )
            )
        )
  
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
    word <- motif
    splitWord <- strsplit(word,"")[[1]] # Split the word into its letters
    splitWordLen <- length(splitWord)
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

    #This section checks for the presence of the motif in the sequence
    # If the reference allele is C, it will check the forward orientation with
    # its sequence
    finder <- sequences %>%
        dplyr::mutate(
            !!motif := ifelse(
              Reference_Allele == "C",
              stringr::str_detect(
                  sequences$seq, strForMotif
              ),
              # If the reference allele is G, it will check the reverse-
              # compliment orientation with its sequence
              ifelse(
                  Reference_Allele == "G",
                  stringr::str_detect(
                      sequences$seq, strRevComp
                  ),
                   # For the rest of mutations it will return "NO"
                   "NO"
              )
            )
        )
    return(finder)
}
