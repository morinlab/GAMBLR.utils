#' @title Fread MAF.
#'
#' @description Read a MAF into R with options for what columns to keep (default all columns will be returned).
#'
#' @details This function lets the user specify a path to a MAF file on disk, this function then reads this MAF into R.
#' The user has the option to keep specific columns from the incoming MAF file (default is to keep all columns).
#' To return all available columns, set `return_cols = TRUE`. For specifying columns of interest, refer to the `select_cols` parameter.
#' For more information, refer to the parameter descriptions as well as function examples.
#'
#' @param maf_file_path Path to maf that is about to be read.
#' @param select_cols Optional parameter for specifying what columns are to be returned. If not specified, all columns will be kept.
#' @param return_cols Optional parameter for returning the names of all available MAF columns. If set to TRUE a character vector of column names will be returned, and no MAF will be read. Default is FALSE.
#'
#' @return A data frame containing MAF data from a MAF file.
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' \dontrun{
#' #read a maf into R with all columns kept
#' my_maf = fread_maf(maf_file_path = "some_directory/this_is_a.maf")
#'
#' #return what columns are available
#' maf_cols = fread_maf(return_cols = TRUE)
#'
#' #read maf with only a selection of columns
#' my_maf = fread_maf(maf_file_path = "some_directory/this_is_a.maf",
#'                    select_cols = c(Hugo_Symbol = "character",
#'                                    Chromosome = "character",
#'                                    Start_Position = "integer",
#'                                    End_Position = "integer",
#'                                    Variant_Type = "character"))
#' }
#'
fread_maf = function(maf_file_path,
                     select_cols,
                     return_cols = FALSE){

  colClasses=c(Hugo_Symbol="character",
               Entrez_Gene_Id="integer",
               Center="character",
               NCBI_Build="character",
               Chromosome="character",
               Start_Position="integer",
               End_Position="integer",
               Strand="character",
               Variant_Classification="character",
               Variant_Type="character",
               Reference_Allele="character",
               Tumor_Seq_Allele1="character",
               Tumor_Seq_Allele2="character",
               dbSNP_RS="character",
               dbSNP_Val_Status="logical",
               Tumor_Sample_Barcode="character",
               Matched_Norm_Sample_Barcode="character",
               Match_Norm_Seq_Allele1="character",
               Match_Norm_Seq_Allele2="character",
               Tumor_Validation_Allele1="logical",
               Tumor_Validation_Allele2="logical",
               Match_Norm_Validation_Allele1="logical",
               Match_Norm_Validation_Allele2="logical",
               Verification_Status="logical",
               Validation_Status="logical",
               Mutation_Status="logical",
               Sequencing_Phase="logical",
               Sequence_Source="logical",
               Validation_Method="logical",
               Score="logical",
               BAM_File="logical",
               Sequencer="logical",
               Tumor_Sample_UUID="logical",
               Matched_Norm_Sample_UUID="logical",
               HGVSc="character",
               HGVSp="character",
               HGVSp_Short="character",
               Transcript_ID="character",
               Exon_Number="character",
               t_depth="integer",
               t_ref_count="integer",
               t_alt_count="integer",
               n_depth="integer",
               n_ref_count="integer",
               n_alt_count="integer",
               all_effects="character",
               Allele="character",
               Gene="character",
               Feature="character",
               Feature_type="character",
               Consequence="character",
               cDNA_position="character",
               CDS_position="character",
               Protein_position="character",
               Amino_acids="character",
               Codons="character",
               Existing_variation="character",
               ALLELE_NUM="integer",
               DISTANCE="numeric",
               STRAND_VEP="numeric", #for some reason, many are -1.0 or 1.0 instead of just -1 or 1. Numeric seems to work
               SYMBOL="character",
               SYMBOL_SOURCE="character",
               HGNC_ID="character",
               BIOTYPE="character",
               CANONICAL="character",
               CCDS="character",
               ENSP="character",
               SWISSPROT="character",
               TREMBL="character",
               UNIPARC="character",
               RefSeq="character",
               SIFT="character",
               PolyPhen="character",
               EXON="character",
               INTRON="character",
               DOMAINS="character",
               AF="character",
               AFR_AF="character",
               AMR_AF="character",
               ASN_AF="logical",
               EAS_AF="character",
               EUR_AF="character",
               SAS_AF="character",
               AA_AF="character",
               EA_AF="character",
               CLIN_SIG="character",
               SOMATIC="character",
               PUBMED="character",
               MOTIF_NAME="character",
               MOTIF_POS="numeric",
               HIGH_INF_POS="character",
               MOTIF_SCORE_CHANGE="logical",
               IMPACT="character",
               PICK="numeric", #for some reason, many are -1.0 or 1.0 instead of just -1 or 1. Numeric seems to work
               VARIANT_CLASS="character",
               TSL="character", # has explicit NA so must be read as character to avoid warnings
               HGVS_OFFSET="character", # has explicit NA so must be read as character to avoid warnings
               PHENO="character",
               MINIMISED="logical",
               GENE_PHENO="numeric",
               FILTER="character",
               flanking_bps="character",
               vcf_id="character",
               vcf_qual="character",
               gnomAD_AF="logical",
               gnomAD_AFR_AF="logical",
               gnomAD_AMR_AF="logical",
               gnomAD_ASJ_AF="logical",
               gnomAD_EAS_AF="logical",
               gnomAD_FIN_AF="logical",
               gnomAD_NFE_AF="logical",
               gnomAD_OTH_AF="logical",
               gnomAD_SAS_AF="logical",
               vcf_pos="integer",
               gnomADg_AF="character",
               blacklist_count="numeric")

  if(missing(select_cols)){
    #get all column names from the variable that defines the classes
    select_cols = names(colClasses)
    if(return_cols){
      return(select_cols)
    }
  }
  maf_df = read_tsv(
    file = maf_file_path
    ) %>%
    dplyr::select(any_of(select_cols))

  return(maf_df)
}
