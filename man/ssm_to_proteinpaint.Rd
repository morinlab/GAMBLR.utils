% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssm_to_proteinpaint.R
\name{ssm_to_proteinpaint}
\alias{ssm_to_proteinpaint}
\title{MAF to ProteinPaint table}
\usage{
ssm_to_proteinpaint(
  maf_data,
  these_sample_ids = NULL,
  these_samples_metadata = NULL,
  this_seq_type = "genome",
  sample_type = "time_point",
  coding_only = FALSE,
  these_genes = GAMBLR.data::lymphoma_genes_comprehensive$Gene,
  debug_flag = FALSE
)
}
\arguments{
\item{maf_data}{A data frame in MAF format.}

\item{these_sample_ids}{A vector of sample IDs to subset the samples of interest from the
input \code{maf_data}. If NULL (the default), samples are not filtered by this parameter.}

\item{these_samples_metadata}{A metadata table, with sample IDs in a column, to subset
the samples of interest from the input \code{maf_data}. Also, important information from the
metadata is added to the output table. If NULL (the default), \code{ssm_to_proteinpaint}
internally creates the metadata table with samples according to the parameters
\code{these_sample_ids} and \code{this_seq_type}.}

\item{this_seq_type}{The seq type to get results for. One of "genome" (default) or
"capture".}

\item{sample_type}{A single sting with the name of the optional column (from \code{maf_data}
or the metadata) that you want to use to distinguish multiple samples from the same
patient. If NULL, the ProteinPaint visualization is split only by pathologies. If not
NULL, the visualization will be split by each combination between pathologies and the
column specified by \code{sample_type}. The default is "time_point".}

\item{coding_only}{Boolean parameter. Set to TRUE to restrict to only coding mutations.
The default is FALSE.}

\item{these_genes}{A vector of strings with name of genes that you want results for.
If NULL, all genes of the input \code{maf_data} are kept. The default is all lymphoma genes.}

\item{debug_flag}{Boolean parameter. Set to TRUE for returning rows from the
incoming MAF that do not contain any values in the required columns. Commonly used for
checking purposes only. Setting this to TRUE, does not produce an output compatible with
Protein paint. The default is FALSE.}
}
\value{
A data frame.
}
\description{
This function takes a MAF-like data frame and convert it to the right format
for visualization with ProteinPaint.
}
\details{
For visualization with ProteinPaint, the output data frame must be saved to a
file and uploaded to https://proteinpaint.stjude.org/.

For a valid ProteinPaint table for visualization, the MAF table must contain the required
columns to be converted to the ProteinPaint format. They are: Hugo_Symbol, RefSeq,
Chromosome, Start_Position, HGVSp_Short, Variant_Classification. Optional columns used
from the MAF table are: Mutation_Status, t_alt_count, t_depth, n_alt_count, n_depth,
Tumor_Sample_Barcode, Reference_Allele. In the case of using functions \code{get_ssm_by_samples},
\code{get_coding_ssm}, or \code{get_ssm_by_patients} to get your MAF table and required/desired
columns are missing, consider using parameter \code{basic_columns = FALSE} to ensure that your
MAF has all needed columns. Other columns are taken from the metadata provided by the
\code{these_samples_metadata} parameter or generated internally.
}
\examples{
library(GAMBLR.data)

# define same somples
my_samples = c("DOHH-2", "OCI-Ly10", "OCI-Ly3", "SU-DHL-10", "SU-DHL-4")

# get maf data frame
my_maf = get_ssm_by_samples(these_sample_ids = my_samples, 
                            this_seq_type = "genome", basic_columns = FALSE)

# convert maf to ProteinPaint format
pp_df = ssm_to_proteinpaint(maf_data = my_maf, this_seq_type = "genome")

# convert only coding mutations to ProteinPaint format
pp_df = ssm_to_proteinpaint(maf_data = my_maf, this_seq_type = "genome",
                            coding_only = TRUE)

}
