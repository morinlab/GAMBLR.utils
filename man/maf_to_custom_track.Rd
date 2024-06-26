% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maf_to_custom_track.R
\name{maf_to_custom_track}
\alias{maf_to_custom_track}
\title{Maf To Custom Track.}
\usage{
maf_to_custom_track(
  maf_data,
  these_sample_ids = NULL,
  these_samples_metadata = NULL,
  this_seq_type = "genome",
  output_file,
  as_bigbed = FALSE,
  colour_column = "lymphgen",
  as_biglolly = FALSE,
  track_name = "GAMBL mutations",
  track_description = "mutations from GAMBL",
  verbose = FALSE,
  padding_size = 0,
  projection = "grch37",
  bedToBigBed_path = "config"
)
}
\arguments{
\item{maf_data}{Either a maf loaded from disk or from the database using a get_ssm function.}

\item{these_sample_ids}{A vector of sample IDs to subset the samples of interest from the input \code{maf_data}. If NULL (the default), all samples in \code{maf_data} are kept.}

\item{these_samples_metadata}{A metadata table to subset the samples of interest from the input \code{maf_data}. If NULL (the default), all samples in \code{maf_data} are kept.}

\item{this_seq_type}{The seq type you want back, default is "genome".}

\item{output_file}{Name for your new bed file that can be uploaded as a custom track to UCSC.}

\item{as_bigbed}{Boolean parameter controlling the format of the returned file. Default is FALSE.}

\item{colour_column}{Set the colouring properties of the returned bed file. Per default, this function will assign colour based on "lymphgen".}

\item{as_biglolly}{Boolean parameter controlling the format of the returned file. Default is FALSE (i.e a BED file will be returned).}

\item{track_name}{Track name. Default is "GAMBL mutations"}

\item{track_description}{Track description. Default is "mutations from GAMBL"}

\item{verbose}{Default is FALSE.}

\item{padding_size}{Optional parameter specifying the padding size in the returned file, default is 0.}

\item{projection}{Specify which genome build to use. Possible values are "grch37" (default) or "hg38". This parameter has
effect only when \code{as_bigbed} or \code{as_biglolly} is TRUE.}

\item{bedToBigBed_path}{Path to your local \code{bedToBigBed} UCSC tool or the string
\code{"config"} (default). If set to \code{"config"}, \code{GAMBLR.helpers::check_config_value}
is called internally and the \code{bedToBigBed} path is obtained from the \code{config.yml}
file saved in the current working directory. This parameter is ignored if both
\code{as_bigbed} and \code{as_biglolly} is set to \code{FALSE}.}
}
\value{
Nothing.
}
\description{
Convert a maf-formatted data frame into a bed custom track file for UCSC.
}
\details{
This function takes an incoming MAF and converts it to a UCSC Genome Browser ready BED (or bigbed/biglolly) file.
Optional parameters available for further customization of the returned file. For more information, refer to the parameter descriptions and function examples.
}
\examples{
library(GAMBLR.data)
region_myc <- dplyr::filter(grch37_lymphoma_genes_bed, hgnc_symbol == "MYC")
my_maf <- get_ssm_by_regions(regions_bed = region_myc, streamlined = FALSE)
maf_to_custom_track(maf_data = my_maf, output_file = "../mutations.bed")

}
