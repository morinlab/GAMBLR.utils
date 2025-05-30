% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_GISTIC_inputs.R
\name{prepare_gistic_inputs}
\alias{prepare_gistic_inputs}
\title{Prepare GISTIC input files}
\usage{
prepare_gistic_inputs(
  these_samples_metadata,
  projection = "grch37",
  output_dir = "/Users/rmorin/Desktop/GISTIC2/",
  flavour = "combined",
  overwrite = FALSE
)
}
\arguments{
\item{these_samples_metadata}{A data frame with metadata.}

\item{projection}{The projection to operate on. Default is grch37.}

\item{output_dir}{The path to the directory to write outputs.}

\item{flavour}{Additional argument.}

\item{overwrite}{Whether to overwrite the existing outputs. Default is FALSE.}
}
\value{
nothing
}
\description{
Prepare GISTIC input files
}
\examples{
\dontrun{
BL_meta = get_gambl_metadata() \%>\% 
  dplyr::filter(pathology=="BL")
prepare_gistic_inputs(BL_meta,
                      "hg38",
                      "/Users/rmorin/Desktop/GISTIC2/BL/",
                      "battenberg")


FL_meta = get_gambl_metadata() \%>\% 
  dplyr::filter(pathology=="FL",seq_type!="mrna")
prepare_gistic_inputs(FL_meta,
                      "grch37",
                      "/Users/rmorin/Desktop/GISTIC2/FL/",
                      "battenberg")
}
}
