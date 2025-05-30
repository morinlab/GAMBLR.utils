% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic_data.R
\name{create_genomic_data}
\alias{create_genomic_data}
\title{Create Genomic Data}
\usage{
create_genomic_data(data, genome_build, subclass = NULL)
}
\arguments{
\item{data}{A data frame containing the genomic data.}

\item{genome_build}{A string specifying the genome build
("grch37" or "hg38").}

\item{subclass}{Optional. A character string specifying
a subclass for the data.
For example, "maf_data", "bed_data", or "seg_data".
Default is \code{NULL}.}
}
\value{
A data frame with class attributes for genomic
data and a genome_build attribute.
}
\description{
This function creates a genomic data object from the given input data frame.
}
\details{
It attaches a \code{genome_build} attribute and assigns a class of
\code{"genomic_data"}.
Optionally, a subclass (e.g., \code{"maf_data"}, \code{"bed_data"},
\code{"seg_data"}, etc.) can be provided to further classify the
data but this handled by the other create_ functions.
}
\keyword{internal}
