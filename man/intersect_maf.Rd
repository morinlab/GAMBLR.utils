% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intersect_maf.R
\name{intersect_maf}
\alias{intersect_maf}
\title{Intersect MAF.}
\usage{
intersect_maf(maf1, maf2, set_returned = "maf1_only")
}
\arguments{
\item{maf1}{First list of MAFs.}

\item{maf2}{Second list of MAFs.}

\item{set_returned}{List of MAFs that doesn't share the same start positions as the other list of MAFs. Accepted commands are; "maf1_only" and "maf2_only", default is "maf1_only".}
}
\value{
Set of MAFs with start positions that don't match the start positions in the other supplied MAF file.
}
\description{
Perform set operations on two MAFs.
}
\details{
Perform set operations on two MAFs.
}
\examples{
\dontrun{
intersected_mafs_l1 = intersect_maf(maf_list1, maf_list2, "maf1_only")
intersected_mafs_l2 = intersect_maf(maf_list1, maf_list2, "maf2_only")
}

}
