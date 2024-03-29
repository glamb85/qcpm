% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_checks.R
\name{get_checks}
\alias{get_checks}
\title{intternal checks}
\usage{
get_checks(data, inner, outer, scheme, tau, ...)
}
\arguments{
\item{data}{matrix or data frame containing the manifest variables.}

\item{inner}{A square (lower triangular) boolean matrix representing 
the inner model (i.e. the path relationships between latent variables).}

\item{outer}{list of vectors with column indices or column names
from \code{data} indicating the sets of manifest variables forming 
each block (i.e. which manifest variables correspond to each block).}

\item{scheme}{string indicating the type of inner weighting
scheme. It is equal to \code{"factorial"} by default.
 Possible values are \code{"centroid"} or  \code{"factorial"}.}

\item{tau}{if sepcifed indicates the specific quantile to be considered}

\item{\dots}{Further arguments passed on to \code{\link{get_checks}}.}
}
\value{
A list containing checked  parameters for internal estimation of 
the qcpm algorithm.
}
\description{
intternal checks
}
\details{
Internal function. \code{get_checks} is called by \code{qcpm}.
}
\keyword{internal}
