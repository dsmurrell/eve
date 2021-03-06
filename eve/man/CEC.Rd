% Generated by roxygen2 (4.0.1): do not edit by hand
\name{CEC}
\alias{CEC}
\title{confidence–error correlation is a measure of the performance of the error estimates against the absolute errors}
\usage{
CEC(x, y)
}
\arguments{
\item{x}{the first variable}

\item{y}{the second variable}
}
\value{
the CEC value
}
\description{
CEC uses the Pearson's correlation between the sigma estimates of the assumed normal error distributions
and the absolute errors. This is normalised by the perfect confidence estimator
defined as the Pearson's correlation between the sorted sigmas and the sorted absolute errors.
CEC = cor(sigmas, |errors|) / cor(sort(sigmas), sort(|errors|))
}
\author{
Daniel Murrell <dsmurrell@gmail.com>
}
\references{
\url{http://dx.doi.org/10.1371/journal.pone.0048723}
}

