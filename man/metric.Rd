% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metric.R
\name{metric}
\alias{metric}
\alias{metric.pframe_1d}
\alias{metric.default}
\title{Compute metrics}
\usage{
metric(obj, ...)

\method{metric}{pframe_1d}(obj, ...)

\method{metric}{default}(obj)
}
\arguments{
\item{obj}{Object of pframe_1d class}
}
\value{
Object of pframe_1d class
}
\description{
Compute metrics from values of object
}
\examples{
metric(resamp(update(project(init_x0(pframe_1d())))))
}
\seealso{
\code{\link{particle}}
}
\author{
Justin Thong \email{justinthong93@gmail.com}
}
