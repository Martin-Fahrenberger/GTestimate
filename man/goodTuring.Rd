% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/goodTuring.R
\name{goodTuring}
\alias{goodTuring}
\title{goodTuring}
\usage{
goodTuring(r, n, conf = 1.96)
}
\arguments{
\item{r:}{vector of observed frequencies}

\item{n:}{vector of frequencies of observed frequencies}

\item{conf:}{confidence factor for internal fit as described in Gale and Sampson (1995)}
}
\value{
Good-Turing frequency estimates for the frequencies given in r
}
\description{
Wrapper function for the Simple Good-Turing algorithm for frequency estimation
as described by Gale and Sampson (1995).
}
