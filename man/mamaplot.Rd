\name{mamaplot}
\alias{mamaplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MA plot function }
\description{
  plot MA from vectors A and M
}
\usage{
mamaplot(A, M, idx, subset = sample(1:length(M), min(c(10000, length(M)))), span = 2/3, family.loess = "gaussian", cex = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ vector of average signal }
  \item{M}{ vector of difference signal }
  \item{idx}{ index for which S/N < 3 }
  \item{subset}{ subset }
  \item{span}{ span }
  \item{family.loess}{ loess fit }
  \item{cex}{ cex value }
  \item{...}{ additional arguments }
}

\value{
  MA plot
}
\references{ bioconductor affy package }
\author{ Yongming Sun }
\note{ Modified from bioconductor affy package }

\seealso{ See Also as \code{\link{mvaPair2}} }
\examples{
##---- Do not call this function DIRECTLY !! ----
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

