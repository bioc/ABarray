\name{lpe.fdr.BH}
\alias{lpe.fdr.BH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform FDR on LPE results }
\description{
  Perform Benjamini and Hochberg FDR adjustment on LPE results
}
\usage{
lpe.fdr.BH(lpe.result, adjp = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpe.result}{ the result from LPE analysis }
  \item{adjp}{ Type of adjustment, default "BH" }
}
\details{
  Do not call this function directly. Called from doLPE
}
\value{
  a matrix with original and ajusted p values
}
\references{ Bioconductor LPE package }
\author{ Yongming Sun }

\examples{
##---- Do not call this function directly !! ----
}
\keyword{ manip }% at least one, from doc/KEYWORDS

