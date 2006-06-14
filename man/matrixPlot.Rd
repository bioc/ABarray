\name{matrixPlot}
\alias{matrixPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ heatmap for matrix }
\description{
  Create heatmap from a matrix
}
\usage{
matrixPlot(x, nrgcols = 50, rlabels = TRUE, clabels = TRUE, rcols = 1, ccols = 1, k = 10, title = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a matrix }
  \item{nrgcols}{ number of colors to use }
  \item{rlabels}{ whether to use row labels }
  \item{clabels}{ whether to use column labels }
  \item{rcols}{ use supplemental row label }
  \item{ccols}{ use supplemental column label }
  \item{k}{ number of tick labels for scale bar }
  \item{title}{ title for the plot }
  \item{...}{ additional argument }
}
\details{
  This function can be used to plot any numberic matrix, e.g.,
  correlation matrix, S/N matrix, signal intensity matrix, etc
}
\value{
  heatmap
}
\author{ Yongming Sun }

\examples{
%- matrixPlot(dataMatrix, title = "Some Interesting Title")
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

