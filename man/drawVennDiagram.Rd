\name{drawVennDiagram}
\alias{drawVennDiagram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw Venn Diagram }
\description{
  Drawing actual Venn diagram
}
\usage{
drawVennDiagram(object, names, mar = rep(0.5, 4), cex = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ VennCounts object produced by \code{VennCounts}, which
	 is numeric matrix with last column "Counts" giving counts for each
	 possible vector outcome }
  \item{names}{ optional character vector giving names for the sets }
  \item{mar}{ numeric vector of length 4 specifying the width of the
	 margins around the plot. This argument is passed to par. }
  \item{cex}{ numerical value giving the amount by which the contrast
	 names should be scaled on the plot relative to the default.plotting
	 text. See par. }
  \item{...}{ any other arguments are passed to plot }
}

\value{
  a plot of Venn Diagram
}
\references{ Bioconductor Limma package }
\author{ Yongming Sun }

\examples{
##---- Do not call this function directly !! ----
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

