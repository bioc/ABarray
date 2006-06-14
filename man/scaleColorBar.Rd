\name{scaleColorBar}
\alias{scaleColorBar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create scale for heatmap }
\description{
  Create a bar for heatmap scales
}
\usage{
scaleColorBar(x, horizontal = FALSE, col = rgcolorsfunc(50), scale = 1:length(x), 
    k = 10, cLen = 9, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of scales need to be plotted }
  \item{horizontal}{ whether the bar is vertical or horizontal }
  \item{col}{ color function }
  \item{scale}{ scale of the bar }
  \item{k}{ number of intervals on scale }
  \item{cLen}{ length of columns }
  \item{...}{ additional arguments }
}

\value{
  none
}

\author{ Yongming Sun }

\examples{
##--- Do not call this function directly !! ----
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

