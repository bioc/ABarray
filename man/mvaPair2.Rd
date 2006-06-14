\name{mvaPair2}
\alias{mvaPair2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ plot MA for each pair of columns }
\description{
  MA plot for each pair of columns
}
\usage{
mvaPair2(x, y = NULL, snThresh = 3, labels = colnames(x), log.it = FALSE, span = 2/3, 
    family.loess = "gaussian", digits = 3, line.col = 2, main = "MA plot", ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ expression matrix }
  \item{y}{ S/N ratio matrix }
  \item{snThresh}{ S/N threshold }
  \item{labels}{ name for the labels }
  \item{log.it}{ should data be log transformed }
  \item{span}{ span of the plot }
  \item{family.loess}{ curve fitting }
  \item{digits}{ number of digits to display }
  \item{line.col}{ size of the line col }
  \item{main}{ title for the MA plot }
  \item{...}{ additional argument }
}
\details{
  If S/N ratio is available, probes with S/N < 3 in both array will be
  colored differently.
}
\value{
  MA plot
}
\author{ Yongming Sun }

\examples{
##---- exprs expression matrix, sn s/n ratio !! ----
%- mvaPairs(exprs, sn)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

