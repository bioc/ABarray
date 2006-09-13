\name{hclusterPlot}
\alias{hclusterPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ heatmap generation }
\description{
  plot clustering heatmap using correlation
}
\usage{
hclusterPlot(expr, title, dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{ matrix of gene expression value }
  \item{title}{ the title for the plot }
  \item{dist}{ whether to use correlation or distance for clustering,
	 default to use Euclidean distance. Use dist = "Correlation" to
	 cluster with correlation coefficient }
}
\details{
  generating heatmap using correlation as distance
}
\value{
  None. heatmap will be generated.
}
\author{ Y Andrew Sun }

\examples{
%- hclusterPlot(exprMatrix, "Some title for plot")
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

