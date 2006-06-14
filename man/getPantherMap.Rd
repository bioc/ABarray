\name{getPantherMap}
\alias{getPantherMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create pie chart for probes involved in Panther Pathway }
\description{
  Given a list of probeID, attempt to find out panther classification
  information
}
\usage{
getPantherMap(probeID, title, figDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probeID}{ a list of probeIDs }
  \item{title}{ the title for the figure to be generated }
  \item{figDir}{ directory for the figures to be placed in }
}

\value{
  None. Several figures will be generated.
}
\author{ Yongming Sun }

\examples{
%- getPantherMap(probeID, "plotName", "figDir")
}
\keyword{ manip }% at least one, from doc/KEYWORDS
