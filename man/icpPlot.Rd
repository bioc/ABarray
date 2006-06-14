\name{icpPlot}
\alias{icpPlot}
\title{ icp plot function }
\description{
  QC plot for internal control probes
}
\usage{
icpPlot(controlData, colProbeID = 1, plotWhat = "Signal", pdfDir, jpgDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlData}{ Signal intensity matrix for icp probes }
  \item{colProbeID}{ the column where probeID is located }
  \item{plotWhat}{ Whether we are plotting signal or S/N }
  \item{pdfDir}{ a directory where pdf files should be produced }
  \item{jpgDir}{ a directory where jpg or bmp files should be produced }
}

\value{
  A series of QC plots
}
\author{ Yongming Sun }

\examples{
##---- Do not call this function DIRECTLY !! ----
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
