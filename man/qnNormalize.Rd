\name{qnNormalize}
\alias{qnNormalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform quantile normalization }
\description{
  Perform quantile normalization between arrays
}
\usage{
qnNormalize(eData, snr, method = 'quantile', snThresh = 3, ties = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eData}{ matrix of gene expression values }
  \item{snr}{ Optional signal/noise ratio. Only used for trimAMean method }
  \item{method}{ The normalization method desired. Default method is quantile }
  \item{snThresh}{ Signal/noise threshold (default = 3) to indicate presence or
    absence of a probe signal }
  \item{ties}{ handle values with same rank }
}

\details{
  This function performs various normalization for the array data. The default is
  quantile normalization method (adapted from Bioconductor limma package). Other
  normalization methods include median, mean, trimMean (trimmed mean), trimAMean
  (mean with absent gene removed).
  
  For the median normalizaiton, the median signal of each array is scaled to the same
  value (this value is calculated to equal to the median of all values in the data).
  The signal values for each array are then adjusted by the scaling factor.
  
  For the mean normalization, the approach is similar to the median normalization
  procedure except that the mean signal of each array is scaled to the same value
  (this value is median of all signals in the data).
  
  For the trimMean normalization, the approach is similar to the mean normalization
  except that the mean for each array is calculated after trimming the top and botton
  5\% of signals (a total of 10\% of values).
  
  For the trimAMean normalization, the signal values for absent probes are not
  considered. If the s/n of a probe is less than snThresh (default = 3), the expression
  of the probe is considered not present (absent). The remaining values are then
  trimmed (top and botton 2.5\%, a total of 5\%), and the mean value for each array
  after trimming is scaled to the same value (median of all values in the data).
}

\value{
  data matrix with quantile normalized data values
}
\references{ bioconductor limma package for quantile normalization }
\author{ Yongming Sun }

\examples{
%- qnNormData <- qnNormalize(data)
}
\keyword{ manip }% at least one, from doc/KEYWORDS

