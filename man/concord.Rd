\name{concord}
\alias{concord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate signal detection concordance }
\description{
  Calculate signal detection concordance between columns
  using S/N threshold (default = 3)
}
\usage{
concord(sn, snThresh = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sn}{ a matrix containing s/n ratio }
  \item{snThresh} {S/N threshold to use, default = 3}
}
\value{
  a matrix with the concordance
}

\author{ Y Andrew Sun }

\examples{
#-concordance <- concord(sn)  ##- sn ratio matrix
}
\keyword{ manip }% at least one, from doc/KEYWORDS
