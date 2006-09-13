\name{calcsn}
\alias{calcsn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate SN summary for each group }
\description{
  Calculate S/N ratio summary for each group
}
\usage{
calcsn(sn, snThresh, pdata, group, grpMember)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sn}{ S/N ratio data }
  \item{snThresh}{ S/N threshold filtering }
  \item{pdata}{ experiment design }
  \item{group}{ which group should be calculated }
  \item{grpMember}{ optional, members of the group }
}

\value{
  data matrix
}

\author{ Y Andrew Sun }

\examples{
%-snSum <- calcsn(sn, 3, pdata, "sampleGroup", c("liver", "muscle"))
}
\keyword{ manip }% at least one, from doc/KEYWORDS

