\name{snSummary}
\alias{snSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create summary information for S/N ratio}
\description{
  Create summary information for S/N ratio for each sample group
}
\usage{
snSummary(eset, snThresh = 3, group, grpMember)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{an \code{ExpressionSet} object}
  \item{snThresh}{S/N ratio threshold to use, default = 3}
  \item{group}{sample group}
  \item{grpMember}{sample group members, optional}
}
\value{
  a matrix containing the number of samples with S/N >=3 for each probe
}
\author{Yongming Sun}
\examples{
%-snSum <- snSummary(eset, "sampleGoup")
%-snSum <- snSummary(eset, "sampleGoup", c("liver", "muscle")
}
\keyword{manip}% at least one, from doc/KEYWORDS
