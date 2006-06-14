\name{getMemberEset}
\alias{getMemberEset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce a sub exprSet given a group and its members }
\description{
  From a group and its member name,  return an exprSet containing just
  these members
}
\usage{
getMemberEset(eset, group, member)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ an \code{exprSet} object }
  \item{group}{ the name of the group which must be in the experiment
	 design file }
  \item{member}{ member name(s) in the above mentioned group }
}

\value{
  an exprSet object
}

\author{ Yongming Sun }

\examples{
%- eset.new <- getMemberEset(eset, "sampleGoup", c("liver", "muscle"))
}
\keyword{ manip }% at least one, from doc/KEYWORDS

