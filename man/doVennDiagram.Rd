\name{doVennDiagram}
\alias{doVennDiagram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create Venn Diagram }
\description{
  Create Venn diagram from lists. 
}
\usage{
doVennDiagram(a, b, c = NULL, names, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ a vector of first list }
  \item{b}{ a vector of second list }
  \item{c}{ a vector of third list, optional }
  \item{names}{ a vector for the name of the set }
  \item{...}{ additional graphical parameter }
}
\details{
  The funciton will create Venn diagram. If two lists (a and b) are
  provided, two-way Venn diagram will produced. If three lists (a, b,
  and c) are provided, three-way Venn diagram will be produced.

  This function depends on some functions of limma package, and is
  derived from limma package.
}
\value{
  A plot of Venn diagram
}
\references{ Bioconductor limma package. }
\author{ Yongming Sun }

\examples{
%-doVennDiagram(list1, list2)
%-doVennDiagram(list1, list2, list3, names = c("method1", "method2", "method3"))
}
\keyword{ hplot }% at least one, from doc/KEYWORDS

