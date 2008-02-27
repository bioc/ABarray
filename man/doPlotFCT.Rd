\name{doPlotFCT}
\alias{doPlotFCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate fold change and t test, the plot }
\description{
  Calculate fold changes and p values from t test, and plot the results
  using preset FDR threshold
}
\usage{
doPlotFCT(eset, group, grpMember, order1 = NULL, order2 = NULL,
detectSample = 0.5, snThresh = 3, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{an \code{ExpressionSet} object}
  \item{group}{which group from experiment design should calculation
	 and plot be performed}
  \item{grpMember}{optional group member within the group}
  \item{order1}{optional, For a pairwise comparison the ordering of
	 the first group of replicates}
  \item{order2}{optional, For a pairwise comparison the ordering of
	 the first group of replicates}
  \item{detectSample}{optional number between 0 and 1 to indicate the
	 percentage of arrays should be above snThresh to include in the t
	 test analysis. Default = 0.5. If the probe is detected in 50\% or
	 more samples on one of the subgroup, the probe is included in the
	 t test, otherwise, it will be excluded in the t test}
  \item{snThresh}{ optional S/N ratio threshold. Default = 3}
  \item{\dots}{Additional argument, currently not implemented}
}
\details{
  Group members are optional. For example, if group name is "tissue",
  and group members in experiment design file include "brain", "liver",
  "lung", "muscle". We could include c("brain", "liver") as group member
  for the parameter, then t test will be performed between "brain" and
  "liver", and "lung" "muscle" will be ignored. However, if we omit
  group member in the arguments, all tissue members will be used for t
  test. In this case, there will be 6 pairwise t test between each
  member of the group.

  If order1 and order2 are specified then a paired sample t-test will be
  conducted between the groups, with the arrays in each group sorted
  according to the ordering specified. For example, if order1 is
  c(1,3,2) and order2 is c(1,2,3), then the sample pairing is a1-b1,
  a3-b2, a2-b3, with a and b are subgroup 1 and subgroup 2 within the
  group.

  The fold changes are difference between averaged subgroup1 expression
  vs averaged subgroup2. If paired t test is performed, the fold changes
  are calculated using each paired difference and take an average of
  paired difference.
}
\value{
  None. But a number of plot and result files will be produced.
}
\author{Y Andrew Sun}
\examples{
#- doPlotFCT(eset, "sampleGroup", c("liver", "muscle"))
#- For a paired t test
#- doPlotFCT(eset, "sampleGroup", c("liver", "muscle"), order1 = c(1,2,3), order2 = c(1,3,2))
}
\keyword{ manip }% at least one, from doc/KEYWORDS
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ htest }% at least one, from doc/KEYWORDS
