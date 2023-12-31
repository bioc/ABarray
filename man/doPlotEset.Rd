\name{doPlotEset}
\alias{doPlotEset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce a number of QA plot plus t and ANOVA test }
\description{
  Produce boxplot, MA plot, scatter plot, correlation,
  S/N detection concordance, CV, and t test, ANOVA test
  if subgroup is more than 2
}
\usage{
doPlotEset(eset, group, name = "", snThresh = 3, test = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{an \code{ExpressionSet} object}
  \item{group}{name of the group from experiment design file}
  \item{name}{a name for use in output files for record purpose}
  \item{snThresh}{threshold of S/N considered detectable, default = 3}
  \item{test}{whether t or ANOVA test should be performed}
  \item{\dots}{Additional arguments, currently not implemented}
}
\details{
  The t test and fold change is performed with function
  \code{fctPlot}. See additional information with \code{fctPlot}. ANOVA
  is performed with \code{doANOVA}.

  If there are more than 2 subgroup in \code{group}, t test and fold
  change will be performed for each pair of subgroup and one way ANOVA
  will be performed. If subgroup is 2, ANOVA will not be performed.
}
\value{
  None. A number of plots and t or ANOVA test result file will be
  produced.
}
\author{Y Andrew Sun}
\examples{
#-doPlotEset(eset, "sampleGroup")
#-doPlotEset(eset, "sampleGroup", name = "perfect")
#-doPlotEset(eset, "sampleGroup", test = FALSE)  ##- t test will be not performed
}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{hplot}% at least one, from doc/KEYWORDS
