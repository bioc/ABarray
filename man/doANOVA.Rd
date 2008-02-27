\name{doANOVA}
\alias{doANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform one way or two way ANOVA}
\description{
  If only one factor is provided in parameter, one way ANOVA is
  performed. If two factors are provided, two way ANOVA is performed.
}
\usage{
doANOVA(eset, group1, group2, snThresh = 3, detectSample = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{An \code{ExpressionSet} object.}
  \item{group1}{A factor name or labels to test on. If eset is an
	 \code{ExpressionSet} object, either name or labels can be used.
	 If \code{eset} is an expression matrix, labels should be used.}
  \item{group2}{A factor name or labels to test on.}
  \item{snThresh}{Using probes detectable for ANOVA analysis, default
	 S/N value is 3 or more to be considered detectable.}
  \item{detectSample}{The percentage of samples the probe is detected
	 in order to be considered in ANOVA analysis.}
}
\details{
  At least one group should be provided. If \code{ExpressionSet} object
  is used, group1 or group2 is the name of the sampleGroup defined in
  experiment design file. If labels are to be used, they can be either
  numeric or text, e.g., c(1,1,2,2,3,3) or c("treat1", "treat1", "treat2",
  "treat2", "treat3", "treat3").

  If the probe is detectable in 50\% (default) or more samples in any
  one of the subgroup, it is included in the ANOVA analysis. 
}
\value{
  a vector if one way ANOVA; a matrix if two way ANOVA
}

\author{Y Andrew Sun}

\examples{
  #- one way ANOVA
#-  anova <- doANOVA(eset, "sampleGroup")

  #- two way ANOVA
#-  anova <- doANOVA(eset, "sampleGroup1", "sampleGoup2")   
}

\keyword{htest}% at least one, from doc/KEYWORDS
