\name{doLPE}
\alias{doLPE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform LPE analysis }
\description{
  The local pooled error test attempts to reduce dependence on the
  within-gene estimates in tests for differential expression, by pooling
  error estimates within regions of similar intensity. Note that with
  the large number of genes there will be genes with low within-gene
  error estimates by chance, so that some signal-to-noise ratios will be
  large regardless of mean expression intensities and fold-change. The
  local pooled error attempts to avert this by combining within-gene
  error estimates with those of genes with similar expression intensity.
}
\usage{
doLPE(eset, group, member, name = "", snThresh = 3, detectSample = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ an exprSet \code{exprSet} object  }
  \item{group}{ which group should LPE be performed }
  \item{member}{ optional. The member names in the group specified above }
  \item{name}{ a prefix name for use when writing output to file }
  \item{snThresh}{ S/N ratio threshold to use to define gene detectability}
  \item{detectSample}{ percetage of samples detectable above snThresh to
	 include in LPE test. The default is 50\%. If the probe is detected
	 in 50\% or more samples in one of the subgroup, it is considered in
	 LPE analysis }
}

\details{
  The LPE test statistic numerator is the difference in medians between
  the two experimental conditions. The test statistic denominator is the
  combined pooled standard error for the two experimental conditions
  obtained by looking up the var.M from each baseOlig.error variance
  function. The conversion to p-values is based on the Gaussian
  distribution for difference if order statistics (medians). The user
  may select bith the smoother degrees of freedom (smaller is smoother)
  and the trim percent to obtain a variance fucntion to suit particular
  issues i.e. variability of genes with low expression intensity.
}

\value{
  Dataframe
}
\references{ Bioconductor LPE package }
\author{ Y Andrew Sun }

\examples{
##---- Some example usage ----

%-result <- doLPE(eset, "tissue", c("liver", "muscle"))
%-result <- doLPE(eset, "tissue", c("liver", "muscle"), detectSample = 0.5)

}
\keyword{ manip }% at least one, from doc/KEYWORDS

