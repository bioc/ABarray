\name{imputeFlag}
\alias{imputeFlag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform imputation for missing values (FLAG > 5000) }
\description{
  Perform imputation for missing values.
}
\usage{
imputeFlag(raw, pd = NULL, group = "", impute = "avg")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{ a matrix containing gene expression with missing values
	 labeled as NA }
  \item{pd}{ phenoData object }
  \item{group}{ which group should average be performed }
  \item{impute}{ choice of impute method, only avg (average) is implemented }
}
\value{
  a list containing a matrix with the imputed values and rows that are
  imputed.
  
}

\author{ Y Andrew Sun }

\examples{
#-imputed <- imputeFlag(raw, pd, group = "tissue", impute = "avg")  ##- sn ratio matrix
}
\keyword{ manip }% at least one, from doc/KEYWORDS
