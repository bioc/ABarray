\name{ABarray}
\alias{ABarray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utility to perform QA, data transformation and statistical analysis }
\description{
  (1) Read output from AB1700 software output; (2) Create raw data QA
  and associated plots including boxplot, control data signal plot;
  (3) Missing value calculation; (4) Create MA, scatter plot; (5)
  Perform quantile normalization; (6) Perform t test and fold change,
  or ANOVA (using separate function if more than 2 subgroups). (7)
  Create heatmap with hierarchical clustering. (8) The results are
  either in graphcis or text files.
}
\usage{
ABarray(dataFile, designFile, group, test = TRUE, impute = "avg", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataFile}{ csv or tab delimit file contain expression
	 measurement that are output from AB1700 software
  }
  \item{designFile}{ Experiment design file, including information for
	 sample type and additional phenotype information }
  \item{group}{ specify which group statistical test will be performed
	 on. The samples will be ordered according the group }
  \item{test} { specify whether to perform t test. By default, t test
	 will be performed using specified group information }
  \item{impute} {Treat flagged value (above 5000) as missing value, and
	 impute the missing value }
  \item{...}{ Additional arguments. Use snThresh and/or detectSample to
	 perform filtering. snThresh is the threshold of S/N value to be
	 considered that the probe is detected (default value = 3, if
	 snThresh is not specified). detectSample is used to determine if a
	 probe should be included in statistical analysis (default value =
	 0.5, ie 50\% of samples in any one subgroup). }
}

\details{
  The function works on AB1700 software export data file. It expects
  certain file format to work. The rows of the file represent
  probes. The columns should contain these headings:
  probeID, geneID, Signal, S/N, Flag, and optionally SDEV, CV,
  AssayNormSignal (these values will be ignored in the process).

  It is optional to have control probes. If they are present, plots will
  be generated for the control probes and they will be removed for
  further analysis.

  It is required to have an experiment design file in certain
  format. The rows of the file are samples or arrays. The first column
  should be sampleName. Perhaps, sampleName should be concise and no
  spaces between characters. Second and third columns maybe assayName
  and arrayName (arrayName is optional). Additional columns should
  specify what type of samples. Note: It is best to have assayName
  the same as in dataFile.

  Group name should be the same as in designFile. The samples will be
  ordered according the group information. The samples within the same
  subgroup will be ordered together. Only one group is accepted.

  If test is TRUE (default), t test and ANOVA (if applicable) results
  will be produced. 

  If impute is knn, the flagged values (flag > 5000) will be
  treated as missing values, and the missing values will be imputed
  using k-nearest neighbor average (using bioconductor package impute)
  as described in Olga Troyanskaya, Michael Cantor, Gavin Sherlock,
  Pat Brown, Trevor Hastie, Robert Tibshirani, David Botstein and Russ
  B. Altman, Missing value estimation methods for DNA microarrays
  BIOINFORMATICS Vol. 17 no. 6, 2001 Pages 520-525.

  If impute is avg (default), the signal values of the flagged probes
  will be imputed from average of the subgroup only if there are 2 or
  more values remaining in the subgroup.

  Even if snThresh is not specified in the argument, snThresh is set to
  3 by default. If a value other than 3 is desired (e.g., 2), put
  'snThresh = 2' in the argument.

  detectSample is also preset to a value = 0.5. This means that if a
  probe is detected in 50\% or more samples in any subgroup within the
  group, it is included in statistical analysis. For example, if the
  group is named 'tissue', and there are 2 subgroups named 'lung' and
  'liver', then, if a probe is detected in 50\% or more samples in
  'lung', it is included in the statistical analysis regardless the
  detectability in the other subgroup ('liver').
}

\value{
  An \code{exprSet} object. The \code{exprs} slot will be populated with
  normalized signals, \code{se.exprs} slot will be populated with
  S/N ratio values, and the \code{phenoData} slot will be populated with
  information from \code{designFile}. Further analysis can be performed
  on the \code{exprSet} object with various R and Bioconductor packages.
}

\author{ Yongming Sun <sunya@appliedbiosystems.com> }

\seealso{ doPlotEset, doPlotFCT, doANOVA, matrixPlot, mvaPair2, doLPE,
  doVennDiagram, hclusterPlot }

\examples{
#- eset <- ABarray(dataFile, designFile, "sampleGroup")
#- eset <- ABarray(dataFile, designFile, "group", detectSample = 0.8)
}
\keyword{ manip }% at least one, from doc/KEYWORDS

