#- $Id: doANOVA.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Perform anova analysis for a given group
doANOVA <- function(eset, group1, group2, snThresh = 3, detectSample = 0.5) {
  anovaFun <- c()
  if(!is(eset, "ExpressionSet")) {
    print("ExpressionSet object must be provided")
    return()
  }
  dataf <- exprs(eset)
  if(missing(group2)) {
    cat("\nPerforming one-way ANOVA for", group1, "...\n")
    flush.console()
    resultDir <- paste("Result_", gsub(" ", "", group1), "/", sep = "")
    if(!file.exists(resultDir)) {
      dir.create(resultDir, showWarnings = FALSE)
    }
    dataDir <- paste(resultDir, "DataResult/", sep = "")
    if(!file.exists(dataDir)) {
      dir.create(dataDir, showWarnings = FALSE)
    }
    fname <- paste(dataDir, "ST_ANOVA_oneway_", group1, ".csv", sep = "")

    factor1 <- as.factor(pData(eset)[, group1])

    snPresent <- TRUE
    if(dim(assayDataElement(eset, "snDetect"))[1] < 2 || 
            sum(is.na(assayDataElement(eset, "snDetect")[,1])) > 1) {
      sn <- array(100, dim(data))
      snPresent <- FALSE
    }
    if(snPresent) {
      snT <- snSummary(eset, snThresh, group1)
      rowUse <- which(snT[, 1] >= detectSample)
      for(i in 2:length(levels(factor1))) {
        rowUse <- union(rowUse, which(snT[, i] >= detectSample))
      }
      dataf <- exprs(eset[rowUse,])
      cat("    Probes with S/N < ", snThresh, " in at least ", (100 - detectSample * 100),
          "% of samples in all subgroups of ", group1, " were not considered.\n", sep = "")
      flush.console()
    }
    
    anovaFun <- function(i) {
      return(summary(aov(dataf[i,] ~ factor1))[[1]][4:5][[2]][1]) 
    }
    anova <- sapply(1:length(dataf[,1]), anovaFun)
    names(anova) <- rownames(dataf)
    
    write.table(data.frame("ProbeID" = rownames(dataf), "p value" = anova), file = fname, sep = ",",
                row.names = FALSE, col.names = TRUE)
    print(paste("One-way ANOVA results (",  length(rownames(dataf)), " probes) ",
                "were written to file: ", fname, sep = ""))    
    invisible(anova)
  }
  else {
    cat("Performing two-way ANOVA for", group1, "and", group2, "...\n")
    flush.console()
    resultDir <- paste("Result_", gsub(" ", "", group1), "_", gsub(" ", "", group2), "/", sep = "")
    if(!file.exists(resultDir)) {
      dir.create(resultDir, showWarnings = FALSE)
    }
    factor1 <- as.factor(pData(eset)[, group1])
    factor2 <- as.factor(pData(eset)[, group2])
    if(is(eset, "ExpressionSet")) {
      snPresent <- TRUE
      if(dim(assayDataElement(eset, "snDetect"))[1] < 2 || 
            sum(is.na(assayDataElement(eset, "snDetect")[,1])) > 1) {
        sn <- array(100, dim(data))
        snPresent <- FALSE
      }
      if(snPresent) {
        snT <- snSummary(eset, snThresh, group1)
        rowUse <- which(snT[, 1] >= detectSample)
        for(i in 2:length(levels(factor1))) {
          rowUse <- union(rowUse, which(snT[, i] >= detectSample))
        }
        snT <- snSummary(eset, snThresh, group2)
        for(i in 1:length(levels(factor2))) {
          rowUse <- union(rowUse, which(snT[, i] >= detectSample))
        }
        dataf <- exprs(eset[rowUse,])
        cat("    Probes with S/N < ", snThresh, " in at least ", (100 - detectSample * 100),
            "% of samples in all subgroups of ", group1, " and ", group2,
            " were not considered.\n", sep = "")
        print(paste("   Number of probes considered:", length(rowUse)))
        flush.console()
      }
    }
        
    anovaFun <- function(i) {
		  return(summary(aov(dataf[i,] ~ factor1 * factor2))[[1]][5][[1]][1:3]) 
    }

    anv <- t(sapply(1:length(dataf[,1]), anovaFun))
                                        #names(anova) <- rownames(dataf)
    rownames(anv) <- rownames(dataf)
    colnames(anv) <- c(paste("p(", group1, ")", sep = ""),
                         paste("p(", group2, ")", sep = ""),
                         paste("p(", group1, "*", group2, ")", sep = ""))

    fname <- paste(resultDir, paste("ST_ANOVA", group1, group2, sep = "_"), ".csv", sep = "")
    write.table(cbind(ProbeID = rownames(dataf), anv), file = fname, sep = ",",
                row.names = FALSE, col.names = TRUE)
    print(paste("Two-way ANOVA results (", length(rownames(anv)), " probes) ",
                "were written to file: ", fname, sep = ""))
    invisible(anv)
  }
}

#-----------------
#-$Log: doANOVA.R,v $
#-Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#-ABarray project converted from ab1700 project
#-
#-Revision 1.6  2006/03/20 23:00:16  sunya
#-Missing values (NA) are kept even after imputation. But most downstream
#-analysis will remove these NA values.
#-The hierarchical clustering probe list and expression values are now
#-exported into csv file.
#-Changed file naming convention. QC plot will begin QC_, and statistical
#-analysis plot will begin ST_ in the figures.
#-Values for control probes are saved into csv file.
#-
#-Revision 1.5  2006/03/14 19:48:30  sunya
#-Changed icp (internal control probe) QC plots.
#-Added function for icp -> icpPlot
#-ANOVA analysis now performs probe filtering, but no FDR is calculated.
#-hclusterPlot now calculate correlation coefficient for probes, previously
#-it used Euclidian distance. The distance between arrays is still Euclidean.
#-
#-Revision 1.4  2005/10/17 21:05:00  sunya
#-Added ab1700gui document. Changed DESCRIPTION version to 1.0.1
#-
#-Revision 1.3  2005/10/07 23:36:47  sunya
#-If eSet and S/N ratio is present, removing probes that are non-detectable
#-in some percentage of samples in all the factor and subgroups, as these
#-probes are non-detectable across all conditions.
#-
