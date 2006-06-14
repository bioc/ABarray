#- $Id: doANOVA.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Perform anova analysis for a given group
doANOVA <- function(dataf, group1, group2, snThresh = 3, detectSample = 0.5) {
  anovaFun <- c()
  if(missing(group2)) {
    cat("\nPerforming one-way ANOVA for", group1, "...\n")
    flush.console()
    resultDir = paste("Result_", gsub(" ", "", group1), "/", sep = "")
    if(! file.exists(resultDir)) {
      dir.create(resultDir, showWarnings = FALSE)
    }
    dataDir = paste(resultDir, "DataResult/", sep = "")
    if(! file.exists(dataDir)) {
      dir.create(dataDir, showWarnings = FALSE)
    }
    fname <- paste(dataDir, "ST_ANOVA_oneway_", group1, ".csv", sep = "")
    if(is(dataf, "exprSet")) {
		if(is.character(group1) & length(group1) < 2) {
        factor1 <- as.factor(as.vector(pData(dataf)[, group1]))
		}
      snPresent = TRUE
      if(dim(se.exprs(dataf))[1] < 2 || sum(is.na(se.exprs(dataf[,1]))) > 1) {
        sn <- array(100, dim(data))
        snPresent <- FALSE
      }
      if(snPresent) {
        snT = snSummary(dataf, snThresh, group1)
        rowUse = which(snT[, 1] >= detectSample)
        for(i in 2:length(levels(factor1))) {
          rowUse = union(rowUse, which(snT[, i] >= detectSample))
        }
        dataf = exprs(dataf[rowUse,])
        cat("    Probes with S/N < ", snThresh, " in at least ", (100 - detectSample * 100),
            "% of samples in all subgroups of ", group1, " were not considered.\n", sep = "")
        flush.console()
      }
      else {
        dataf <- exprs(dataf)
      }
    }
    else {
      factor1 <- as.factor(group1) 
    }
    
    anovaFun <- function(i) {
		return(summary( 
                     aov(dataf[i,] ~ factor1) 
                     ) 
             [[1]][4:5][[2]][1]) 
    }
    anova <- sapply(1:length(dataf[,1]), anovaFun)
    names(anova) <- rownames(dataf)
    
    write.table(data.frame("ProbeID" = rownames(dataf), "p value" = anova), file = fname, sep = ",",
                row.names = F, col.names = T)
    print(paste("One-way ANOVA results (",  length(rownames(dataf)), " probes) ",
                "were written to file: ", fname, sep = ""))    
    invisible(anova)
  }
  else {
    factor1 <- group1
    factor2 <- group2
    cat("Performing two-way ANOVA for", group1, "and", group2, "...\n")
    flush.console()
    resultDir = paste("Result_", gsub(" ", "", group1), "_", gsub(" ", "", group2), "/", sep = "")
    if(! file.exists(resultDir)) {
      dir.create(resultDir, showWarnings = FALSE)
    }

    if(is(dataf, "exprSet")) {
		if(is.character(group1) & length(group1) < 2) {
        group1 <- pData(dataf)[, group1]
		}
      if(is.character(group2) & length(group2) < 2) {
        group2 <- pData(dataf)[, group2]
      }
      snPresent = TRUE
      if(dim(se.exprs(dataf))[1] < 2 || sum(is.na(se.exprs(dataf[,1]))) > 1) {
        snPresent <- FALSE
      }
      if(snPresent) {
        snT = snSummary(dataf, snThresh, factor1)
        rowUse = which(snT[, 1] >= detectSample)
        for(i in 2:length(levels(as.factor(group1)))) {
          rowUse = union(rowUse, which(snT[, i] >= detectSample))
        }
        snT = snSummary(dataf, snThresh, factor2)
        for(i in 1:length(levels(as.factor(group2)))) {
          rowUse = union(rowUse, which(snT[, i] >= detectSample))
        }
        dataf = exprs(dataf[rowUse,])
        cat("    Probes with S/N < ", snThresh, " in at least ", (100 - detectSample * 100),
            "% of samples in all subgroups of ", factor1, " and ", factor2,
            " were not considered.\n", sep = "")
        flush.console()
      }
      else {
        dataf <- exprs(dataf)
      }
    }
    else {
      factor1 <- "factor1"
      factor2 <- "factor2"
    }
    group1 <- as.factor(group1)
    group2 <- as.factor(group2)
        
    anovaFun <- function(i) {
		return(summary( 
                     aov(dataf[i,] ~ group1 * group2) 
                     )
             [[1]][5][[1]][1:3]) 
    }

    anova <- t(sapply(1:length(dataf[,1]), anovaFun))
                                        #names(anova) <- rownames(dataf)
    rownames(anova) <- rownames(dataf)
    colnames(anova) <- c(paste("p(", factor1, ")", sep = ""),
                         paste("p(", factor2, ")", sep = ""),
                         paste("p(", factor1, "*", factor2, ")", sep = ""))

    fname <- paste(resultDir, paste("ST_ANOVA", factor1, factor2, sep = "_"), ".csv", sep = "")
    write.table(cbind(ProbeID = rownames(dataf), anova), file = fname, sep = ",",
                row.names = F, col.names = T)
    print(paste("Two-way ANOVA results (", length(rownames(anova)), " probes) ",
                "were written to file: ", fname, sep = ""))
    invisible(anova)
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
