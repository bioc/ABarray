#- $Id: doLPE.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Perform LPE analysis on the eset data.

"doLPE" <-
function(eset, group, member, name = "", snThresh = 3, detectSample = 0.5) {
	require(LPE)

  cat("Performing LPE analysis for", group, "...\n")
  flush.console()
	pd <- pData(eset)
	grpMember <- pd[, colnames(pd) == group]

	if(missing(member)) {
		member <- levels(grpMember)
	}

	idx.a <- which(is.element(grpMember, member[1]))
	idx.b <- which(is.element(grpMember, member[2]))

	idx.prob <- 1:dim(exprs(eset))[1]

  snDetectDim <- dim(assayDataElement(eset, "snDetect"))
	if(snDetectDim[1] > 1 & sum(is.na(assayDataElement(eset[,1], "snDetect"))) < 1) {
		snT <- snSummary(eset, snThresh = snThresh, group = group, grpMember = member)
		idx1 <- which(snT[,1] >= detectSample)
		idx2 <- which(snT[,2] >= detectSample)
		idx.prob <- union(idx1, idx2)
	}

	var.a <- baseOlig.error(exprs(eset[idx.prob, idx.a]), q = 0.02)
	var.b <- baseOlig.error(exprs(eset[idx.prob, idx.b]), q = 0.02)

	lpe.var <- data.frame(lpe(exprs(eset[idx.prob, idx.b]), exprs(eset[idx.prob, idx.a]),
									var.b, var.a, probe.set.name = geneNames(eset[idx.prob,])))

	lpe.fdr <- lpe.fdr.BH(lpe.var)

	idx.med.diff <- length(idx.a) + length(idx.b) + 5
	idx.rawp <- length(idx.a) + length(idx.b) + 8
	idx.adjp <- length(idx.a) + length(idx.b) + 9

  result <- cbind(rownames(lpe.fdr), lpe.fdr[ ,c(idx.med.diff, idx.rawp, idx.adjp)])
  colnames(result) <- c("ProbeID",
                         paste("Median log2(FC) (", member[2], "/", member[1], ")", sep = ""),
                         "LPE pVal", "LPE FDR(BH)")
  resultDir <- paste("Result_", gsub(" ", "", group), "/", sep = "")
  if(! file.exists(resultDir)) {
     dir.create(resultDir, showWarnings = FALSE)
	}
  dataDir = paste(resultDir, "DataResult/", sep = "")
  if(! file.exists(dataDir)) {
     dir.create(dataDir, showWarning = FALSE)
  }   
   
  fname <- paste(dataDir, "LPE", name, "_", group, "_", member[2], "_", member[1], ".csv", sep = "")
  write.table(result, fname, sep = ",", col.names = T, row.names = F)
  print(paste("LPE results were written to file", fname))
	invisible(lpe.fdr)
}

#################################################
#- $Log: doLPE.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.3  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
