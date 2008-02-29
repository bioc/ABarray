#- $Id: snSummary.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Calculate S/N summary for each subgroup.

"snSummary" <-
function(eset, snThresh = 3, group, grpMember) {
	pd <- pData(eset)
	if(missing(grpMember)) {
		grpMember <- levels(factor(pd[, colnames(pd) == group]))
	}
	subgrpCount <- length(grpMember)

	sn <- assayDataElement(eset, "snDetect")

	snT <- matrix(nrow = dim(sn)[1], ncol = subgrpCount) 
	colnames(snT) <- grpMember
	rownames(snT) <- rownames(sn)

	for(i in 1:subgrpCount) {
		idx.member <- which(pd[, colnames(pd) == group] == grpMember[i])
		memberCount <- length(idx.member)
		snSum <- apply(sn[, idx.member] >= snThresh, 1, sum, na.rm = TRUE)

		snT[, i] <- snSum / memberCount 
	}

	return(snT)
}

#####################################################
#- $Log: snSummary.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.2  2006/03/14 19:48:31  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
