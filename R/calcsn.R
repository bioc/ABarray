#- $Id: calcsn.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

"calcsn" <-
function(sn, snThresh, pdata, group, grpMember) {
	pd <- pData(pdata)
	if(missing(grpMember)) {
		grpMember <- levels(pd[, colnames(pd) == group])
	}
	groupCount <- length(grpMember)

	snT <- matrix(nrow = dim(sn)[1], ncol = groupCount) 
	colnames(snT) <- grpMember
	rownames(snT) <- rownames(sn)

	for(i in 1:groupCount) {
		memberCount <- length(which(pd[, colnames(pd) == group] == grpMember[i]))
		idx.mc <- which(pd[, colnames(pd) == group] == grpMember[i])

		snSum <- apply(sn[, idx.mc] >= snThresh, 1, sum, na.rm = T)

		snT[, i] <- snSum / memberCount 
	}

	return(snT)
}

#####################################
#- $Log: calcsn.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.3  2006/03/21 23:00:17  sunya
#- Improved memory usage.
#- Control probes in the exported control signal value file are now sorted.
#- Improved handling of NA values. The statistical analysis will ignore NA.
#- If there is only one member of a subgroup, this subgroup will not appear
#- in t test.
#-
#- Revision 1.2  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
