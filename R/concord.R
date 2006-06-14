#- $Id: concord.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Calculate detection concordance.

"concord" <-
function(sn, snThresh = 3) {
	ncol <- ncol(sn)
	nrow <- nrow(sn)
	data <- array(1, dim=c(ncol, ncol))
	rownames(data) <- colnames(sn)
	colnames(data) <- colnames(sn)

	for(i in 1:ncol) {
	   for(j in 1: ncol) {
			sn.com <- length(intersect(which(sn[,i] >= snThresh), which(sn[,j] >= snThresh)))
			#sn.union <- length(union(which(sn[,i] >= 3), which(sn[,j] >= 3)))
			sn.ave <- (sum(sn[,i] >= snThresh, na.rm = T) + sum(sn[,j] >= snThresh, na.rm = T)) / 2
			#data[i,j] <- sn.com / sn.union
			data[i,j] <- sn.com / sn.ave
		}
	}
	return(data)
}

####################################################################
#- $Log: concord.R,v $
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
