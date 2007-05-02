#- $Id: concord.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Calculate detection concordance.

"concord" <-
function(sn, snThresh = 3) {
	ncol <- ncol(sn)
	nrow <- nrow(sn)
	concor <- array(1, dim=c(ncol, ncol))
	rownames(concor) <- colnames(sn)
	colnames(concor) <- colnames(sn)
	
	snr <- sn >= snThresh
	detect <- colSums(snr, na.rm=T)
	for(i in 1:ncol) {
    for(j in 1:ncol) {
      tmp <- cbind(snr[, i], snr[, j])
      sn.com <- sum(rowSums(tmp, na.rm=T) >= 2)
      sn.ave <- (detect[i] + detect[j]) / 2
      concor[i, j] <- sn.com / sn.ave
    }
  }
  return(concor)
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
