#- $Id: mamaplot.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Perform MA plot

"mamaplot" <-
function (A, M, idx, subset = sample(1:length(M), min(c(10000, length(M)))), 
    span = 2/3, family.loess = "gaussian", cex = 2, ...) 
{
	 #- idx: index of probes that are S/N >= 3
    fn.call <- list(...)
    mean <- median(M, na.rm = TRUE)
    if (!is.element("ylim", names(fn.call))) {
        yloc <- max(M, na.rm = TRUE)
    }
    else {
        yloc <- max(fn.call$ylim, na.rm = TRUE)
    }
    if (!is.element("xlim", names(fn.call))) {
        xloc <- max(A, na.rm = TRUE)
    }
    else {
        yloc <- max(fn.call$xlim, na.rm = TRUE)
    }
    aux <- loess(M[subset] ~ A[subset], degree = 1, span = span, 
        family = family.loess)$fitted

	 plot(A[-idx], M[-idx], col = "lightblue", ...)
	 points(A[idx], M[idx], col = "blue", pch = ".")

    #o <- order(A[subset])
    #A <- A[subset][o]
    #M <- aux[o]
    #o <- which(!duplicated(A))
    #lines(approx(A[o], M[o]), col = "red")
    abline(0, 0, col = "blue")
}

################################################
#- $Log: mamaplot.R,v $
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
