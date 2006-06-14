#- $Id: matrixPlot.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Produce heatmap with a given matrix

"matrixPlot" <-
function (x, nrgcols = 50, rlabels = TRUE, clabels = TRUE, 
    rcols = 1, ccols = 1, k = 10, title = "", ...) 
{
	# layout(matrix(c(1,2), 1, 2), 5, 1)
	 layout(matrix(c(1,1,1,1,1,1,1,2), 2,8, byrow=TRUE))
	
	 color <- rgcolorsfunc(nrgcols)
    n <- nrow(x)
    p <- ncol(x)

	 rlabels <- rownames(x)
	 clabels <- colnames(x)
	 #- Check the largest character length for the label and adjust for margins
	 titleHeight = 4
	 rMax <- max(nchar(rlabels)) + 1
	 cMax <- max(nchar(clabels)) + titleHeight + 2
	 if (cMax > 18) {
		cMax = 18
	 }
	 if (rMax > 15) {
		rMax = 15
	 }
	 
	 op <- par(no.readonly = TRUE) # the whole list of settable par's.
	 par(mar = c(2, rMax, cMax, 2))

    image(1:p, 1:n, t(x[n:1, ]), col = color, 
        axes = FALSE, xlab = "", ylab = "", ...)
    if (length(ccols) == 1) {
        axis(3, at = 1:p, labels = clabels, las = 2, cex.axis = 1, 
            col.axis = ccols)
    }
    if (length(ccols) == p) {
        cols <- unique(ccols)
        for (i in 1:length(cols)) {
            which <- (1:p)[ccols == cols[i]]
            axis(3, at = which, labels = clabels[which], las = 2, 
                cex.axis = 0.6, col.axis = cols[i])
        }
    }
    if (length(rcols) == 1) {
        axis(2, at = n:1, labels = rlabels, las = 2, cex.axis = 1, 
            col.axis = rcols)
    }
    if (length(rcols) == n) {
        cols <- unique(rcols)
        for (i in 1:length(cols)) {
            which <- (1:n)[rcols == cols[i]]
            axis(2, at = (n:1)[which], labels = rlabels[which], 
                las = 2, cex.axis = 0.6, col.axis = cols[i])
        }
    }
    mtext(title, side = 3, line = cMax - titleHeight, cex = 2)
    box()

	 scaleColorBar(x, col = color, k = k, cLen = cMax)

	 #- return to default margins
	 par(mar = c(5, 4, 4, 2) + 0.1)
	 par(op)  #- reset par to oldpar
}

##########################################################
#- $Log: matrixPlot.R,v $
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
