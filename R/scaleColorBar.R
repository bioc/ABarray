#- $Id: scaleColorBar.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Scale bar heatmap for matrixplot.

"scaleColorBar" <-
function (x, horizontal = FALSE, col = rgcolorsfunc(50), scale = 1:length(x), 
    k = 10, cLen = 9, ...) 
{
	 x <- unique(sort(as.vector(x)))
	 scale <- 1:length(x)	

    if (is.numeric(x)) {
        #- x <- x
		  #- Let's use just a sequence from min to max, not just real value of x
		  x <- seq(min(x), max(x), (max(x) - min(x)) / length(x))
        colmap <- col
    }
    else {
        colmap <- x
        low <- range(scale)[1]
        high <- range(scale)[2]
        x <- seq(low, high, length = length(x))
    }

    if (length(x) > k) 
        x.small <- seq(x[1], x[length(x)], length = k)
    else x.small <- x

	 #- Create appropriate margins
	 par(mar = c(2, 2, cLen, 5))

    if (horizontal) {
        image(x, 1, matrix(x, length(x), 1), axes = FALSE, xlab = "", 
            ylab = "", col = colmap, ...)
        axis(1, at = rev(x.small), labels = signif(rev(x.small), 
            2), srt = 270)
    }
    if (!horizontal) {
        image(1, x, matrix(x, 1, length(x)), axes = FALSE, xlab = "", 
            ylab = "", col = colmap, ...)
        par(las = 1)
        #axis(2, at = rev(x.small), labels = signif(rev(x.small), 2))
        axis(4, at = rev(x.small), labels = signif(rev(x.small), 3))
        par(las = 0)
    }
    box()
}

##############################################
#- $Log: scaleColorBar.R,v $
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
