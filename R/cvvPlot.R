#- $Id: cvvPlot.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Plot the CV result as a function of mean signal intensity.

"cvvPlot" <-
function(data, name = "") {
	sd <- apply(data, 1, function(x) {sd(x, na.rm = T)})
	mean <- rowMeans(data, na.rm = T)
	cvValue <- sd / abs(mean)

	#- Let's create cv plot here
	title = paste("CV plot ", name, sep = "")

	y.max <- max(cvValue, na.rm = T)

	#- Calculate position for text information
	txt <- c(paste("Mean  = ", format(mean(cvValue, na.rm = T) * 100, digits = 2), "%", sep = ""), 
            paste("Median= ", format(median(cvValue, na.rm = T) * 100, digits = 2), "%", sep=""),
            paste("Max   = ", format(max(cvValue, na.rm = T) * 100, digits = 4), "%", sep = ""))

	#- Check to see if mean is in log scale
	if(mean(mean, na.rm = T) > 100) {
		mean <- log2(mean)
	}

	plot(mean, cvValue, col = "blue", pch = ".", main = title, 
			xlab = "mean signal intensity (log2)", ylab = "CV")
	abline(v = 10, col = "green")
	
	#- Draw a line at 25%, 50%, and 75% quantile cv
	cv.q <- quantile(cvValue, na.rm = T)
	abline(h = cv.q[3], col = "red")
	abline(h = cv.q[4], col = "purple")

   legend("topright", txt)
   legend("right", lwd = 2, col = c("red", "purple"), legend = c("50%tile", "75%tile"))
}

############################################################
#- $Log: cvvPlot.R,v $
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
