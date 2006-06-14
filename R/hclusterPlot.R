#- $Id: hclusterPlot.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Perform hierarchical clustering and produce heatmap on expression matrix

"hclusterPlot" <-
  function(expr, title = "", dist = "Euclidean") {
    if(dim(expr)[1] < 2) {
      print("There is only one probe. No need for clustering")
      return()
    }
    if(dim(expr)[2] < 2) {
		print("Only one sample is available. No need for clustering")
		return()
    }
    
    n <- 255                 #- Use 255
                                        # convert to rgb
    low <- col2rgb("green")
    mid <- col2rgb("black")
    high <- col2rgb("red")
    
                                        # determine length of each component
    lower <- floor(n/2)
    upper <- n - lower
    
    lowcolor  <- c(
                   seq(low[1,1], mid [1,1], length=lower),
                   seq(mid[1,1], high[1,1], length=upper)
                   )/n
    
    highcolor <- c(
                   seq(low[3,1], mid [3,1], length=lower),
                   seq(mid[3,1], high[3,1], length=upper)
                   )/n
    
    midcolor <- c(
                  seq(low[2,1], mid [2,1], length=lower),
                  seq(mid[2,1], high[2,1], length=upper)
                  )/n
    
    
    hcolor <- rgb(lowcolor, midcolor, highcolor)

    if(any(grep("cor", dist, ignore.case = T)) > 0) {
      hv = heatmap(expr, distfun = function(x) (as.dist(1-abs(cor(t(x), use = "complete.obs")))),
          col = hcolor, main = paste(title, dist), margins = c(10,5),
          cex.main = 0.5, cex.axis = 0.6, keep.dendro = T)
    }
    else {
      rv = as.dendrogram(hclust(as.dist(1-cor(t(expr), use = "complete.obs"))))
      hv = heatmap(expr, Rowv = rv, col = hcolor, main = paste(title, dist), margins = c(10,5),
              cex.main = 0.5, cex.axis = 0.6, keep.dendro = T)
    }
    invisible(hv)
  }

###################################################################
#- $Log: hclusterPlot.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.7  2006/03/21 23:00:17  sunya
#- Improved memory usage.
#- Control probes in the exported control signal value file are now sorted.
#- Improved handling of NA values. The statistical analysis will ignore NA.
#- If there is only one member of a subgroup, this subgroup will not appear
#- in t test.
#-
#- Revision 1.6  2006/03/20 23:00:16  sunya
#- Missing values (NA) are kept even after imputation. But most downstream
#- analysis will remove these NA values.
#- The hierarchical clustering probe list and expression values are now
#- exported into csv file.
#- Changed file naming convention. QC plot will begin QC_, and statistical
#- analysis plot will begin ST_ in the figures.
#- Values for control probes are saved into csv file.
#-
#- Revision 1.5  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
