#- $Id: drawVennDiagram.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Drawing the VennDiagram. Part of the source code is modified
#- from related Venn diagram functions in limma package

"drawVennDiagram" <-
  function (object, names, mar = rep(1, 4), cex = 1.5,  ...) 
{
  colors <- c("blue", "cyan", "orange")
  if (!is(object, "VennCounts")) 
    object <- vennCounts(object)
  nsets <- ncol(object) - 1
  if (nsets > 3) 
    stop("Can't plot Venn diagram for more than 3 sets")
  if (missing(names)) 
    names <- colnames(object)[1:nsets]
  counts <- object[, "Counts"]
  totalCounts = sum(counts)
  theta <- 2 * pi * (1:360)/360
  xcentres <- list(0, c(-1, 1), c(-1, 1, 0))[[nsets]]
  #ycentres <- list(0, c(0, 0), c(1/sqrt(3), 1/sqrt(3), -2/sqrt(3)))[[nsets]]
  ycentres <- list(0, c(0, 0), c(1/sqrt(5), 1/sqrt(5), -2/sqrt(5)))[[nsets]]
  r <- c(1.5, 1.5, 1.5)[nsets]
  xtext <- list(-1.2, c(-1.2, 1.2), c(-1.2, 1.2, 0))[[nsets]]
  ytext <- list(1.8, c(1.8, 1.8), c(2.4, 2.4, -3))[[nsets]]
  old.par <- par(mar = mar)
  on.exit(par(old.par))
  plot(x = 0, y = 0, type = "n", xlim = c(-4, 4),
       ylim = c(-4, 4), xlab = "", ylab = "", axes = FALSE, ...)
  for (circle in 1:nsets) {
    lines(xcentres[circle] + r * cos(theta), ycentres[circle] + 
          r * sin(theta), col = colors[circle], lwd = 3)
    text(xtext[circle], ytext[circle], names[circle], cex = cex)
  }
  switch(nsets, {
    rect(-3, -2.5, 3, 2.5)
    text(2.3, -2.1, totalCounts, cex = cex)
    text(0, 0, counts[2], cex = cex)
  }, {
    rect(-3, -2.5, 3, 2.5)
    text(2.3, -2.1, totalCounts, cex = cex)
    text(1.5, 0.1, counts[2], cex = cex)
    text(-1.5, 0.1, counts[3], cex = cex)
    text(0, 0.1, counts[4], cex = cex)
  }, {
    rect(-3, -3.5, 3, 3.3)
    text(2.5, -3, totalCounts, cex = cex)
    text(0, -1.7, counts[2], cex = cex)
    text(1.5, 1, counts[3], cex = cex)
    text(0.75, -0.35, counts[4], cex = cex)
    text(-1.5, 1, counts[5], cex = cex)
    text(-0.75, -0.35, counts[6], cex = cex)
    text(0, 0.9, counts[7], cex = cex)
    text(0, 0, counts[8], cex = cex)
  })
  invisible()
}

######################################################################
#- $Log: drawVennDiagram.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.4  2006/03/27 23:09:04  sunya
#- Added 'no imputation' option in GUI interface. Additional cosmetic changes.
#-
