#- $Id: mvaPair2.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

"mvaPair2" <-
function (x, y = NULL, snThresh = 3, labels = colnames(x), log.it = FALSE, span = 2/3, 
    family.loess = "gaussian", digits = 3, line.col = 2, main = "MA plot",  ...) 
{
  J <- dim(x)[2]
  rowna <- c()
  for(i in seq(J)) {
    rowna <- union(rowna, which(is.na(x[, i])))
    rowna <- union(rowna, which(is.infinite(x[,i])))
  }
  if(length(rowna) > 0) {
    x <- x[-rowna,]
  }
  
  if (log.it) 
    x <- log2(x)

  probeCount <- dim(x)[1]
  cexVal <- 1
  if (J < 8) { cexVal <- 2 }
  else { 
    if (J < 12) { cexVal <- 1}
    else { cexVal <- 0.5 }
  }

    #- Convert S/N ratio data to see if they are above cutoff
  if (!is.null(y) ) {
    if(length(rowna) > 0) {
      y <- y[-rowna,]
    }
    y <- y >= snThresh
  }
  
  frame()
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mfrow = c(J, J), mgp = c(0, 0.2, 0), mar = c(1, 1, 1, 1), oma = c(1, 1.4, 2, 1))
  for (j in 1:(J - 1)) {
    par(mfg = c(j, j))
    plot(1, 1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(1, 1, labels[j], cex = 1)
    for (k in (j + 1):J) {
      par(mfg = c(j, k))
      
      idx <- 1:dim(x)[1]
      if ( ! is.null(y)) {
        idx <- which((y[,j] + y[, k]) > 1)
      }
      yy <- x[, j] - x[, k]
      xx <- (x[, j] + x[, k]) / 2
      
      xx[xx == -Inf] < NA
            
	  ##- Calculate correlation coefficient for signals S/N >= snThresh
      r <- format(cor(x[idx, j], x[idx, k], use = "complete.obs"), digits = digits)
      
      subset <- sample(1:length(x), min(c(10000, length(x))))
      mamaplot(xx, yy, idx, tck = 0, subset = subset, 
               pch = ".",  xlab = "", ylab = "", tck = 0, ylim = c(-3.5, 3.5),
               xlim = c(min(xx, na.rm = TRUE) - 0.5, max(xx, na.rm = TRUE) + 0.5), ...)
      
      abline(1, 0, col = "red")
      abline(-1, 0, col = "red")
      par(mfg = c(k, j))
      
      plot(c(0, 1), c(0, 1), type = "n", ylab = "", xlab = "", 
           xaxt = "n", yaxt = "n")
      snText <- paste("S/N", snThresh, ": ", 
                      format(100 * length(idx) / probeCount, digits = digits), "%", sep = "" )
      
      text(0.5, 0.5, paste(snText, paste("r=", r), sep = "\n"), cex = 1)
      
    }
  }
  par(mfg = c(J, J))
  plot(1, 1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  text(1, 1, labels[J], cex = 1)
  mtext("A", 1, outer = TRUE, cex = 1)
  mtext("M", 2, outer = TRUE, cex = 1, las = 1)
  mtext(main, 3, outer = TRUE, cex = 1.5)
  invisible()
}

########----------------------------------------------------------
#- $Log: mvaPair2.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.4  2005/10/07 18:05:43  sunya
#- 1) Remove manufacturing control probes from output
#- 2) Remove txt files for the control probes
#- 3) Change pVal, FC to volcano plot
#- 4) Bitmapped images (pdf and bmp format output)
#- 5) Header for quantile normalized csv file
#- 6) FC (which vs which)
#- 7) Add probeID for the first column
#- 8) Added progress messages
#- 9) Changed output directories
#- 10)Changed legend for cvplots
#-
#- Revision 1.3  2005/08/03 23:16:42  sunya
#- Remove -Inf (treat as NA) in MA plot.
#-
#- Revision 1.2  2005/07/29 21:45:59  sunya
#- Allow NA in data to plot as NA will be ignored in plot.
#-
