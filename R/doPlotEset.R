#- $Id: doPlotEset.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Main function to perform QC and analysis on exprSet data.

"doPlotEset" <-
function(eset, group, name = "", snThresh = 3, test = TRUE, ...) {
  subgrpCount <- 1
  members <- 1
  color <- rep(c("blue", "cyan", "green", "lightblue", "purple", "pink", "brown", "magenta", "orange"), 20)
  color.label <- c()
  pd <- pData(eset)

  kWidth <- 1024  #- default figure size
  kHeight <- 768

  snPresent <- TRUE
  showHiddenPlot = FALSE
  rawData = FALSE
  scatter <- TRUE
  detectSample = 0.5
  
     #- Process additional arguments for sweave and hidden variables, ...
  addParam = list(...)
  if(is.element("showHiddenPlot", names(addParam))) {
    showHiddenPlot = addParam$showHiddenPlot
  }
  if(is.element("rawData", names(addParam))) {
    rawData = addParam$rawData
  }
  if(is.element("scatter", names(addParam))) {
    scatter = addParam$scatter
  }
  if(is.element("detectSample", names(addParam))) {
    detectSample = addParam$detectSample
  }
  
  arrayCount = dim(exprs(eset))[2]
  arrayPairCount = (arrayCount * (arrayCount - 1)) / 2
  pixelPerArray = as.integer(300 / log2(arrayCount))
  
  data <- exprs(eset)
  sn <- se.exprs(eset)
  
  #- If there is no S/N data, use all values. fill sn with 100 for each value
  if(dim(se.exprs(eset))[1] < 2 || sum(is.na(se.exprs(eset[,1]))) > 1) {
    sn <- array(100, dim(data))
    snPresent <- FALSE
  }
  
  #- reorder the column based on group. Member of the same group should stay together
  idx.order <- c()
  
  if ( ! is.null(pd) ) {
    if( ! any(colnames(pd) == group)) {
      print("The group you provided does not match any in the eset")
      print("No processing will be provided")
      return()
    }
    members <- levels(pd[, colnames(pd) == group])	
    subgrpCount <- length(members)
    
    idx.memberCount <- 1
    for(i in 1:subgrpCount) {
      memberCount <- length(which(pd[, colnames(pd) == group] == members[i]))
      color.label <- c(color.label, rep(color[i], memberCount))
      idx.mc <- which(pd[, colnames(pd) == group] == members[i])
      idx.order <- c(idx.order, idx.mc)
    }
    pd <- pd[idx.order,]
    pData(eset) <- pd
    data <- data[, idx.order]
    sn <- sn[, idx.order]
  }

  if(max(data, na.rm = T) > 2000) {
    data <- log2(data)
  }
	
  #- Create a directory to put all the figures in
	resultDir <- paste("Result_", gsub(" ", "", group), "/", sep = "")
	if(! file.exists(resultDir)) {
		dir.create(resultDir, showWarnings = FALSE)
	}
   dataDir = paste(resultDir, "DataResult/", sep = "")
   pdfDir = paste(resultDir, "pdfFigure/", sep = "")
   jpgDir = paste(resultDir, "jpgFigure/", sep = "")
   if(! file.exists(dataDir)) {
     dir.create(dataDir, showWarning = FALSE)
   }
   if(! file.exists(pdfDir)) {
     dir.create(pdfDir, showWarning = FALSE)
   }
   if(! file.exists(jpgDir)) {
     dir.create(jpgDir, showWarning = FALSE)
   }
   
  
  prefix <- paste(gsub(" ", "", name), "_", gsub(" ", "", group), sep = "")
  
	#- Create boxplot
  print(paste("Creating boxplot for", prefix, " ..."))
  flush.console()
  fname <- paste("QC_SignalBoxplot_", prefix, sep = "")
  #bmp(filename = paste(fname, ".bmp", sep = ""), width = kWidth, height = kHeight)
  pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = kWidth / 100, height = kHeight / 100)
  dev.control("enable")
  chLength = max(nchar(colnames(data)) * 1.2, 5)
  par(mar = c(chLength, 5, 4, 2))
  boxplot(split(data, col(data)), names = colnames(data), col = color.label, las = 2,
          main = paste("Box plot", gsub("_", " ", prefix)))
  savejpg(paste(jpgDir, fname, ".jpg", sep = ""), width = kWidth, height = kHeight)
  dev.off()
   
	#- MA plot, only plot if not raw unprocessed data
  if(! rawData) {
    if(arrayCount > 35 || arrayCount < 2) {
      print(paste(". Will not create MA plot for all samples: ", arrayPairCount, " pairs", sep = ""))
      flush.console()
    }
    else {
      print(paste("Creating MA plot for", prefix))
      flush.console()
		#- MA plot
      fname <- paste("QC_MA_Plot_all_", prefix, sep = "")
      fWidth = arrayCount * pixelPerArray
      fHeight = arrayCount * pixelPerArray
      if(fWidth < kWidth) {fWidth = kWidth}
      if(fHeight < kHeight) {fHeight = kHeight}
      bmp(filename = paste(jpgDir, fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
      #dev.control("enable")
      mvaPair2(data, sn, snThresh = snThresh, main = paste("MA plot", gsub("_", " ", prefix)))
      dev.off()
    }

	#- MA plot for subgroup
    if (subgrpCount > 1) {
      for(i in 1:subgrpCount) {
        idx <- which(pd[, colnames(pd) == group] == members[i])
        if(length(idx) > 1) {
          title = paste("MA Plot ", prefix, ": ", members[i], sep = "")
          print(paste("Creating ", title, " ...", sep = ""))
          flush.console()
          fname <- paste("QC_MA_Plot_", prefix, "_", members[i], sep = "")
          fWidth = pixelPerArray * length(idx)
          fHeight = pixelPerArray * length(idx)
          if(fWidth < kWidth) {fWidth = kWidth}
          if(fHeight < kHeight) {fHeight = kHeight}
          pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth /100, height = fHeight/100)
          dev.control("enable")
          #bmp(filename = paste(jpgDir, fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
          mvaPair2(data[, idx], sn[, idx], snThresh = snThresh, main = gsub("_", " ", title))
          savejpg(paste(jpgDir, fname, sep = ""), width = fWidth, height = fHeight)
          dev.off()
         
			#- Let's create cv plot here too
          fname <- paste("QC_CVplot_", prefix, "_", members[i], sep = "")
          print(paste("Creating ", fname, " ...", sep = ""))
          flush.console()
          fWidth = kWidth * 0.8; fHeight = kHeight * 0.8
          pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth/100, height = fHeight/100)
          dev.control("enable")
          cvvPlot(data[, idx], paste(members[i], "_", prefix, sep = ""))
          savejpg(paste(jpgDir, fname, sep = ""), width = fWidth, height = fHeight)
          dev.off()
        }
      }
    }
  }


  if(! rawData & scatter) {
	#- Scatter plot without filtering
    nonFilterScatter = FALSE
    if(nonFilterScatter) {
      print(paste("Creating scatter plot for", prefix, " ..."))
      flush.console()
      fWidth = arrayCount * pixelPerArray
      fHeight = arrayCount * pixelPerArray
      if(fWidth < kWidth) {fWidth = kWidth}
      if(fHeight < kHeight) {fHeight = kHeight}
      fname <- paste("QC_ScatterPlot_", prefix, sep = "")
      #pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth/100, height = fHeight/100)
      bmp(filename = paste(fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
      #dev.control("enable")
      pairs(data, lower.panel = panel.scatter, upper.panel= panel.cor,
          main = paste("Scatter plot", gsub("_", " ", prefix)))
      #savejpg(paste(jpgDir, fname, sep = ""), width = fWidth, height = fHeight)
      dev.off()
    }
    
	#- Scatter plot filtering S/N >= snTresh
    print(paste("Creating scatter plot filtering S/N <", snThresh," for", prefix, " ..."))
    flush.console()
    fWidth = arrayCount * pixelPerArray
    fHeight = arrayCount * pixelPerArray
    if(fWidth < kWidth) {fWidth = kWidth}
    if(fHeight < kHeight) {fHeight = kHeight}
    fname <- paste("QC_ScatterPlot_FilterSN", snThresh, "_", prefix, sep = "")
    #pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth/100, height = fHeight/100)
    bmp(filename = paste(jpgDir, fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
    #dev.control("enable")
    dataFilter <- data
    dataFilter[sn < snThresh] <- NA
    pairs(dataFilter, lower.panel = panel.scatter, upper.panel= panel.cor, 
          main = paste("Scatter plot SN >=", snThresh, gsub("_", " ", prefix)))
    #savejpg(paste(jpgDir, fname, sep = ""), width = fWidth, height = fHeight)
    dev.off()
    dataFilter = NULL
  }

	#- correlation matrix and plot
  print(paste("Creating correlation matrix and plot", prefix, " ..."))
  flush.console()
  data.cor <- cor(data, use = "complete.obs")
  fname <- paste("QC_Signal_Correlation_", prefix, sep = "")
  pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = kWidth/100, height = kHeight/100)
  #bmp(file = paste(fname, ".bmp", sep = ""), width = 800, height = 600)
  dev.control("enable")
  matrixPlot(data.cor, title = paste("Signal Correlation", gsub("_", " ", prefix)))
  savejpg(paste(jpgDir, fname, sep = ""), width = kWidth, height = kHeight)
  dev.off()
  data.cor = NULL
  
	#- S/N concordance
  if(snPresent) {
    print(paste("Creating detection concordance plot", prefix, "..."))
    flush.console()
    fname <- paste("QC_DetectConcordance_", prefix, sep = "")
    pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = kWidth/100, height = kHeight/100)
    #bmp(filename = paste(fname, ".bmp", sep = ""), width = 800, height = 600)
    dev.control("enable")
    matrixPlot(concord(sn, snThresh), title = paste("Detection Concordance (S/N >= ", snThresh, ")", sep = ""))
    savejpg(paste(jpgDir, fname, sep = ""), width = kWidth, height = kHeight)
    dev.off()
  }
  else {
    print(paste("S/N has no value for ", prefix, ". Skip concordance plot ..."))
    flush.console()
  }

  gc()
  
  if(test) {
      #- perform t test (if more than 2 groups, pairwise test)
    doPlotFCT(eset, group, snThresh = snThresh, ...)
    
		#- Get ANOVA and write anova result into file
    if(subgrpCount > 2) {
      doANOVA(eset, group, snThresh = snThresh, detectSample = detectSample)
    }
  }
}

###############################################################################
#- $Log: doPlotEset.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.8  2006/03/21 23:00:17  sunya
#- Improved memory usage.
#- Control probes in the exported control signal value file are now sorted.
#- Improved handling of NA values. The statistical analysis will ignore NA.
#- If there is only one member of a subgroup, this subgroup will not appear
#- in t test.
#-
#- Revision 1.7  2006/03/20 23:00:16  sunya
#- Missing values (NA) are kept even after imputation. But most downstream
#- analysis will remove these NA values.
#- The hierarchical clustering probe list and expression values are now
#- exported into csv file.
#- Changed file naming convention. QC plot will begin QC_, and statistical
#- analysis plot will begin ST_ in the figures.
#- Values for control probes are saved into csv file.
#-
#- Revision 1.6  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-