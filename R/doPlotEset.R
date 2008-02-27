#- $Id: doPlotEset.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Main function to perform QC and analysis on ExpressionSet data.

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
  showHiddenPlot <- FALSE
  rawData <- FALSE
  scatter <- TRUE
  maPlot <- TRUE
  detectSample <- 0.5
  
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
  if(is.element("maPlot", names(addParam))) {
    maPlot = addParam$maPlot
  }
  if(is.element("detectSample", names(addParam))) {
    detectSample = addParam$detectSample
  }
  
  arrayCount = dim(exprs(eset))[2]
  arrayPairCount = (arrayCount * (arrayCount - 1)) / 2
  pixelPerArray = as.integer(300 / log2(arrayCount))
  
  sigData <- assayDataElement(eset, "exprs")
  snDetect <- assayDataElement(eset, "snDetect")
  
  #- If there is no S/N data, use all values. fill snDetect with 100 for each value
  if(dim(snDetect)[1] < 2 || max(rowSums(is.na(snDetect))) > 20000) {
    snDetect <- array(100, dim(sigData))
    snPresent <- FALSE
  }
  
  #- reorder the column based on group. Member of the same group should stay together
  idx.order <- c()
  
  if(!is.null(pd) ) {
    if(!any(colnames(pd) == group)) {
      print("The group you provided does not match any in the eset")
      print("No processing will be provided")
      return()
    }
    members <- levels(pd[, group])	
    subgrpCount <- length(members)
    
    idx.memberCount <- 1
    for(i in 1:subgrpCount) {
      memberCount <- length(which(pd[, group] == members[i]))
      color.label <- c(color.label, rep(color[i], memberCount))
      idx.mc <- which(pd[, group] == members[i])
      idx.order <- c(idx.order, idx.mc)
    }
    pd <- pd[idx.order,]
    pData(eset) <- pd
    sigData <- sigData[, idx.order]
    snDetect <- snDetect[, idx.order]
  }

  if(max(sigData, na.rm = TRUE) > 2000) {
    sigData <- log2(sigData)
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
  chLength = max(nchar(colnames(sigData)) * 1.2, 5)
  par(mar = c(chLength, 5, 4, 2))
  boxplot(split(sigData, col(sigData)), names = colnames(sigData), col = color.label, las = 2,
          main = paste("Box plot", gsub("_", " ", prefix)))
  savejpg(paste(jpgDir, fname, sep = ""), width = kWidth, height = kHeight)
  dev.off()
   
	#- MA plot, only plot if not raw unprocessed data
  if(! rawData && maPlot) {
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
      ##-bmp(filename = paste(jpgDir, fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
      pdf(file=paste(pdfDir, fname, ".pdf", sep=""), width=fWidth/100, height=fHeight/100)
      dev.control("enable")
      mvaPair2(sigData, snDetect, snThresh = snThresh, main = paste("MA plot", gsub("_", " ", prefix)))
      savejpg(paste(jpgDir, fname, sep=""), width=fWidth, height=fHeight)
      dev.off()
    }

	#- MA plot for subgroup
    if (subgrpCount > 1 && maPlot) {
      for(i in 1:subgrpCount) {
        idx <- which(pd[, colnames(pd) == group] == members[i])
        if(length(idx) > 1) {
          title <- paste("MA Plot ", prefix, ": ", members[i], sep = "")
          print(paste("Creating ", title, " ...", sep = ""))
          flush.console()
          fname <- paste("QC_MA_Plot_", prefix, "_", members[i], sep = "")
          fWidth <- pixelPerArray * length(idx)
          fHeight <- pixelPerArray * length(idx)
          if(fWidth < kWidth) {fWidth = kWidth}
          if(fHeight < kHeight) {fHeight = kHeight}
          pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth /100, height = fHeight/100)
          dev.control("enable")
          #bmp(filename = paste(jpgDir, fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
          mvaPair2(sigData[, idx], snDetect[, idx], snThresh = snThresh, main = gsub("_", " ", title))
          savejpg(paste(jpgDir, fname, sep = ""), width = fWidth, height = fHeight)
          dev.off()
         
			#- Let's create cv plot here too
          fname <- paste("QC_CVplot_", prefix, "_", members[i], sep = "")
          print(paste("Creating ", fname, " ...", sep = ""))
          flush.console()
          fWidth <- kWidth * 0.8
		  fHeight <- kHeight * 0.8
          pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth/100, height = fHeight/100)
          dev.control("enable")
          cvvPlot(sigData[, idx], paste(members[i], "_", prefix, sep = ""))
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
      fWidth <- arrayCount * pixelPerArray
      fHeight <- arrayCount * pixelPerArray
      if(fWidth < kWidth) {fWidth = kWidth}
      if(fHeight < kHeight) {fHeight = kHeight}
      fname <- paste("QC_ScatterPlot_", prefix, sep = "")
      #pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth/100, height = fHeight/100)
      bmp(filename = paste(fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
      #dev.control("enable")
      pairs(sigData, lower.panel = panel.scatter, upper.panel= panel.cor,
          main = paste("Scatter plot", gsub("_", " ", prefix)))
      #savejpg(paste(jpgDir, fname, sep = ""), width = fWidth, height = fHeight)
      dev.off()
    }
    
	#- Scatter plot filtering S/N >= snTresh
    print(paste("Creating scatter plot filtering S/N <", snThresh," for", prefix, " ..."))
    flush.console()
    fWidth <- arrayCount * pixelPerArray
    fHeight <- arrayCount * pixelPerArray
    if(fWidth < kWidth) {fWidth = kWidth}
    if(fHeight < kHeight) {fHeight = kHeight}
    fname <- paste("QC_ScatterPlot_FilterSN", snThresh, "_", prefix, sep = "")
    #pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = fWidth/100, height = fHeight/100)
    bmp(filename = paste(jpgDir, fname, ".bmp", sep = ""), width = fWidth, height = fHeight)
    #dev.control("enable")
    dataFilter <- sigData
    dataFilter[snDetect < snThresh] <- NA
    pairs(dataFilter, lower.panel = panel.scatter, upper.panel= panel.cor, 
          main = paste("Scatter plot SN >=", snThresh, gsub("_", " ", prefix)))
    #savejpg(paste(jpgDir, fname, sep = ""), width = fWidth, height = fHeight)
    dev.off()
    dataFilter <- NULL
  }

	#- correlation matrix and plot
  print(paste("Creating correlation matrix and plot", prefix, " ..."))
  flush.console()
  data.cor <- cor(sigData, use = "complete.obs")
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
    matrixPlot(concord(snDetect, snThresh), title = paste("Detection Concordance (S/N >= ", snThresh, ")", sep = ""))
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