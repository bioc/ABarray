#- $Id: icpPlot.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Create QC plots for control probes

icpPlot = function(controlData, colProbeID = 1, plotWhat = "Signal", pdfDir, jpgDir) {
  controls <- c("Hybridization_Control", "Negative_Control",
    "IVT_Kit_Control_BIOB", "IVT_Kit_Control_BIOC", "IVT_Kit_Control_BIOD",
    "RT_Kit_Control_DAP", "RT_Kit_Control_LYS", "RT_Kit_Control_PHE")
  controlLabel <- gsub("_Control", "", controls)
  controlLabel <- gsub("_Kit", "", controlLabel)

  kWidth <- 800
  kHeight <- 600
  kArrayPerRow <- 20
  if(missing(pdfDir)) {
    pdfDir <- "figTest/"
  }
  if(missing(jpgDir)) {
    jpgDir <- "figTest/"
  }
  
  #- sort probe names
  probeSort <- sort(controlData[, colProbeID], index = TRUE)$ix
  controlData <- controlData[probeSort,]
  
  for(ctl in seq(along = controls)) {
    rowCtl <- grep(controls[ctl], controlData[, colProbeID])
    if(length(rowCtl) > 0) {
      theHeight <- 12 * length(rowCtl)
      if(theHeight < kHeight) {
        theHeight <- kHeight
      }
      plotData <- as.matrix(controlData[rowCtl, -colProbeID])
        
      rownames(plotData) <- gsub("_Control", "", controlData[rowCtl, colProbeID])
      rownames(plotData) <- gsub("_Cp", "", rownames(plotData))
      controlShortNames <- sort(unique(gsub("_\\d+", "", rownames(plotData), perl = TRUE)))     
      
      arrayCount <- dim(plotData)[2]
      rFactor <- arrayCount %/% kArrayPerRow + 1
      eachRow <- arrayCount %/% rFactor + 1
      if(controls[ctl] == "Negative_Control") {
        rFactor <- 1
        eachRow <- arrayCount
      }
      
      fWidth <- 12 * (arrayCount + 5)
      if(arrayCount < 41) {
        fWidth <- 24 * (arrayCount + 5)
      }
      if(fWidth < kWidth) {
        fWidth <- kWidth
      }
      fHeight <- kHeight
      if(rFactor > 2) {
        fHeight <- rFactor * kHeight / 1.5
      }

      fname <- paste("QC_control", plotWhat, "_", controlLabel[ctl], sep = "")
      jpeg(file = paste(jpgDir, fname, "_Heatmap.jpg", sep = ""), width = fWidth, height = fHeight)
                                        #dev.control("enable")                          
      matrixPlot(plotData, k = 20, title = paste(plotWhat, controls[ctl]))             
      dev.off()
      
      print(paste("Creating plot for", plotWhat, controls[ctl], "..."))
      flush.console()

      pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = kWidth / 100, height = fHeight / 100)
      dev.control("enable")

      par(mfrow = c(rFactor, 1), mar = c(5,5,1,1), oma = c(1,0,4,0))
      for(row in seq(rFactor)) {
        rightCount <- min(arrayCount, row * eachRow)
        leftCount <- (row -1) * eachRow + 1        
        for(r in seq(along = controlShortNames)) {
          rowIndIdx <- grep(controlShortNames[r], rownames(plotData))
          thePlotData <- plotData[rowIndIdx,]
          uniqControlName <- sort(unique(rownames(thePlotData)))
          controlCount <- length(uniqControlName)        
          plotData2 <- NULL
          for(i in leftCount:rightCount) {                    
            for(j in seq(controlCount)) {
              rowIdx <- which(rownames(thePlotData) == uniqControlName[j])
              plotData2 <- c(plotData2, median(thePlotData[rowIdx, i]))
            }
            plotData2 <- c(plotData2, NA)
          }
          uniqControlName <- sub("IVT_Kit_", "", uniqControlName)
          uniqControlName <- sub("RT_Kit_", "", uniqControlName)
          uniqControlName <- sub("Hybridization_", "Hyb", uniqControlName)
        #names(plotData2) = rep(c(uniqControlName, " "), arrayCount)
          color.label = rainbow(controlCount+1)
          if(controlShortNames[r] == "Negative") {
            ylimm <- quantile(thePlotData, probs = c(0.05, 0.95))
            boxplot(split(thePlotData, col(thePlotData)), outline = TRUE, varwidth = TRUE,
                    medcol = "red", medlwd = 5, las = 2, ylim = ylimm,
                    col = "blue", outpch = 16, outcol = "green", names = colnames(plotData))
          }
          else {
            barplot(plotData2, las = 2, col = color.label, xaxt = "n")
            #legend("topright", uniqControlName, lty = 1, lwd = 2, col = color.label)
            midPoint <- controlCount %/% 2 + 1
            atPos <- midPoint + (1.2 * (1+controlCount)) *(0:(rightCount - leftCount))
            axis(1, at = atPos, las = 2,
                 labels = colnames(plotData)[leftCount:rightCount])
            box()
          }
        }
      }
      title(main = paste(controls[ctl], plotWhat), outer = TRUE)
      savejpg(paste(jpgDir, fname, sep = ""), fWidth, fHeight)
      dev.off()
    }
  }
}

################################################
#- $Log: icpPlot.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.3  2006/03/21 23:00:17  sunya
#- Improved memory usage.
#- Control probes in the exported control signal value file are now sorted.
#- Improved handling of NA values. The statistical analysis will ignore NA.
#- If there is only one member of a subgroup, this subgroup will not appear
#- in t test.
#-
#- Revision 1.1  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
