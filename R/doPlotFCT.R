#- $Id: doPlotFCT.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Main function to perform fold change and t test.

"doPlotFCT" <-
function(eset, group, grpMember, order1=NULL,order2=NULL, detectSample = 0.5, snThresh = 3, ...) {
  require(Biobase)
  require(multtest)

   #- snThresh S/N threshold value for filtering
   #- detectSample Percentage of samples that meet snThresh for t test analysis.
  paired <- FALSE
  filePrefix <- paste("DetectPct", detectSample * 100, sep = "")
  if(!is.null(order1) && !is.null(order2)) {
    filePrefix <- paste("paired_", filePrefix, sep = "")
    paired <- TRUE
    if(length(order1) != length(order2)) {
      cat("Number of samples for order1 and order2 is not equal,\n",
          "paired t test cannot be performed.\n")
      return()
    }
  }

  probeLimit <- 800
  binColor <- c("cyan", "green", "blue", "black", "red")
  fcBin <- list("FC 1.0-1.2" = c(1.0, 1.2), "FC 1.2-1.6" = c(1.2, 1.6),
                "FC 1.6-2.0" = c(1.6, 2.0),
                "FC 2.0-4.0" = c(2.0, 4.0), "FC > 4" = c(4.0, Inf))
  pValPreset <- c(0.01, 0.05)
  fdrPreset <- c(0.01, 0.05, 0.1, 0.25)
  
  #- Process additional arguments for sweave and hidden variables, ...
  addParam <- list(...)
  if(is.element("probeLimit", names(addParam))) {
    probeLimit <- addParam$probeLimit
  }
  if(is.element("binColor", names(addParam))) {
    if(length(addParam$binColor) == length(fcBin)) {
      binColor <- addParam$binColor
    }
  }
  if(is.element("pValPreset", names(addParam))) {
    pValPreset <- addParam$pValPreset
  }
  if(is.element("fdrPreset", names(addParam))) {
    fdrPreset <- addParam$fdrPreset
  }

  
  kWidth <- 800  #- default figure size
  kHeight <- 600
  
  fdpPreset <- c(pValPreset, fdrPreset)
  fdpName <- c(rep("pVal", length(pValPreset)), rep("FDR", length(fdrPreset)))
  fdpResult <- rep(0, length(fdpPreset))
   
  pd <- pData(eset)
  if(missing(grpMember)) {
    grpMember <- levels(factor(pd[, colnames(pd) == group]))
  }
  subgrpCount <- length(grpMember)

  if(subgrpCount <= 1) {
    print("Only one group member available, no t test will be performed")
    return()
  }

  eset <- getMemberEset(eset, group, grpMember)

  #- Check to see if there is S/N ratio values
  snPresent <- TRUE
  sn <- assayDataElement(eset, "snDetect")
  if(dim(sn)[2] <= 1 || sum(is.na(sn[,])) > 1) {
    snPresent <- FALSE
  }

	#- Create a directory to put all the figures in
  resultDir <- paste("Result_", gsub(" ", "", group), "/", sep = "")
  if(! file.exists(resultDir)) {
    dir.create(resultDir, showWarnings = FALSE)
  }
  dataDir <- paste(resultDir, "DataResult/", sep = "")
  pdfDir <- paste(resultDir, "pdfFigure/", sep = "")
  jpgDir <- paste(resultDir, "jpgFigure/", sep = "")
  if(!file.exists(dataDir)) {
    dir.create(dataDir, showWarning = FALSE)
  }
  if(!file.exists(pdfDir)) {
    dir.create(pdfDir, showWarning = FALSE)
  }
  if(!file.exists(jpgDir)) {
    dir.create(jpgDir, showWarning = FALSE)
  }
  
  #-print(paste("Count:", subgrpCount))
  for(grp1 in 1:(subgrpCount - 1)) {
    if(sum(pd[, group] == grpMember[grp1]) > 1) {
      for(grp2 in (grp1 + 1):subgrpCount) {
        if(sum(pd[, group] == grpMember[grp2]) > 1) {
          member <- c(grpMember[grp1], grpMember[grp2])
          memberName <- paste(member[1], member[2], sep = "-")
          cat("\nPerforming t test between", member[1], "and", member[2], "...\n")
          flush.console()
          
          idxProbe <- 1:dim(exprs(eset))[1]
          if(snPresent) {
            snT <- snSummary(eset, snThresh, group, member)
            idx1 <- which(snT[,1] >= detectSample)
            idx2 <- which(snT[,2] >= detectSample)
            idxProbe <- union(idx1, idx2)
          }      
          probeCount <- length(idxProbe)
          esetUse <- eset[idxProbe,]
          colMember1 <- which(pData(esetUse)[, group] == member[1])
          colMember2 <- which(pData(esetUse)[, group] == member[2])
          sd1 <- sd(t(exprs(esetUse[, colMember1])), na.rm = TRUE)
          sd2 <- sd(t(exprs(esetUse[, colMember2])), na.rm = TRUE)
          rowSd1 <- union(which(is.na(sd1)), which(sd1 == 0))
          rowSd2 <- union(which(is.na(sd2)), which(sd2 == 0))
          rowSd <- union(rowSd1, rowSd2)
          if(length(rowSd) > 0) {
            cat(length(rowSd), "probes had sdev = 0, they were removed in the t test\n")
            print(featureNames(esetUse)[rowSd])
            esetUse <- esetUse[-rowSd,]
            idxProbe <- idxProbe[-rowSd]
            probeCount <- length(idxProbe)        
          }
           #- Calculate fold change and t test using simpleaffy
           #-fct <- get.fold.change.and.t.test(esetUse, group, member, a.order=order1, b.order=order2)
          tRes <- apply(exprs(esetUse), 1, function(x) {
            if(paired) {
              t.test(x[colMember1[order1]], x[colMember2[order2]], na.rm = TRUE, paired = TRUE)$p.value
            }
            else {
              t.test(x[colMember1], x[colMember2], na.rm = TRUE)$p.value
            }
          })
          fcRes <- apply(exprs(esetUse), 1, function(x) {
            if(paired) {
              mean(x[colMember1[order1]] - x[colMember2[order2]], na.rm = TRUE)
            }
            else {
              mean(x[colMember1], na.rm = TRUE) - mean(x[colMember2], na.rm = TRUE)
            }
          })
      
       #- Volcano plot
          print("Creating volcano plot ...")
          for(i in seq(along = pValPreset)) {
            fileName <- paste("ST_", filePrefix, "VolcanoPlot_p", sub("0.", "", pValPreset[i]), memberName, sep = "")
            pdf(file = paste(pdfDir, fileName, ".pdf", sep = ""), width = kWidth/100, height = kHeight/100)
            dev.control("enable")
            idxPval <- which(tRes <= pValPreset[i])
            idxFC <- which(abs(fcRes) >= 1)
            idxUse <- intersect(idxPval, idxFC)
            plot(fcRes[-idxUse], -log10(tRes[-idxUse]), yaxt = "n", pch = 16, col = "blue",
                 xlab = paste("log2 of FC", memberName), ylab = "p values",
                 main = paste("Volcano Plot ", memberName))
      
            points(fcRes[idxUse], -log10(tRes[idxUse]), col = "red", pch = 16)
            axis(2, 1:5, labels = c(0.1, 0.01, 0.001, 0.0001, 0.00001))
            abline(h = -log10(pValPreset[i]), lty = 2)
            abline(v = c(-1, 1), col = "blue", lty = 2)
            savejpg(paste(jpgDir, fileName, sep = ""), width = kWidth, height = kHeight)
            dev.off()
          }
       
       #- Apply FDR
          result.adj <- mt.rawp2adjp(tRes, "BH")
          adjpFDR <- result.adj$adjp[order(result.adj$index),]

                                        #- write fc, tt, FDR to csv file
          fcBinInfo <- vector(mode = "character", length = probeCount)
          for(bin in seq(along = fcBin)) {
            idxBin <- intersect(which(abs(fcRes) > log2(fcBin[[bin]][1])),
                                which(abs(fcRes) <= log2(fcBin[[bin]][2])))
            fcBinInfo[idxBin] <- names(fcBin)[bin]
          }

          output <- cbind("ProbeID" = featureNames(esetUse),  "FC" = 2 ^ fcRes, "Log2(FC)" = fcRes,
                          "p value" = adjpFDR[,1], "FDR (BH)" = adjpFDR[,2], "FC Bin" = fcBinInfo)
          colnames(output) = c("ProbeID", paste("FC (", member[1], "/", member[2], ")", sep = ""),
                    "Log2(FC)", "p value", "FDR (BH)", "FC Bin")
          outName <- paste(dataDir, "ST_", filePrefix, "FCpvalFDR_", memberName, ".csv", sep = "")
          print(paste("Writing t test result to file: ", outName))
          flush.console()
          write.table(output, file = outName, sep = ",", row.names = FALSE, col.names = TRUE)
          output <- NULL
       
       #- Preset an interval to check how many probes pass FDR at specified level
          fdp <- data.frame(preset = fdpPreset, name = fdpName, result = fdpResult)
          cat("\nThe t test was performed if the probe shows S/N >= ", snThresh, " in ",
              detectSample *100, "% or more samples\n",
              "    in either ", member[1], " or ", member[2], "\n", sep = "")
          cat("The number of probes after filtering is:", probeCount, "\n\n")
          flush.console()

          for(i in 1:dim(fdp)[1]) {
            idxProbe <- c()
            fdr <- paste(fdp$name[i], fdp$preset[i]) # what and which level
            if(fdp$name[i] == "pVal") {
              pvCount <- sum(adjpFDR[,1] <= fdp$preset[i], na.rm = TRUE)
              fdp$result[i] <- pvCount
              idxProbe <- which(adjpFDR[,1] <= fdp$preset[i])
              cat("p=", fdp$preset[i], ", significant probes: ", pvCount, "\n", sep = "")
            }
            else {
              fdrCount <- sum(adjpFDR[,2] <= fdp$preset[i], na.rm = TRUE)
              fdp$result[i] <- fdrCount
              idxProbe <- which(adjpFDR[,2] <= fdp$preset[i])
              cat("FDR=", fdp$preset[i], ", significant probes: ", fdrCount, "\n", sep = "")
            }
         
       #- Make plots for no FDR (p = 0.01, and 0.05) first,
       #- then, make plots after FDR adjustment
            if(length(idxProbe) > 1) {
                 #fctUse <- fct[idxProbe,]
              tResUse <- tRes[idxProbe]
              fcResUse <- fcRes[idxProbe]
              rowMean <- rowMeans(exprs(esetUse[idxProbe,]), na.rm = TRUE)
                                        #- Divide fc into 5 bins:
                                        #- fc 1 - 1.2;  1.2 - 1.6;  1.6 - 2.0;  2.0 - 4.0; > 4.0
              idxBinProbe <- list()
              for(bin in seq(along = fcBin)) {
                theIdx <- intersect(which(abs(fcResUse) > log2(fcBin[[bin]][1])),
                  which(abs(fcResUse) <= log2(fcBin[[bin]][2])))
                idxBinProbe[[bin]] <- theIdx
              }
                                        #- Create Fold Change plot
              print(paste("Creating Fold Change plot for:", memberName, fdr, " ..."))
              flush.console()
              fileName <- paste("ST_", filePrefix, "fct_", memberName, "_", sub(" 0.", "", fdr), sep = "")
              pdf(file = paste(pdfDir, fileName, ".pdf", sep = ""), width = kWidth/100, height = kHeight/100)
              dev.control("enable")
              plot(min(rowMean), 0, type = "n", xlim = c(min(rowMean), max(rowMean)), ylim = c(-4, 4),
                   xlab = "Average of Log2 Signal", ylab = "Fold Change (Log2)",
                   main = paste("Average FC ", memberName, " (", length(idxProbe), " probes) ",
                     fdr, sep = ""))
              legendText <- rep("", length(fcBin))
              for(bin in seq(along = fcBin)) {
                points(rowMean[idxBinProbe[[bin]]], fcResUse[idxBinProbe[[bin]]],
                       col = binColor[bin], pch = 21, bg = binColor[bin])
                legendText[bin] = paste(names(fcBin)[bin], ", ", length(idxBinProbe[[bin]]), "    ", sep = "")
              }
              legend("topright", legendText, col = binColor, pt.bg = binColor, pch = 21, pt.cex = 2)
          
              savejpg(paste(jpgDir, fileName, sep = ""), width = kWidth, height = kHeight)
              dev.off()

                                        #- Create cluster
              if(length(idxProbe) < probeLimit & length(idxProbe) > 2) { 
                exps <- exprs(esetUse[idxProbe,])
                dist <- c("Euclidean")
                for(step in 1:length(dist)) {
                  print(paste("Creating cluster heatmap for ", fdr, " using ", dist[step], " ...", sep = ""))
                  flush.console()
                  filename <- paste("ST_", filePrefix, "hcluster_", memberName, "_", sub(" 0.", "", fdr),
                                    dist[step], sep = "")
                  pdf(file = paste(pdfDir, filename, ".pdf", sep = ""), width = kWidth/100, height = kWidth/100)
                  dev.control("enable")
                  hv <- hclusterPlot(exps, paste("Cluster ", fdr, group), dist = dist[step])
                  savejpg(paste(jpgDir, filename, sep = ""), width = kWidth, height = kWidth)
                  dev.off()
              
              #- Export the list of probes
                  fname = paste(dataDir, filename, "ProbeList.csv", sep = "")
                  print(paste("Exporting lists to", fname))
                  write.csv(exps[rev(hv$rowInd), hv$colInd], file = fname)
                }
              }
              else {
                print(paste("Clustering was not performed, because there are", length(idxProbe), "probes."))
                print(paste("    The limit of probe counts for clustering is:", probeLimit))
              }
            }
          }
        }
      }
    }
  }
}

##############################################################
#- $Log: doPlotFCT.R,v $
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
