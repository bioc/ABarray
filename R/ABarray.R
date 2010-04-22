#- $Id: ABarray.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Main preprocessing functions for ab1700 data

"ABarray" <-
function(dataFile, designFile, group, test = TRUE, impute = "avg", normMethod = "quantile", ...) {
  if(missing(designFile)) {
    print("Did you forget to provide an experiment design file?")
    return()
  }
  
  if(missing(dataFile)) {
    print("Did you provide the name of your data file?")
    return()
  }
  
  if(missing(group)) {
    print("What group do you want to perform statistics?")
    return()
  }

  kWidth <- 800  #- default figure size
  kHeight <- 600
  snThresh <- 3   #- Default S/N filtering threshold
  showHiddenPlot <- FALSE   #- There are some hidden plots normally not needed
  showControlOnly <- FALSE  #- Just plot controls no variability plots
  showControlSN <- FALSE
  esetOnly <- FALSE
  
  normUse <- "quantile"
  if (normMethod == "all") {
    normMethod <- c(normUse, "mean", "median", "trimMean", "trimAMean")
  } else {
    ifelse(normMethod == "quantile", (normUse <- "quantile"),
    ifelse(normMethod == "mean", (normUse <- "mean"),
    ifelse(normMethod == "median",  (normUse <- "median"),
    ifelse(normMethod == "trimMean", (normUse <- "trimMean"),
    ifelse(normMethod == "trimAMean", (normUse <- "trimAMean"),
    ifelse(normMethod == "none", (normUse <- "none"),
    (normMethod <- c(normUse, "mean", "median", "trimMean", "trimAMean"))))))))
  }
   
   #- Process additional arguments for sweave and hidden variables, ...
  addParam <- list(...)
  if(is.element("showHiddenPlot", names(addParam))) {
    showHiddenPlot <- addParam$showHiddenPlot
  }
  if(is.element("showControlOnly", names(addParam))) {
    showControlOnly <- addParam$showControlOnly
  }
  if(is.element("snThresh", names(addParam))) {
    snThresh <- addParam$snThresh
  }
  if(is.element("showControlSN", names(addParam))) {
    showControlSN <- addParam$showControlSN
  }
  if(is.element("esetOnly", names(addParam))) {
    esetOnly <- addParam$esetOnly
  }
  
  ##- Find out what is the field sepeartor, tab or comma
  sep <- "\t"
  if(any(grep(".csv", designFile))) {
    sep <- ","
  }
  ##-pd <- read.phenoData(designFile, sep = sep, colClasses = "factor", header = TRUE)
  pd <- read.table(designFile, sep=sep, colClasses="factor", header=TRUE)
  sampleCount <- dim(pd)[1]
  
  sep <- "\t"
  if(any(grep(".csv", dataFile))) {
    sep <- ","
  }
  ##- Preprocess data before reading the entire data set
  data <- read.table(dataFile, sep = sep, nrows = 5)
  dataColNames <- t(data[1,])
  dataColCount <- dim(data)[2]
  colClass <- rep("numeric", dataColCount)  ##- for reading data file
  colClass[1:2] <- "character"
  
  colProbe <- grep("probe", dataColNames, ignore.case = TRUE)
  colGene <- grep("gene", dataColNames, ignore.case = TRUE)
  colName <- grep("name", dataColNames, ignore.case = TRUE)
  colID <- grep("id", dataColNames, ignore.case = TRUE)

  colProbeID <- 0    #- which column contains probeID (or probeName)
  if(any(colName)) {
    if(any(intersect(colProbe, colName))) {
      colProbeID <- intersect(colProbe, colName)
    }
  }
  if(any(colID)) {
    if(any(intersect(colProbe, colID))) {
      colProbeID <- intersect(colProbe, colID)
    }
  }
  if(colProbeID == 0) {
    cat("I cannot find which column contains probeID or probeName.\n")
    cat("Processing will not start. Make sure probeID or probeName header is like this:\n")
    print(paste("ProbeID", "Probe ID", "Probe_ID", sep = " or "))
    print(paste("ProbeName", "Probe Name", "Probe_Name", sep = " or "))
    cat("Either upper case or lower case is OK\n")
    return()
  }

  colGeneID <- 0
  if(any(colName)) {
    if(any(intersect(colGene, colName))) {
      colGeneID <- intersect(colGene, colName)
    }
  }
  if(any(colID)) {
    if(any(intersect(colGene, colID))) {
      colGeneID <- intersect(colGene, colID)
    }
  }      

  ##- Find out which columns contain signals, which contains S/N and flags
  dataColNames <- gsub("assay_normalized_signal", "assay_norm_sig", dataColNames,
                       ignore.case = TRUE)
  dataColNames <- gsub("Assay Normalized Signal", "assay_norm_sig", dataColNames,
                       ignore.case = TRUE)
  dataColNames <- gsub("Assay.Normalized.Signal", "assay_norm_sig", dataColNames,
                       ignore.case = TRUE)

  colSignal <- grep("signal", dataColNames, ignore.case = TRUE)
  colSN <- grep("S/N", dataColNames, ignore.case = TRUE)
  if(length(colSN) < 1) {
    colSN <- grep("S.N", dataColNames, ignore.case = TRUE)
  }
  if(length(colSN) < 1) {
    colSN <- grep("S2N", dataColNames, ignore.case = TRUE)
  }
  colFlag <- grep("FLAG", dataColNames, ignore.case = TRUE)
  
  if(length(colSignal) > 0 & length(colSignal) != length(colSN)) {
    print(paste("Number of signal column: ", colSignal))
    print(paste("Number of S/N column: ", colSN))
    cat("Number of columns for signal is not the same as S/N ratio\n")
    cat("Processing will not be performed.\n")
    return()
  }
  if(length(colSignal) == 0) {
    cat("Which columns contain signals?\n")
    cat("No signals found. Processing stopped.\n")
    return()
  }
  col.assayNorm <- grep("assay_norm_sig", dataColNames, ignore.case = TRUE)
  col.sdev <- grep("SDEV", dataColNames, ignore.case = TRUE)
  col.cv <- grep("CV", dataColNames, ignore.case = TRUE)

  #- Figure out what is the sample,assay, arrayname
  idx.arrayName <- grep("arrayName", colnames(pd), ignore.case = TRUE)
  idx.assayName <- grep("assayName", colnames(pd), ignore.case = TRUE)
  idx.sampleName <- grep("sampleName", colnames(pd), ignore.case = TRUE)

  #- What ID to use, sampleName, or arrayName or assayName
  idUse <- c()
  idType <- "sampleName"
  if(any(idx.arrayName)) {
    idUse <- as.vector(pd[,idx.arrayName])
    idType <- "arrayName"
  }
  else if(any(idx.assayName)) {
    idUse <- as.vector(pd[,idx.assayName])
    idType <- "assayName"
  }
  else {
    idUse <- as.vector(pd[,idx.sampleName])
  }

  cat("Using", idType, "to match experiment with signal in file:", dataFile, "\n")

  ##- group samples for boxplot
  color <- rep(c("blue", "cyan", "green", "lightblue", "purple", "pink", "brown", "magenta", "orange"), 20)
  color.label <- c()
   
  ordSampleMember <- c()
  ##-pd.matrix <- pData(pd)
  grpMember <- levels(pd[, group])
  for(i in 1:length(grpMember)) {
    idx.member <- which(pd[, group] == grpMember[i])
    ordSampleMember <- c(ordSampleMember, idx.member)
    color.label <- c(color.label, rep(color[i], length(idx.member)))
  }
  pd <- pd[ordSampleMember,]
  sampleName <- as.vector(pd[,1])
  rownames(pd) <- sampleName
  
  cat("The sample names are:\n")
  print(sampleName)
   
  idUse <- idUse[ordSampleMember]  #- idUse is the name of the pData columns to use
												#- May not be same as in data file

  assayName.inData <- gsub("_SIGNAL", "", dataColNames[colSignal], ignore.case = TRUE)
  assayName.inData <- gsub("Signal ", "", assayName.inData, ignore.case = TRUE)
  assayName.inData <- gsub("Signal.", "", assayName.inData, ignore.case = TRUE)
  
  cat("AssayNames in dataFile:", dataFile, "\n")
  print(assayName.inData)
  if(length(idUse) > length(assayName.inData)) {
    cat("You have more arrays in designFile than dataFile.\n")
    cat("Processing cannot go further.\n")
    return()
  }

 	#- Sometimes, there could be more assays in data file than we actually want to analyze.
	#- For example, we decided to ignore one or two of the arrays, because the images
	#- are bad and there are reasons to exlude these arrays. But instead to delete the
	#- related columns in data file, we simply change the experimental design file to just
	#- include those arrays we need to analyze.

   #- Figure out the order of signal, sn, and flag according to sampleName order
   #- grouped by group in the above process
  ordSig <- c()
  ordSN <- c()
  ordFlag <- c()
  for(i in seq(along = idUse)) {
    theOrder <-  grep(paste("\\b", idUse[i], sep = ""), assayName.inData)
    if(length(theOrder) > 1) {
      cat("\nOne of the name is part of another name in your assay:\n")
      print(assayName.inData[theOrder])
      cat("\nProcessing stopped.\n")
      return()
    }
    ordSig <- c(ordSig, theOrder)
    theOrder <- grep(idUse[i], dataColNames[colSN])
    if(length(theOrder) > 1) {
      cat("\nOne of the S/N name is part of another S/N name:\n")
      print(dataColNames[colSN][theOrder])
      cat("\nPlease fix it. Process stopped.\n")
      return()
    }
    ordSN <- c(ordSN, theOrder)
    theOrder <- grep(idUse[i], dataColNames[colFlag])
    if(length(theOrder) > 1) {
      cat("\nOne of the FLAG name is part of another FLAG name:\n")
      print(dataColNames[colFlag][theOrder])
      cat("\nPlease fix it. Process stopped.\n")
      return()
    }
    ordFlag <- c(ordFlag, theOrder)
  }
  
  if(length(assayName.inData[-ordSig]) > 0) {
    cat("\nThe following arrays were excluded from analysis:\n")
    print(assayName.inData[-ordSig])
  }
   #- Check to see if the number of arrays match to the design file
  if(length(ordSig) < length(sampleName)) {
    cat("The number of arrays in the design is not the same as in data file\n")
    cat("    ", length(sampleName), "vs", length(ordSig), "\n")
    cat("Or the names of assays do not match.\n")
    cat("Processing is now stopped.\n")
    return(assayName.inData)
  }

   #- Let's read the data.
  cat("Reading data from", dataFile, "......\n    This may take several minutes....\n")
  flush.console()
  ##-data <- read.table(dataFile, sep = sep, header = TRUE, as.is = TRUE)
  ##- Check to see if data contains comma as thousand seperator (not decimal point)
  colRead <- rep("NULL", dataColCount)
  colRead[colSignal[1]] <- "character"
  data <- read.table(dataFile, header=TRUE, sep=sep, nrow=5000, colClasses=colRead)

  if(any(grep(",", data[, 1]))) {
    colRead[c(colProbeID, colSignal, colSn, colFlag)] <- "character"
    if(colGeneID > 0) {
      colRead[colGeneID] <- "character"
    }
    data <- read.table(dataFile, header=TRUE, sep=sep, nrow=36000, colClasses=colRead, as.is=TRUE,
        check.names=FALSE, comment.char="", na.string=c("NA", "MultipleValues", "Multiple Values"))
    cat("Checking data format:")
    flush.console()
    data = apply(data, 2, function(x) {
      cat(":"); flush.console()
      as.numeric(gsub(",", "", x))})
    cat("||\n")
  }
  else {
    colRead[c(colSignal, colSN, colFlag)] <- "numeric"
    colRead[colProbeID] <- "character"
    if(colGeneID > 0) {
      colRead[colGeneID] <- "character"
    }
    data <- read.table(dataFile, header=TRUE, sep=sep, nrow=36000, colClasses=colRead,
       check.names=FALSE, comment.char = "", na.string = c("NA", "MultipleValues", "Multiple Values"))
  }
  
  ##- Since some columns are skipped, we need to find which column contains signal, etc
  ##-
  dataColNames <- colnames(data)
  colProbe <- grep("probe", dataColNames, ignore.case=TRUE)
  colGene <- grep("gene", dataColNames, ignore.case=TRUE)
  colSig <- grep("signal", dataColNames, ignore.case=TRUE)
  colSnr <- grep("S.N", dataColNames, ignore.case=TRUE)
  colFlg <- grep("flag", dataColNames, ignore.case=TRUE)
  
  colSigOrd <- c()
  colSnrOrd <- c()
  colFlgOrd <- c()
  for(assayID in idUse) {
    colSigOrd <- c(colSigOrd, intersect(colSig, grep(assayID, dataColNames, ignore.case=TRUE)))
    colSnrOrd <- c(colSnrOrd, intersect(colSnr, grep(assayID, dataColNames, ignore.case=TRUE)))
    colFlgOrd <- c(colFlgOrd, intersect(colFlg, grep(assayID, dataColNames, ignore.case=TRUE)))
   }    
    
  cat("Finished data reading.\n")
  ##- Create a directory to put all the figures in
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

  cat("The results will be in the folder:", resultDir, "\n")
  flush.console()
  ##- Process controls
  controlType <- c("Buffer_Blank", "CL_ControlLadder", "CLFL_GridLandmark", "FL_ControlLadder", "FL_Fiducial_Cp",
                   "Hybridization_Control", "ICP_FLOnly_Control", "IVT_Kit_Control", "Manufacturing_Test_Control",
                   "Negative_Control", "Positive_Control", "RT_Kit_Control", "Blank")
  
  controlName <- c("Buffer_Blank", "CL_Ladder", "CLFL_Grid", "FL_Ladder", "FL_Fiducial",
                   "Hybridization", "ICP_FLOnly", "IVT_Kit", "Manufacturing",
                   "Negative", "Positive", "RT_Kit", "Blank")
  controlPlot = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, showHiddenPlot, TRUE, TRUE, TRUE, FALSE)
  names(controlType) <- controlName
  rowControl = c()
  if(colGeneID != 0) {
    geneID <- as.vector(data[, colGeneID])
    rowControl <- grep("NULL-GENE", data[,colGeneID])
    rowControl <- c(rowControl, which(data[, colGeneID] == ""))
  }
  else {
    geneID <- as.vector(data[, colProbeID])
    for(i in seq(along = controlType)) {
      rowControl <- c(rowControl, grep(controlType, data[, colProbeID], ignore.case = TRUE))
    }
  }
  
  controlData <- data[rowControl,]
  ##- Remove control data for processing
  if(any(rowControl)) {
    data <- data[-rowControl,]
    geneID <- geneID[-rowControl]
  }
  
  ##rawSig <- as.matrix(data[, colSignal[ordSig]])
  ##sn <- as.matrix(data[, colSN[ordSN]])
  ##flag <- as.matrix(data[, colFlag[ordFlag]])
  rawSig <- as.matrix(data[, colSigOrd])
  sn <- as.matrix(data[, colSnrOrd])
  flag <- as.matrix(data[, colFlgOrd])
  
  colnames(rawSig) <- sampleName
  rownames(rawSig) <- data[, colProbe]
  colnames(sn) <- paste("SN", sampleName)
  rownames(sn) <- data[, colProbe]
  colnames(flag) <- paste("Flag", sampleName)
  rownames(flag) <- data[, colProbe]
  data <- NULL
  
  ##- Processing control data and plot
  conData <- controlData[, c(colProbe, colSigOrd)]
  colnames(conData) <- c("ProbeID", sampleName)

  controlToExport <- c("Hybridization_Control", "IVT_Kit_Control", "RT_Kit_Control")
  idxControlToExport <- c()
  for(i in seq(along = controlToExport)) {
    idxControlToExport <- c(idxControlToExport, grep(controlToExport[i], conData[, 1]))
  }
  probeSortToExport <- sort(conData[idxControlToExport, colProbe], index = TRUE)$ix
  write.table(conData[idxControlToExport[probeSortToExport],],
    file=paste(dataDir, "ControlRawData.csv", sep = ""), row.names=FALSE, col.names=TRUE, sep=",")
   
  icpPlot(conData, pdfDir = pdfDir, jpgDir = jpgDir)
  if(showControlSN) {
    conData <- controlData[, c(colProbe, colSnrOrd)]
    colnames(conData) <- c("ProbeID", sampleName)
    icpPlot(conData, plotWhat = 'SN', pdfDir = pdfDir, jpgDir = jpgDir)
  }
  conData <- NULL
  controlData <- NULL
  gc()
   
  write.csv(rawSig, file = paste(dataDir, "ExpressionRawData.csv", sep = ""))
   
  cat("\nPerform basic analysis for non-normalized data ...\n")
                                        #- probes detectable
  probeDetect <- apply(sn, 2, function(x) {sum(x >= snThresh, na.rm = TRUE)})
  names(probeDetect) <- sampleName
  ymin <- min(probeDetect, na.rm = TRUE)
  ymax <- max(probeDetect, na.rm = TRUE)
  ydiff <- ymax - ymin
  fname <- paste("QC_DetectableProbeSN", snThresh, sep = "")
  
  chLength <- max(nchar(names(probeDetect)) * 1.2, 5)
  print(paste("Creating barplot for probes detectable ...", fname))
  probeCount = dim(sn)[1]
  arrayCount = dim(sn)[2]
  pctLabel = format(100 * probeDetect / probeCount, digits = 2)
  pdf(file = paste(pdfDir, fname, ".pdf", sep = ""), width = 10, height = 8)
  dev.control("enable")
  par(mar = c(chLength, 5, 4, 2))
  barplot(probeDetect, ylim = c(ydiff, ymax + ydiff), xpd = FALSE, las = 2,
          col = color.label, main = paste("Probe Detectable (S/N >=", snThresh, ")"))
  text(seq(arrayCount) * 1.2 - 0.5, probeDetect + ydiff / 5, labels = pctLabel)
  savejpg(paste(jpgDir, fname, sep = ""), width = 800, height = 600)
  dev.off()
  probeDetect <- NULL
  
  ##- Create new ExpressionSet
  pheno <- as(pd, "AnnotatedDataFrame")
  esetUse <- new("ExpressionSet", phenoData=pheno, exprs=rawSig, snDetect=sn, flags=flag)

  if(! showControlOnly) {
    if(!esetOnly) {
      doPlotEset(esetUse, group, name = "Raw", test = FALSE, rawData = TRUE, snThresh = snThresh)
    }
    esetUse <- NULL
    for(i in seq(along = normMethod)) {
      if(normMethod[i] == "trimAMean") {
        normData <- log2(qnNormalize(rawSig, snr = sn, method = normMethod[i], snThresh = snThresh))
      }
      else if (normMethod[i] == "none") {
        normData <- log2(rawSig)
      }
      else {
        normData <- log2(qnNormalize(rawSig, method = normMethod[i]))
      }
    
      if(impute == "avg" || impute == "knn") {
        flagSet <- which(flag > 5000)
        rawSig[flagSet] <- NA

        imputedValue <- imputeFlag(normData, pd = pd, group = group, impute = impute)
        normData <- imputedValue$imputedData
        flagCount <- imputedValue$flagCount
        if(any(flagCount > 0)) {
          flagged <- which(flagCount > 0)
          if(impute == "knn") {
            print(paste("Array ", sampleName[flagged], " has ", flagCount, " FLAGS > 5000", sep = ""))
            cat("The signals of each flagged probe were treated as missing values.\n")
            cat("Missing values were estimated using knn imputation algorithm as described in\n")
            cat("Troyandskaya et al (2001) Bioinformatics 17(6):520-525\n")
            flush.console()
          }
          else {
            rowImputed <- imputedValue$rowImputed
            if(length(rowImputed) > 0) {
              print(paste("Array ", sampleName[flagged], " has ", flagCount, " FLAGS > 5000", sep = ""))             
              cat("The signals of each flagged probe were replaced with average signals\n",
                  "   from replicate arrays within the same subgroup (",
                  paste(grpMember, collapse = ", "), ").", sep = "")
              filename <- paste(dataDir, "ProbeImputed_", group, ".csv", sep = "")
              write.table(rownames(rawSig)[rowImputed], file = filename, sep = ",", row.names = FALSE, col.names = FALSE)
              cat("\nThe imputed probes were written to file:", filename, "\n")
            }
            else {
              cat("Imputation was not performed. There is no enough data to perform imputation.\n")
            }
            flush.console()
          }
        }
      }
      colnames(normData) <- paste(normMethod[i], colnames(normData))
    
                                        #- write quantile normalized data to file
      csvName <- paste("Expression_", normMethod[i], "Normalized.csv", sep = "")
      normExp <-  cbind(ProbeID = rownames(normData), GeneID = geneID, normData, sn, flag)
      write.table(normExp, file = paste(dataDir, csvName, sep = ""), sep = ",",
                  col.names = TRUE, row.names = FALSE)
      cat("\n", normMethod[i], "normalized data was saved to:\n    ", paste(dataDir, csvName, sep = ""), "\n")
      flush.console()
      normExp <- NULL
      rm("normExp")
      gc()

      colnames(normData) = colnames(rawSig)
      eset <- new("ExpressionSet", phenoData=pheno, exprs=normData, snDetect=sn, flags=flag)
      if(normMethod[i] == normUse) {
        esetUse <- eset
      }
      esetFile <- sub(".txt", "", designFile);  esetFile = sub(".csv", "", esetFile)
      esetFile <- paste(esetFile, "_eset_", normMethod[i], format(Sys.time(), "%d%b%Y"), ".Rdata", sep = "")
      save(eset, file = esetFile)
      cat("\nThe expression data object was saved to R workspace file:\n    ", esetFile, "\n")
    }
    
    cat("\nPerform data analysis on", normUse, "normalized data ...\n")
    flush.console()
    data <- NULL
    rawSig <- NULL
    normData <- NULL
    sn <- NULL
    flag <- NULL
    geneID <- NULL
    gc()
    if(!esetOnly) {
      doPlotEset(esetUse, group, name = normUse, test = test, ...)
    }
  }  
  return(esetUse)	
}

#############################################################################
#- $Log: ABarray.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.1  2006/06/06 17:14:36  sunya
#- Changed ab1700 to ABarray in all releted files to be more compatible
#- in name convention to Bioconductor.
#-
#- Revision 1.16  2006/05/03 16:29:59  sunya
#- In phenoData, rownames of the pData is now replaced with sampleNames, rather
#- than the order it was read previously to avoid the conflict with Biobase in
#- Bioconductor release 1.8.
#-
#- Revision 1.15  2006/04/17 21:06:34  sunya
#- Correcting typos in the calling qnNormalize function in ab1700 and some
#- typos in qnNormalize.Rd document.
#-
#- Revision 1.14  2006/04/17 19:04:56  sunya
#- Added vignette files. Added normalization methods for mean, median, trimmed mean.
#-
#- Revision 1.13  2006/03/27 23:09:04  sunya
#- Added 'no imputation' option in GUI interface. Additional cosmetic changes.
#-
#- Revision 1.11  2006/03/20 23:00:16  sunya
#- Missing values (NA) are kept even after imputation. But most downstream
#- analysis will remove these NA values.
#- The hierarchical clustering probe list and expression values are now
#- exported into csv file.
#- Changed file naming convention. QC plot will begin QC_, and statistical
#- analysis plot will begin ST_ in the figures.
#- Values for control probes are saved into csv file.
#-
#- Revision 1.10  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
#- Revision 1.9  2005/10/25 23:48:56  sunya
#- Changed imputation occurance. Imputation now occurs after normalization.
#- Added function imputeFlag, previously in ab1700.R itself.
#-
#- Revision 1.8  2005/10/17 21:05:00  sunya
#- Added ab1700gui document. Changed DESCRIPTION version to 1.0.1
#-
#- Revision 1.7  2005/10/14 22:38:33  sunya
#- Added support for additional filtering options. snThresh and detectSample
#- can be specified in the command line and also in GUI interface.
#-
#- Revision 1.6  2005/10/13 22:06:54  sunya
#- Updated documents. Matched pdf and jpg figure sizes.
#-
#- Revision 1.5  2005/10/07 18:05:43  sunya
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
#- Revision 1.4  2005/08/22 22:43:55  sunya
#- Fixed avg imputation error. The imputed values were put into wrong columns.
#- Updated pdf documents.
#-
#- Revision 1.3  2005/08/02 23:40:30  sunya
#- Output median signals for control probes.
#- Impute missing value function takes parameter, if knn, use impute.knn from Bioconductor library impute
#- function. If ave, use average value computed within subarrays, if there are more than 2 arrays without
#- missing values. ave if the default.
#-
#- Revision 1.2  2005/07/29 17:34:37  sunya
#- Added support for probe detectable graph; removing probes with no expression values
#-
