#- $Id: ABarrayGUI.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Open up a GUI interface to take the parameters and pass to ABarray function.
ABarrayGUI <- function() {
  options(show.error.messages = TRUE)
  
  tt <- tktoplevel()
  tkwm.title(tt, "AB Genome Survey Array Data Analysis")

  designFile <- ""
  designFileVar <- tclVar("")
  dataFile <- ""
  dataFileVar <- tclVar("")

  group <- ""
  snThresh <- ""
  detectSample <- ""
  groupNames <- c("Which Group")
  groupSpecVar <- tclVar(1)
  ttest <- FALSE
  ttestVar <- tclVar(1)
  impute <- "None"
  imputeVar <- tclVar(1)
  snVar <- tclVar("3")
  detectSampleVar <- tclVar("0.5")

  titleFont <- "Helvetica 14"
  normalFont <- "Helvetica 10"

  #########################################################################
  #######################################################################

  postMsg <- function(msg) {
    tkconfigure(messageText, state = "normal")
    tkinsert(messageText, "end", msg)
    tkconfigure(messageText, state = "disabled")
    flush.console()
  }

  findDesignFile <- function() {
    tclvalue(designFileVar) <- tclvalue(tkgetOpenFile())
    designFile <- tclvalue(designFileVar)
    sep <- "\t"
    if(any(grep(".csv", designFile))) {
		sep <- ","
    }
    design <- read.table(designFile, header = TRUE, sep = sep)
    designColname <- colnames(design)
    colExclude <- c()
    nameExclude <- c("sample", "assay", "array", "file")
    for(i in seq(along = nameExclude)) {
      idx <- grep(nameExclude[i], designColname, ignore.case = TRUE)
      if(any(idx)) {
        colExclude = c(colExclude, idx)
      }
    }    
    if(length(designColname[-colExclude]) < 1) {
      postMsg("Your designFile does not contain grouping column.\n")
    }
    else {
      groupNames <- designColname[-colExclude]
      postMsg("Please make sure t test group is the one you intended.\n")
      tkdelete(groupBox, 0)
      for(i in 1:length(groupNames)) {
        tkinsert(groupBox, "end", groupNames[i])
      }
      tkselection.set(groupBox, 0) 
    }
    setwd(dirname(designFile))
  }
  
  findDataFile <- function() {
    tclvalue(dataFileVar) = tclvalue(tkgetOpenFile())
    dataFile <- tclvalue(dataFileVar)
  }

  getImputeFunc <- function() {
    if(tclvalue(imputeVar) == 1) {
      tkconfigure(imputeAvgRbtn, state = "normal")
    }
    else if(tclvalue(imputeVar) == 2) {
      tkconfigure(knnRbtn, state = "normal")
    }
  }

  getTestFunc <- function() {
    ttest <- tclvalue(ttestVar)
  }
  
  getGroupFunc <- function() {
    ##-group = tclvalue(groupVar)
    ##-if(tclvalue(groupVar) == 1) {
    ##-  tkconfigure(groupEntry, state = "normal")
    ##-}
    group <- get(groupNames)[as.numeric(tkcurselection(groupBox))+1]
  }

  getSNfunc <- function() {
    sn <- tclvalue(snVar)
  }

  getDetectSampleFunc <- function() {
    detectSample <- tclvalue(detectSampleVar)
  }
  
  executeFunc <- function() {
    postMsg("The analysis is in progress ...\nPlease wait ...\n")

    designFile <- tclvalue(designFileVar)
    dataFile <- tclvalue(dataFileVar)
    ##-if(tclvalue(groupSpecVar) == 1) {
    ##-  group = tclvalue(groupVar)
    ##-}
    ##-size <- as.numeric(tksize(groupBox))
    ##-groupNames <- c()
    ##-for(i in 1:size) {
    ##- groupNames[i] <- tclvalue(tkget(groupBox, (i - 1)))
    ##-}
    ##group <- groupNames[as.numeric(tkcurselection(groupBox))+1]
    group <- tclvalue(tkget(groupBox, as.numeric(tkcurselection(groupBox))))
    
    if(tclvalue(ttestVar) == 1) {
      ttest <- TRUE
    }
    if(tclvalue(imputeVar) == 1) {
      impute <- "avg"
    }
    else if (tclvalue(imputeVar) == 2) {
      impute <- "knn"
    }
    else {
      impute <- "None"
    }

    snThresh <- as.numeric(tclvalue(snVar))
    detectSample <- as.numeric(tclvalue(detectSampleVar))

    ##-postMsg(paste("GroupName:", groupNames, "\n"))
    postMsg(paste("Group:", group, "\n"))
    postMsg(paste("DesignFile:", designFile, "\n"))
    postMsg(paste("DataFile:", dataFile, "\n"))
    postMsg(paste("ttest:", ttest, "\n"))
    postMsg(paste("Impute:", impute, "\n"))
    postMsg(paste("snThresh:", snThresh, "\n"))
    postMsg(paste("detectSample:", detectSample, "\n"))

    setwd(dirname(designFile))
    eset <- ABarray(dataFile, designFile, group, test = ttest, impute = impute,
                   snThresh = snThresh, detectSample = detectSample)

    cat("\n\n******Analysis completed******\n")
    postMsg("The analysis is completed.\n")
    tkdestroy(tt)
  }

  topFrame <- tkframe(tt, borderwidth = 2)

  usage <- paste("\nTo perform data analysis on output from Applied Biosystems Genome Survey",
    "Array, please complete the following\n")
  usageFrame <- tkframe(topFrame, relief = "raised", bd = 2)
  tkpack(tklabel(usageFrame, text = usage), anchor = "w")
  #tkgrid(tklabel(tt, text = usage))
  tkpack(usageFrame, fill = "x", expand = TRUE)
  
  #- frame for design file
  designFrame = tkframe(topFrame, relief = "raised", bd = 2)
  tkpack(tklabel(designFrame, text = "Please choose design file:", font = titleFont),
         anchor = "w")
  designFileFrame = tkframe(designFrame, relief = "groove", bd = 2)
  designFileAskFrame = tkframe(designFileFrame)
  designFileAskLabel = tklabel(designFileAskFrame, text = "DesignFile:", font = normalFont)
  designFileAskEntry = tkentry(designFileAskFrame, textvariable = designFileVar, font = normalFont,
    justify = "center")
  tkpack(designFileAskLabel, side = "left")
  tkpack(designFileAskEntry, side = "right", fill = "x", expand = TRUE)
  tkpack(designFileAskFrame, fill = "x", expand = TRUE)
  fButtonFrame = tkframe(designFileFrame)
  fileButton = tkbutton(fButtonFrame, text = "Browse File", font = normalFont,
    command = findDesignFile)
  tkgrid(fileButton)
  tkpack(fButtonFrame, anchor = "e")
  tkpack(designFileFrame, fill = "x")
  tkpack(designFrame, fill = "x", expand = TRUE)
  
  #- frame for data file
  dataFrame = tkframe(topFrame, relief = "raised", bd = 2)
  tkpack(tklabel(dataFrame, text = "Please choose data file:", font = titleFont),
         anchor = "w")
  dataFileFrame = tkframe(dataFrame, relief = "groove", bd = 2)
  dataFileAskFrame = tkframe(dataFileFrame)
  dataFileAskLabel = tklabel(dataFileAskFrame, text = "DataFile:", font = normalFont)
  dataFileAskEntry = tkentry(dataFileAskFrame, textvariable = dataFileVar, font = normalFont,
    justify = "center")
  tkpack(dataFileAskLabel, side = "left")
  tkpack(dataFileAskEntry, side = "right", fill = "x", expand = TRUE)
  tkpack(dataFileAskFrame, fill = "x", expand = TRUE)
  fButtonFrame = tkframe(dataFileFrame)
  fileButton = tkbutton(fButtonFrame, text = "Browse File", font = normalFont,
    command = findDataFile)
  tkgrid(fileButton)
  tkpack(fButtonFrame, anchor = "e")
  tkpack(dataFileFrame, fill = "x")
  tkpack(dataFrame, fill = "x", expand = TRUE)
  
  #- frame for analysis options
  optionFrame = tkframe(topFrame, relief = "raised", bd = 2)
  tkpack(tklabel(optionFrame, text = "Analysis Options:", font = titleFont), anchor = "w")
  optionInsertFrame = tkframe(optionFrame, relief = "groove", bd = 2)
  
  #- what impute funtion to use
  imputeLabelFrame = tkframe(optionInsertFrame)
  tkpack(tklabel(imputeLabelFrame, text = "Choose impute function:",
                 font = normalFont), anchor = "w")
  tkpack(imputeLabelFrame, fill = "x", expand = TRUE)
  imputeAvgFrame = tkframe(optionInsertFrame, padx = 10)
  imputeAvgRbtn = tkradiobutton(imputeAvgFrame, text = "Subgroup average", font = normalFont,
    value = 1, variable = imputeVar, command = getImputeFunc)
  tkpack(imputeAvgRbtn, side = "left", anchor = "w")
  tkpack(imputeAvgFrame, fill = "x", expand = TRUE)

  ##- Knn imputation
  ##knnFrame = tkframe(optionInsertFrame, padx = 10)
  ##knnRbtn = tkradiobutton(knnFrame, text = "KNN imputation", font = normalFont,
  ##  value = 2, variable = imputeVar, command = getImputeFunc)
  ##tkpack(knnRbtn, side = "left", anchor = "w")
  ##tkpack(knnFrame, fill = "x", expand = TRUE)

  #- No impute selection
  noImpFrame = tkframe(optionInsertFrame, padx = 10)
  noImpRbtn = tkradiobutton(noImpFrame, text = "No imputation", font = normalFont,
    value = 3, variable = imputeVar, command = getImputeFunc)
  tkpack(noImpRbtn, side = "left", anchor = "w")
  tkpack(noImpFrame, fill = "x", expand = TRUE)

  #- Which group to use for comparison
  groupFrame = tkframe(optionInsertFrame)
  groupBox <- tklistbox(groupFrame, height=4, selectmode="single", background="white")
  tkpack(tklabel(groupFrame, text="Which group for analysis"), anchor="w")
  for(i in 1:length(groupNames)) {
    tkinsert(groupBox, "end", groupNames[i])
  }
  tkselection.set(groupBox, 1)
  
  ##-groupCbtn = tkradiobutton(groupFrame, text = "Specify array grouping:", font = normalFont,
  ##-  variable = groupSpecVar, value = 1, command = getGroupFunc)
  ##-groupEntry = tkentry(groupFrame, textvariable = groupVar, font = normalFont,
  ##-  width = 20, state = "normal", justify = "center")
  ##-tkpack(groupCbtn, side = "left", anchor = "w")
  ##-tkpack(groupEntry, side = "left")
  tkpack(groupBox, expand=TRUE)
  tkpack(groupFrame, fill = "x", expand = TRUE)

    #- whether to use t test
  ttestFrame = tkframe(optionInsertFrame)
  ttestCbtn = tkcheckbutton(ttestFrame, text = "Use t test", font = normalFont,
    variable = ttestVar, command = getTestFunc)
  ttestLabel = "Should t test be performed within the specified group?"
  tkpack(ttestCbtn, side = "left", anchor = "w")
  tkpack(tklabel(ttestFrame, text = ttestLabel), side = "left", anchor = "w")
  tkpack(ttestFrame, fill = "x", expand = TRUE)

  tkpack(optionInsertFrame, fill = "x", expand = TRUE)
  tkpack(optionFrame, fill = "x", expand = TRUE)

  #- Filtering options
  filterFrame = tkframe(optionInsertFrame)
  tkpack(tklabel(filterFrame, text = "Specify filtering options:", font = normalFont),
                 anchor = "w")
  tkpack(filterFrame, fill = "x", expand = TRUE)
  snFrame = tkframe(filterFrame, padx = 10)
  tkpack(tklabel(snFrame, text = "S/N threshold >="), side = "left", anchor = "w")
  snEntry = tkentry(snFrame, textvariable = snVar, font = normalFont, width = 8,
    state = "normal", justify = "center")
  tkpack(snEntry, side = "left")
  tkpack(snFrame, fill = "x", expand = TRUE)
  
  detectSampleFrame = tkframe(filterFrame, padx = 10)
  tkpack(tklabel(detectSampleFrame, text = "% Detect Samples:"), side = "left", anchor = "w")
  detectSampleEntry = tkentry(detectSampleFrame, textvariable = detectSampleVar,
    font = normalFont, width = 8, state = "normal", justify = "center")
  tkpack(detectSampleEntry, side = "left")
  tkpack(detectSampleFrame, fill = "x", expand = TRUE)
  tkpack(filterFrame, fill = "x", expand = TRUE)
  
  #- Action button
  actionFrame = tkframe(topFrame, relief = "raised", bd = 2)
  tkpack(tklabel(actionFrame, text = "Perform Data Analysis:", font = titleFont), anchor = "w")
  actionInsertFrame = tkframe(actionFrame, relief = "groove", bd = 2)
  executeBtn = tkbutton(actionInsertFrame, text = "Perform Analysis", font = titleFont,
    command = executeFunc)
  tkpack(executeBtn, side = "right", anchor = "e")
  tkpack(actionInsertFrame, fill = "x", expand = TRUE)
  tkpack(actionFrame, fill = "x", expand = TRUE)

  messageFrame = tkframe(topFrame, relief = "raised", bd = 2)
  messageText = tktext(messageFrame, bg = "white", height = 5, width = 5)
  messageScr = tkscrollbar(messageFrame, command = function(...) tkyview(messageText, ...))
  tkconfigure(messageText, yscrollcommand = function(...) tkset(messageScr, ...))
  tkpack(messageText, side = "left", fill = "x", expand = TRUE)
  tkpack(messageScr, side = "right", fill = "y")
  tkpack(messageFrame, fill = "x")
  
  tkpack(topFrame)
  tkwm.focusmodel(tt, "active")
}

####################################################################
#- $Log: ABarrayGUI.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.1  2006/06/06 17:14:37  sunya
#- Changed ab1700 to ABarray in all releted files to be more compatible
#- in name convention to Bioconductor.
#-
#- Revision 1.3  2006/03/27 23:09:04  sunya
#- Added 'no imputation' option in GUI interface. Additional cosmetic changes.
#-
#- Revision 1.2  2006/03/14 19:48:30  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
