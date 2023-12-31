#- $Id: imputeFlag.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

imputeFlag = function(rawSig, pd = NULL, group = "", impute = "avg") {
  flagCount <- rep(0, dim(rawSig)[2])
  imputeValue <- list(flagCount = flagCount)

  rowImputed <- c()
  flagSet <- is.na(rawSig)    #- which flag
  flagCount <- apply(flagSet, 2, sum)
  imputeValue$flagCount <- flagCount
  
  #- knn impute Removed, since R 2.4, as it will not pass R check or build
  #- without a funcitonal impute package.
  ##if(impute == "knn") {
  ##  cat("\nPerform KNN imputation for probes with FLAG > 5000\n")
  ##  flush.console()

  ##  if( ! any(grep("impute", installed.packages()))) {
  ##    source("http://www.bioconductor.org/getBioC.R")
  ##    getBioC(pkgs = "impute")
  ##  }

  ##  library(impute)
  ##  imputed <- impute.knn(rawSig)
  ##  imputeValue$imputedData <- imputed$data
  ##  imputeValue$rowImputed <- NULL
  ##}
  if(impute == "knn") {
    cat("\nknn imputation is NOT implemented, using default subgroup average imputation\n")
  }
  ##else if(impute == "avg") {
    cat("Perform Average imputation for probes with FLAG > 5000\n")
    flush.console()

     #- Average impute
    grpMember <- levels(pd[, group])
    for(i in 1:dim(rawSig)[1]) {
      if(sum(flagSet[i,]) > 0) {
        for(j in 1:length(grpMember)) {
          memberHere <- which(pd[, group] == grpMember[j])
          if(sum(flagSet[i, memberHere], na.rm = TRUE) > 0) {
            flagged <- which(flagSet[i, memberHere])
            if(length(memberHere) - length(flagged) > 1) {
                                        #avg <- mean(rawSig[i, memberHere][-flagged])
              avg <- mean(rawSig[i, memberHere[-flagged]], na.rm = TRUE)
              if(! is.na(avg)) {
                rawSig[i, memberHere[flagged]] <- avg
                rowImputed <- c(rowImputed, i)
              }
            }
          }
        }
      }
    }
    imputeValue$imputedData <- rawSig
    imputeValue$rowImputed <- rowImputed    
  ##}
  #-write.table(rawSig, file = "ImputedValue.csv", sep = ",")
  return(imputeValue)
}

########################################################################
#- $Log: imputeFlag.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.3  2006/03/27 23:09:04  sunya
#- Added 'no imputation' option in GUI interface. Additional cosmetic changes.
#-
#- Revision 1.2  2006/03/20 23:00:16  sunya
#- Missing values (NA) are kept even after imputation. But most downstream
#- analysis will remove these NA values.
#- The hierarchical clustering probe list and expression values are now
#- exported into csv file.
#- Changed file naming convention. QC plot will begin QC_, and statistical
#- analysis plot will begin ST_ in the figures.
#- Values for control probes are saved into csv file.
#-
#- Revision 1.1  2005/10/25 23:48:56  sunya
#- Changed imputation occurance. Imputation now occurs after normalization.
#- Added function imputeFlag, previously in ab1700.R itself.
#-
