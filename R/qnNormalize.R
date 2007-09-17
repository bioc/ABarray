#- $Id: qnNormalize.R,v 1.1.1.1 2006/06/06 22:06:37 sunya Exp $

#- Perform quantile normalization. Adpated from limma package.  Part of the source
#- code is modified version of the limma package related function for
#- normalizeBetweenArrays.

"qnNormalize" <-
function (eData, snr, method = "quantile", snThresh = 3, ties = TRUE)
{
  eDim <- dim(eData)
  if (is.null(eDim)) 
    return(eData)
  if (eDim[2] == 1) 
    return(eData)

  theScale = rep(1, eDim[2])
  targetMedian = median(as.vector(eData), na.rm = T)
  
  ##- No normalization
  if(method == "none") {
    return(eData)
  }
  
    #- Perform quantile normalizaiton
  if(method == "quantile") {
    O <- S <- array(, eDim)
    if (ties) 
      R <- O
    nobs <- rep(eDim[1], eDim[2])
    i <- (0:(eDim[1] - 1))/(eDim[1] - 1)
    for (j in 1:eDim[2]) {
      Si <- sort(eData[, j], method = "quick", index.return = TRUE)
      if (ties) 
        R[, j] <- rank(eData[, j])
      nobsj <- length(Si$x)
      if (nobsj < eDim[1]) {
        nobs[j] <- nobsj
        isna <- is.na(eData[, j])
        S[, j] <- approx((0:(nobsj - 1))/(nobsj - 1), Si$x, 
                         i, ties = "ordered")$y
        O[!isna, j] <- ((1:eDim[1])[!isna])[Si$ix]
      }
      else {
        S[, j] <- Si$x
        O[, j] <- Si$ix
      }
    }
    m <- rowMeans(S)
    for (j in 1:eDim[2]) {
      if (nobs[j] < eDim[1]) {
        isna <- is.na(eData[, j])
        if (ties) 
          eData[!isna, j] <- approx(i, m, (R[!isna, j] - 1)/
                                    (nobs[j] - 1), ties = "ordered")$y
        else eData[O[!isna, j], j] <-
          approx(i, m, (0:(nobs[j] - 1))/(nobs[j] - 1), ties = "ordered")$y
      }
      else {
        if (ties) 
          eData[, j] <- approx(i, m, (R[, j] - 1)/(eDim[1] - 1), 
                               ties = "ordered")$y
        else eData[O[, j], j] <- m
      }
    }
    return(eData)
  }
  
  else if(method == "median") {
    theMedian = apply(eData, 2, median, na.rm = T)
    theScale = targetMedian / theMedian
  }
  else if(method == "mean") {
    theMean = colMeans(eData, na.rm = T)
    theScale = targetMedian / theMean
  }
  else if(method == "trimMean") {
    theMean = apply(eData, 2, function(x) {
      lowUp = quantile(x, probs = c(0.05, 0.95), na.rm = T)
      idxUse = intersect(which(x > lowUp[1]), which( x < lowUp[2]))
      mean(x[idxUse], na.rm = T)
    })
    theScale = targetMedian / theMean
  }
  else if(method == 'trimAMean') {
    if(missing(snr)) {
      print("S/N ratio not provided. Regular trimmed mean normalization will be provided.")
      theMean = apply(eData, 2, function(x) {
        lowUp = quantile(x, probs = c(0.05, 0.95), na.rm = T)
        idxUse = intersect(which(x > lowUp[1]), which( x < lowUp[2]))
        mean(x[idxUse], na.rm = T)
      })
      theScale = targetMedian / theMean
    }
    else {
      theMean = NULL
      for(i in 1:eDim[2]) {
        idxA = which(snr[, i] < snThresh)
        lowUp = quantile(eData[-idxA,i], probs = c(0.025, 0.975), na.rm = T)
        idxUse = intersect(which(eData[-idxA,i] > lowUp[1]), which(eData[-idxA, i] < lowUp[2]))
        theMean[i] = mean(eData[idxUse, i], na.rm = T)
      }
      theScale = targetMedian / theMean
    }
  } 
  
  retData = eData
  for(i in seq(along = theScale)) {
    retData[, i] = eData[, i] * theScale[i]
  }
  return(retData)
}

########################################################
#- $Log: qnNormalize.R,v $
#- Revision 1.1.1.1  2006/06/06 22:06:37  sunya
#- ABarray project converted from ab1700 project
#-
#- Revision 1.4  2006/04/17 19:04:56  sunya
#- Added vignette files. Added normalization methods for mean, median, trimmed mean.
#-
#- Revision 1.3  2006/03/27 23:09:04  sunya
#- Added 'no imputation' option in GUI interface. Additional cosmetic changes.
#-
#- Revision 1.2  2006/03/14 19:48:31  sunya
#- Changed icp (internal control probe) QC plots.
#- Added function for icp -> icpPlot
#- ANOVA analysis now performs probe filtering, but no FDR is calculated.
#- hclusterPlot now calculate correlation coefficient for probes, previously
#- it used Euclidian distance. The distance between arrays is still Euclidean.
#-
