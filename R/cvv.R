"cvv" <-
function(data) {
   sd <- apply(data, 1, function(x) {sd(x, na.rm = TRUE)})
   mean <- rowMeans(data, na.rm = TRUE)
   cvValue <- sd / abs(mean)
   return(cvValue)
}

