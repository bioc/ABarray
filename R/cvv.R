"cvv" <-
function(data) {
   sd <- apply(data, 1, function(x) {sd(x, na.rm = T)})
   mean <- rowMeans(data, na.rm = T)
   cvValue <- sd / abs(mean)
   return(cvValue)
}

