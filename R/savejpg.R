"savejpg" <-
function(x, width = 1024, height = 768) {
  if(tolower(.Platform$OS.type) == "windows") {
    dev.copy(bmp, file = paste(x, ".bmp", sep = ""), width = width, height = height)
    dev.off()
  }
  else {
    dev.copy(jpeg, quality=100, file = paste(x, ".jpg", sep = ""), width = width, height = height)
    dev.off()
  }
}

