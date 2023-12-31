"panel.cor" <-
function(x, y, digits=3, prefix="", cex.cor)
	{
     usr <- par("usr"); on.exit(par(usr))
     par(usr = c(0, 1, 0, 1))
     r <- abs(cor(x, y, use = "complete.obs"))
     txt <- format(c(r, 0.123456789), digits=digits)[1]
     txt <- paste(prefix, txt, sep="")
     if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
     text(0.5, 0.5, txt, cex = cex.cor)
   }

