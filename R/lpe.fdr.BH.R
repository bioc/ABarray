"lpe.fdr.BH" <-
function(lpe.result, adjp = "BH")
  {
	x.location <- grep("^x", names(lpe.result))
	y.location <- grep("^y", names(lpe.result))
   x <- lpe.result[, x.location]
   y <- lpe.result[, y.location]

   pnorm.diff <- pnorm(lpe.result$median.diff, mean = 0, 
            sd = lpe.result$pooled.std.dev)
   p.out <- 2 * apply(cbind(pnorm.diff, 1 - pnorm.diff), 
            1, min)
   #p.adj <- mt.rawp2adjp(p.out, proc = adjp)
   p.adj <- p.adjust(p.out, method = adjp)
        
	data.out <- data.frame(x = x, median.1 = lpe.result$median.1, 
            std.dev.1 = lpe.result$std.dev.1, y = y, median.2 = lpe.result$median.2, 
            std.dev.2 = lpe.result$std.dev.2, median.diff = lpe.result$median.diff, 
            pooled.std.dev = lpe.result$pooled.std.dev, abs.z.stats = abs(lpe.result$z.stats), 
            p.raw = p.out, p.adj = p.adj)

	aa <- cbind(data.out, z.real = data.out$abs.z.stats)
	#aa <- aa[order(aa[, 2], decreasing = TRUE), ]
	return(aa)
}

####################################
#- $Log$

