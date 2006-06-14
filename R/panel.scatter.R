"panel.scatter" <-
function (x, y, col = "blue", bg = NA, pch = ".", 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
#    ok <- is.finite(x) & is.finite(y)
#    if (any(ok)) 
#        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
#            col = col.smooth, ...)

	 #abline(0, 1, col = col.smooth)
	 abline(1, 1, col = col.smooth)
	 abline(-1, 1, col = col.smooth)
}

