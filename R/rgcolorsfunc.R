"rgcolorsfunc" <-
function (n = 50) 
{
    k <- round(n/2)
    r <- c(rep(0, k), seq(0, 1, length = k))
    g <- c(rev(seq(0, 1, length = k)), rep(0, k))
	 b <- c(seq(0, 1, length = k), rev(seq(0, 1, length = k)))
    res <- rgb(r, g, b)
    #res <- rgb(r, g, rep(0, 2 * k))
    res
}

