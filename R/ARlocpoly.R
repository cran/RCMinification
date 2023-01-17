ARlocpoly <- function(z, deg = 1, h, ...) {
    n <- length(z)
    xn <- z[-1]
    xn1 <- z[-n]
    if (missing(h)) h <- dpill(xn1, xn)
    x.lp <- locpoly(xn1, xn, bandwidth = h, degree = deg, ...)
    x.lp[[3]] <- h; names(x.lp)[3] <- "h"
    return(x.lp)
}    
