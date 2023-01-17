sharpARlocpoly <- function(z, deg = 1, h, ...) {
    z.lp <- ARlocpoly(z, deg=deg, h, ...)
    n <- length(z)
    fittedvalues <- approx(z.lp$x, z.lp$y, xout=z[-n], rule=2)$y
    residuals <- z[-1] - fittedvalues
    zsharp <- z[-1] + residuals
    if (missing(h)) h <- z.lp$h
    z.shp <- locpoly(z[-n], zsharp, bandwidth = h, degree = deg, ...)
    return(z.shp)
}

