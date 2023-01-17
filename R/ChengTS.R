ChengTS <- function(z, degree = 1, hopt, ...) {
    x <- z[-length(z)]
    y <- z[-1]
    Cheng <- function(x, y, degree=1, hopt, ...) {
        if (missing(hopt)) hopt <- dpill(x, y)
        hs <- hopt*seq(1, 2, length=11)
        xy.lp <- vector(length(hs), mode="list")
        Yvals <- matrix(0, nrow=401, ncol=length(hs))
        for (i in 1:length(hs)) {
            h <- hs[i]
            xy.lp[[i]] <- locpoly(x, y, degree = degree, bandwidth = h, ...)
            Yvals[,i] <- xy.lp[[i]]$y
        }
        H2 <- hs^2
        H2bar <- mean(H2)
        Ybar <-  (apply(Yvals, 1, mean))
        numer <- (H2 - H2bar)*(t(Yvals - Ybar))
        beta1hat <- apply(numer, 2, sum)/sum((H2 - H2bar)^2)
        beta0hat <- Ybar - beta1hat*H2bar
        return(list(x = xy.lp[[1]]$x, y = beta0hat))
    }
    out <- Cheng(x, y, degree=degree, hopt=hopt, ...)
    return(list(x = out$x, y = out$y))
}

