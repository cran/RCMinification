rRCMT <- function(n, p, p.delta, p.eps, lambda, gamma, mu, sigma, RCMTobj) {
# simulator of random coefficient minification process based on 
# exponential or Weibull tail distributions (as defined in Littlejohn, 
# 1994)
# Random coefficients are lognormal with parameters mu and sigma
    if (missing(n)) {
         p <- RCMTobj$p
         p.delta <- RCMTobj$p.delta
         p.eps <- RCMTobj$p.eps
         lambda <- RCMTobj$lambda
         gamma <- RCMTobj$gamma
         mu <- RCMTobj$mu
         sigma <- RCMTobj$sigma
    }
    x <- 0
    y <- as.numeric(n)
    if (x==0) {
        y[1] <- rET(1, p.delta, gamma)
    } else {
        z <- rET(1, p.eps, lambda)
        alpha <- rlnorm(1, meanlog=mu, sdlog=sigma)
        k <- 1+1/alpha
        y[1] <- min(z, alpha*x)*k
    }
    for (i in 2:n) {
        if (y[i-1]==0) {
            y[i] <- rET(1, p.delta, gamma)
        } else {
            z <- rET(1, p.eps, lambda)
            alpha <- rlnorm(1, meanlog=mu, sdlog=sigma)
            k <- 1+1/alpha
            y[i] <- min(z, alpha*y[i-1])*k
        }
    }  
    y^p
}


