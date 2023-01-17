nonlinearAR1.sim <- function(n, g, ...) { 
    eps <- rnorm(n, ...)
    x <- g(rnorm(1)) + eps[1]
    for (i in 2:n) {
        x[i] <- g(x[i-1]) + eps[i]
    }
    return(x)
}
