rET <- function(n, prob, rate) {
    B <- rbinom(n, 1, 1-prob)
    X <- rexp(n, rate)
    B*X
}

