RCMTmle <- function(y) {
p00hat <- function(y) { # p.delta
    if (p10hat(y) > 1e-6) {
        n <- length(y)
        zeros <- which(y[-n]==0)
        return(mean(y[zeros+1]==0))
    } else {
        return(1e-7)
    }
}

p10hat <- function(y) { # p.varepsilon
    n <- length(y)
    nonzeros <- which(y[-n] > 0)
    phat <- mean(y[nonzeros+1]==0)
    if (phat < 1e-7) {
       phat <- 1e-7
    }
    phat
}

powerestimator <- function(x, interval=c(0.2,2.5)) {
   x <- x[x > 0]
# this estimator assumes independence
    powerfun <- function(power, x) {
        power + mean(log(x)) - sum(x^(1/power)*log(x))/sum(x^(1/power))
    }
    power <- uniroot(powerfun, interval, x=x)$root
    list(power=power, rate=1/mean(x^(1/power)))
}

RCalphaEst5 <- function(y, P, muStart) {
    x <- y^(1/P)
    ACF <- acf(x[x>0], plot=FALSE)$acf[2:3]
    Z <- rnorm(1e+6)
    eX <- function(mu, sigma) {
        exp(sigma*Z + mu)
    }
    j <- 1:2
    dataset <- data.frame(j=j, ACF=ACF)
    MCACF <- function(mu, sigma, j) {
        a1 <- eX(mu, sigma)
        out1 <- mean(1/(a1+1))
        out2 <- out1^2
        c(out1, out2)
    }
    nls(ACF ~ MCACF(mu, 1, j), start = c(mu = muStart), data = dataset, control=nls.control(tol=1e-3), weights=c(2:1)^2)
}

logLikelihood <- function(x, lambda, p.eps, p.delta, gamma, f.RC, ...) {
    n <- length(x)
    logL <- 0
    for (j in 2:n) {
    logL <- logL +
        log(conditionaldensity(x[j], x[j-1], lambda,  p.eps, p.delta, gamma,
             f.RC, ...))
    }
    logL
}

integrand1 <- function(alpha, y, x, lambda, f.RC, ...) {
    exp(-lambda*alpha*y/(alpha+1))*f.RC(alpha, ...)
}

integrand2 <- function(alpha, y, x, lambda, f.RC, ...) {
    exp(-lambda*alpha*y/(alpha+1))*alpha/(alpha+1)*f.RC(alpha, ...)
}

conditionaldensity <- function(y, x, lambda, p.eps, p.delta, gamma, f.RC, ...) {
    if (x > 0 & y > 0) {
        cpdf <- lambda*integrate(function(alpha)
        integrand2(alpha, y, x, lambda, f.RC, ...), max((y/x)-1, 0), Inf)$value +
        integrand1(y/x-1, y, x, lambda, f.RC, ...)/x
        cpdf <- (1-p.eps)*cpdf
    } else {
        if (y > 0) {
            cpdf <- (1-p.delta)*gamma*exp(-gamma*y)
        } else {
            cpdf <- p.eps*(x > 0) + p.delta*(x==0)
        }
    }
    return(cpdf)
}
RCMTmle.p.mu.sigma <- function(y, p=.91, mu = -2, sigma = 1) {
    n <- length(y)
    p.eps <- p10hat(c(0,y))
    p.delta <- p00hat(c(0,y))
    ests <- list(n = n, p = p, p.eps = p.eps, p.delta = p.delta, mu = mu, sigma = sigma)
    xn <- y[-1]
    xn1 <- y[-n]
    if (any(xn1==0 & xn > 0)) {
        gamma <- 1/mean(xn[xn1==0 & xn > 0]^(1/p))
    } else {
        gamma <- 1/xn[1]
    }
    lambda <- 1/mean(y[y>0]^(1/p))
    par <- c(lambda, gamma)
    fun <- function(par, data, est) {
        -logLikelihood(data^(1/p), par[1],  est$p.eps, est$p.delta, par[2], f.RC=dlnorm, meanlog = est$mu, sdlog = est$sigma)
    }
    out <- optim(par=par, fun, lower=c(lambda/10, gamma/10), upper=c(lambda*5, gamma*5), method="L-BFGS-B", data = y, est = ests)
    ests$lambda <- out$par[1]; ests$gamma <- out$par[2]
    ests$like <- out$value; ests$y <- y
    ests$lambda <- out$par[1]; ests$gamma <- out$par[2]
    ests$like <- out$value; ests$y <- y
    ests
}
P <- powerestimator(y[y>0])$power
mu <- as.vector(coef(RCalphaEst5(y, P, 0)))
sigma <- 1
ests <- RCMTmle.p.mu.sigma(y,  P,  mu, sigma)
return(ests)
}
