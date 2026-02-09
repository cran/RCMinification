\name{AR1Typesim}
\alias{AR1Typesim}
\title{Simulated Autoregressive Order 1 Time Series with Random Coefficient}
\usage{data(AR1Typesim)}
\description{
This is simulated AR(1) data with uniformly distributed parameter on the interval (0, 1.6).  
The noise is normally distributed with variance 1.
}
\format{
  A numeric vector
}
\examples{
    acf(AR1Typesim)
    pacf(AR1Typesim)
    AR1Typesim.fit <- arima(AR1Typesim, order = c(1, 0, 0))
    AR1sim.fit <- arima(AR1sim, order = c(1, 0, 0))
    AR1Typesim.fit
    boxplot(resid(AR1sim.fit), resid(AR1Typesim.fit))
    qqnorm(resid(AR1sim.fit))
    qqnorm(resid(AR1Typesim.fit))
}
\keyword{datasets}
