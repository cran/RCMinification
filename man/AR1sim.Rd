\name{AR1sim}
\alias{AR1sim}
\title{Simulated Autoregressive Order 1 Time Series}
\usage{data(AR1sim)}
\description{
This is simulated AR(1) data with parameter 0.8.  The noise is normally distributed with variance 1.
}
\format{
  A numeric vector
}
\examples{
    acf(AR1sim)
    pacf(AR1sim)
    AR1sim.fit <- arima(AR1sim, order = c(1, 0, 0))
    AR1sim.fit
}
\keyword{datasets}
