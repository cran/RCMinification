\name{rRCMT}
\alias{rRCMT}
\title{
Tailed Exponential and Weibull Random Coefficient Minification Process Simulator
}
\description{
This function simulates sequences of tailed exponential and
Weibull random coefficient minification process variates. 
Random coefficients are lognormal distributed with parameters
mu and sigma. 
}
\usage{
rRCMT(n, p, p.delta, p.eps, lambda, gamma, mu, sigma, RCMTobj)
}
\arguments{
  \item{n}{number of observations.}
  \item{p}{power for transformation from exponential to Weibull.}
  \item{p.delta}{tailed exponential probability parameter when preceding observation is 0}
  \item{p.eps}{tailed exponential probability parameter when preceding observation is nonzero.}
  \item{lambda}{tailed exponential rate parameter when preceding observation is nonzero.}
  \item{gamma}{tailed exponential rate parameter when preceding observation is 0.} 
  \item{mu}{mu parameter for lognormal distribution used to simulated random coefficients.}
  \item{sigma}{sigma parameter for lognormal distribution used to simulate random coefficients.}
  \item{RCMTobj}{list containing elements n, p, p.delta, p.eps, lambda and gamma}
}
\author{
L. Han
}
\references{
Han, L., Braun, W.J. and Loeppky (2018) Random Coefficient Minification Processes. Statistical Papers, pp 1-22.  
}
\keyword{statistics}
