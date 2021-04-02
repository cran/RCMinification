\name{RCMTmle}
\alias{RCMTmle}
\title{
Tailed Exponential and Weibull Random Coefficient Minification Maximum Likelihood Estimation
}
\description{
This function estimates parameters for 
 tailed exponential and
Weibull random coefficient minification process models from a 
nonnegative time series. 
}
\usage{
RCMTmle(y)
}
\arguments{
  \item{y}{numeric vector of nonnegative observations.}
}
\value{A list containing
  \item{n}{the number of time series observations.}
  \item{p}{estimated power for transformation from exponential to Weibull.}
  \item{p.eps}{estimated tailed exponential probability parameter when preceding observation is nonzero.}
  \item{p.delta}{estimated tailed exponential probability parameter when preceding observation is 0}
  \item{mu}{estimated mu parameter for lognormal distribution used to simulated random coefficients.}
  \item{sigma}{estimated sigma parameter for lognormal distribution used to simulate random coefficients.}
  \item{lambda}{estimated tailed exponential rate parameter when preceding observation is nonzero.}
  \item{gamma}{estimated tailed exponential rate parameter when preceding observation is 0.} 
  \item{like}{maximum value of likelihood.}
  \item{y}{original observations}
}
\author{
L. Han
}
\references{
Han, L., Braun, W.J. and Loeppky (2018) Random Coefficient Minification Processes. Statistical Papers, pp 1-22.  
}
\keyword{statistics}
