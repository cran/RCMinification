\name{ARlocpoly}
\alias{ARlocpoly}
\title{
Fit a nonlinear AR1 model using local polynomial regression
}
\description{
This function uses local polynomial regression to nonparametrically estimate 
the autoregression function in a nonlinear AR1 model.
}
\usage{
ARlocpoly(z, deg = 1, h, ...)
}
\arguments{
  \item{z}{numeric vector of time series observations.}
  \item{deg}{numeric, degree of local polynomial fit.}
  \item{h}{numeric, bandwidth for local polynomial estimate.}
  \item{...}{any other arguments taken by \code{locpoly}.}
}
\value{A list containing
  \item{x}{numeric vector of evaluation points.}
  \item{y}{numeric vector of nonparametric estimates at the values in \code{x}.}
  \item{h}{numeric, bandwidth}
}
\author{
L. Han and S. Snyman
}
\references{
Fan, J. and Yao, Q. (2008) Nonlinear Time Series: Nonparametric and 
Parametric Methods.  Springer. 
}
\examples{
x <- nonlinearAR1.sim(100, g = function(x) x*sin(x),  sd = 1.5) # simulated data 
ARlocpoly(x, deg = 0, h = 0.5)
}
\keyword{statistics}
