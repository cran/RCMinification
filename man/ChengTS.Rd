\name{ChengTS}
\alias{ChengTS}
\title{
Fit a nonlinear AR1 model using local polynomial regression via 
the method of Cheng et al.  
}
\description{
This function uses local polynomial regression to nonparametrically estimate 
the autoregression function in a nonlinear AR1 model using Cheng's bias
reduction method.
}
\usage{
ChengTS(z, degree = 1, hopt, ...)
}
\arguments{
  \item{z}{numeric vector of time series observations.}
  \item{degree}{numeric, degree of local polynomial fit.}
  \item{hopt}{numeric, base bandwidth for local polynomial estimate.}
  \item{...}{any other arguments taken by \code{locpoly}.}
}
\value{A list containing
  \item{x}{numeric vector of evaluation points.}
  \item{y}{numeric vector of nonparametric estimates at the values in \code{x}.}
}
\author{
L. Han and S. Snyman
}
\references{
Cheng, M., Huang, R., Liu, P. and Liu, H. (2018) Bias reduction
for nonparametric and semiparametric regression models.  Statistica
Sinica 28(4):2749-2770.  
}
\examples{
x <- nonlinearAR1.sim(100, g = function(x) x*sin(x),  sd = 1.5) # simulated data 
ChengTS(x, degree = 1, hopt = 0.5)

x <- nonlinearAR1.sim(100, g = function(x) x*sin(x),  sd = 0.5) # simulated data 
degree <- 1; xrange <- c(-.5, .5); n <- length(x)
h <- thumbBw(x[-n], x[-1], deg = degree, kernel=gaussK)
x.lp <- ARlocpoly(x, deg = degree, h = h, range.x = xrange)
x.shp <- sharpARlocpoly(x, deg = degree, range.x = xrange, h = x.lp$h*n^(4/45))
x.cheng <- ChengTS(x, degree = degree,  hopt = h, range.x = xrange)
lag.plot(x, do.lines=FALSE)
lines(x.lp)
lines(x.shp, col=2)
lines(x.cheng, col=4)
}
\keyword{statistics}
