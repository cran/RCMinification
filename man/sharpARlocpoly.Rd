\name{sharpARlocpoly}
\alias{sharpARlocpoly}
\title{
Fit a nonlinear AR1 model using local polynomial regression and data sharpening
}
\description{
This function uses local polynomial regression to nonparametrically estimate 
the autoregression function in a nonlinear AR1 model, after employing 
data sharpening on the responses.
}
\usage{
sharpARlocpoly(z, deg = 1, h, ...)
}
\arguments{
  \item{z}{numeric vector of time series observations.}
  \item{deg}{numeric, degree of local polynomial fit.}
  \item{h}{numeric, bandwidth for local polynomial estimate.}
  \item{...}{any other arguments taken by \code{ARlocpoly}.}
}
\value{A list containing
  \item{x}{numeric vector of evaluation points.}
  \item{y}{numeric vector of nonparametric estimates at the values in \code{x}.}
}
\author{
L. Han and S. Snyman
}
\references{
Choi, E., Hall, P. and Rousson, V. (2000) Data Sharpening Methods
for Bias Reduction in Nonparametric Regression.  Annals of Statistics
28(5):1339-1355. 
}
\examples{
x <- nonlinearAR1.sim(100, g = function(x) x*sin(x),  sd = 1.5) # simulated data 
sharpARlocpoly(x, deg = 0, h = 0.5)
}
\keyword{statistics}
