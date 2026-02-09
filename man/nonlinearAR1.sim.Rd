\name{nonlinearAR1.sim}
\alias{nonlinearAR1.sim}
\title{
Nonlinear AR1 Simulator
}
\description{
This function simulates sequences of variates follow a nonlinear autoregressive order 1 process
of the form z_n = g(z_n-1) + epsilon. A normal distribution is assumed for the innovations.
}
\usage{
nonlinearAR1.sim(n, g, ...)
}
\arguments{
  \item{n}{number of observations.}
  \item{g}{autoregression function.}
  \item{...}{any parameters that are taken by \code{rnorm}} 
}
\author{
L. Han and S. Snyman
} 
\examples{
    x <- nonlinearAR1.sim(50, g = function(x) x*sin(x), sd = 2.5)
    ts.plot(x)
}
\keyword{statistics}
