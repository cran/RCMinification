\name{rET}
\alias{rET}
\title{
Tailed Exponential Random Number Generator
}
\description{
This function simulates sequences of tailed exponential variates which
have survivor function P(X > x) = (1-p)exp(-lambda x), for x > 0 and 
P(X = 0) = p.   
}
\usage{
rET(n, prob, rate)
}
\arguments{
  \item{n}{number of observations.}
  \item{prob}{vector of probabilities.}
  \item{rate}{vector of exponential rate parameters.} 
}
\author{
L. Han
}
\references{
Littlejohn, R.P. (1994) A Reversibility Relationship for Two Markovian
Time Series Models with Stationary Exponential Tailed Distribution. Journal
of Applied Probability.  31 pp 575-581.  
}
\keyword{statistics}
