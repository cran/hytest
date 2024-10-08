% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gumbel_errorI.R
\name{gumbel_errorI}
\alias{gumbel_errorI}
\title{Empirical Error Type I Associated with a Gumbel Distribution}
\usage{
gumbel_errorI(c = 1, n = 150, theta0 = 0, sigma = 1, R = 15000)
}
\arguments{
\item{c}{numeric, represents a positive value that defines a critical region. Default value is 1.}

\item{n}{numeric, represents the size of the sample. Default value is 100.}

\item{theta0}{numeric, represents the location parameter under the null hypothesis of a sample from a Gumbel distribution. Default value is 0.}

\item{sigma}{numeric, represents the scale parameter of a Gumbel distribution. It is assumed known and its default value is 1.}

\item{R}{numeric, represents the number of replicates. Default value is 15000.}
}
\value{
A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
associated with a likelihood ratio statistic.
}
\description{
\code{gumbel_errorI} is used to obtain an empirical error type I when we use a random sample from a Gumbel distribution.
}
\examples{
# Error type I when we use a random sample of size 150 from a Gumbel distribution,
# a critical value c = 0.5 and R = 500 to test H_0: theta = 3 vs H_1: theta != 3
library(gamlss.dist)
gumbel_errorI(0.5,n=150,theta0=3,R=500)
}
\references{
Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.

Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
