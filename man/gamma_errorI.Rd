% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma_errorI.R
\name{gamma_errorI}
\alias{gamma_errorI}
\title{Empirical Error Type I Associated with a Gamma Distribution}
\usage{
gamma_errorI(c = 1, n = 150, theta0 = 1, beta = 1, R = 15000)
}
\arguments{
\item{c}{numeric, represents a positive value that defines a critical region. Default value is 1.}

\item{n}{numeric, represents the size of the sample. Default value is 100.}

\item{theta0}{numeric, represents the shape parameter under the null hypothesis of a sample from a Gamma distribution. Default value is 1.}

\item{beta}{numeric, represents the scale parameter of a Gamma distribution. It is assumed known and its default value is 1.}

\item{R}{numeric, represents the number of replicates. Default value is 15000.}
}
\value{
A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
associated with a likelihood ratio statistic.
}
\description{
\code{gamma_errorI} is used to obtain an empirical error type I when we use a random sample from a Gamma distribution.
}
\examples{
# Error type I when we use a random sample of size 120 from a Gamma distribution,
# a critical value c = 0.5 and R = 200 to test H_0: theta = 1.5 vs H_1: theta != 1.5
gamma_errorI(0.5,n=120,theta0=1.5,R=200)
}
\references{
Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.

Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
