#' Empirical Error Type I Associated with a Weibull Distribution
#'
#'\code{weibull_errorI} is used to obtain an empirical error type I when we use a random sample from a Weibull distribution.
#' @param c numeric, represents a positive value that defines a critical region. Default value is 1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the location parameter under the null hypothesis of a sample from a Weibull distribution. Default value is 1.
#' @param sigma numeric, represents the scale parameter of a Weibull distribution. It is assumed known and its default value is 1.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @return A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Error type I when we use a random sample of size 150 from a Weibull distribution,
#' # a critical value c = 0.5 and R = 500 to test H_0: theta = 3 vs H_1: theta != 3
#' library(gamlss.dist)
#' weibull_errorI(0.5,n=150,theta0=3,R=500)
#' @importFrom gamlss.dist rWEI dWEI WEI
#' @importFrom gamlss gamlss gamlss.control
#' @export weibull_errorI

weibull_errorI <- function(c=1,n=150,theta0=1,sigma=1,R=15000){
  counter <- 0
  for(j in 1:R){
    muestra <- rWEI(n,mu=theta0,sigma=sigma)
    result <- gamlss(muestra ~ 1,family = WEI,control=gamlss.control(trace = FALSE))
    mu_est <- exp(result$mu.coefficients)
    ln_ratio <- sum(dWEI(muestra,mu=theta0,sigma=sigma,log=TRUE) - dWEI(muestra,mu=mu_est,sigma=sigma,log=TRUE))
    if(ln_ratio < log(c)){
      counter <- counter + 1
    }
  }
  return(list(Rep=R,n=n,c=c,errorI_sim = counter/R))
}
