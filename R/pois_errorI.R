#' Empirical Error Type I Associated with a Poisson Distribution
#'
#'\code{pois_errorI} is used to obtain an empirical error type I when we use a random sample from a Poisson distribution.
#' @param c numeric, represents a positive value that defines a critical region. Default value is 1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the location parameter under the null hypothesis of a sample from a Poisson distribution. Default value is 1.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @return A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Error type I when we use a random sample of size 200 from an Poisson distribution,
#' # a critical value c = 0.85 and R = 20000 to test H_0: theta = 2 vs H_1: theta != 2
#' pois_errorI(0.85,n=100,theta0=2,R=20000)
#' @importFrom stats rpois dpois
#' @export pois_errorI

pois_errorI <- function(c=1,n=100,theta0=1,R=15000){
  counter <- 0
  for(j in 1:R){
    muestra <- rpois(n,theta0)
    ln_ratio <- sum(dpois(muestra,theta0,log=TRUE) - dpois(muestra,mean(muestra),log=TRUE))
    if(ln_ratio < log(c)){
      counter <- counter + 1
    }
  }
  return(list(Rep=R, n=n, c=c, theta0=theta0, errorI_sim = counter/R))
}
