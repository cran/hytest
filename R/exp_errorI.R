#' Empirical Error Type I Associated with an Exponential Distribution
#'
#'\code{exp_errorI} is used to obtain an empirical error type I when we use a random sample from an Exponential distribution.
#' @param c numeric, represents a positive value that defines a critical region. Default value is 1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the location parameter under the null hypothesis of a sample from an Exponential distribution. Default value is 1.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @return A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Error type I when we use a random sample of size 200 from an Exponential distribution,
#' # a critical value c = 0.24 and R = 20000 to test H_0: theta = 2 vs H_1: theta != 2
#' exp_errorI(c=0.24,n=200,theta0=2,R=20000)
#' @importFrom stats rexp
#' @export exp_errorI

exp_errorI <- function(c=1,n=100,theta0=1,R=15000){
  counter <- 0
  for(j in 1:R){
    muestra <- rexp(n,1/theta0)
    ratio <- mean(muestra)/theta0
    Lambda_cal  <- (ratio^n)*exp(n*(1 - ratio))
    if(Lambda_cal < c){
      counter <- counter + 1
    }
  }
  return(list(Rep=R,n=n,c=c,errorI_sim = round(counter/R,5)))
}
