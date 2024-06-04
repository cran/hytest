#' Empirical Error Type I Associated with a Normal Distribution
#'
#'\code{norm_errorI} is used to obtain an empirical error type I when we use a random sample from a Normal distribution.
#' @param c numeric, represents a positive value that defines a critical region. Default value is 1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the location parameter under the null hypothesis of a sample from a Normal distribution. Default value is 0.
#' @param sd numeric, represents the scale parameter of a Normal distribution. It is assumed known and its default value is 1.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @return A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Error type I when we use a random sample of size 70 from an Normal distribution,
#' # a critical value c = 0.65 and R = 20000 to test H_0: theta = 0 vs H_1: theta != 0
#' norm_errorI(0.65,70,theta0=0,sd=1,R=20000)
#' @importFrom stats rnorm
#' @export norm_errorI

norm_errorI <- function(c=1,n=100,theta0=0,sd=1,R=15000){
  counter <- 0
  for(j in 1:R){
    muestra <- rnorm(n,mean=theta0,sd)
    lambda_cal <- exp(-0.5*n*(mean(muestra) - theta0)^2)
    if(lambda_cal < c){
      counter <- counter + 1
    }
  }
  return(list(Rep=R, n=n, c=c, theta0=theta0, errorI_sim = counter/R))
}

