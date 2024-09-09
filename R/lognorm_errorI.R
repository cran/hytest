#' Empirical Error Type I Associated with a Log Normal Distribution
#'
#'\code{lognorm_errorI} is used to obtain an empirical error type I when we use a random sample from a Log Normal distribution.
#' @param c numeric, represents a positive value that defines a critical region. Default value is 1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the natural logarithm of location parameter under the null hypothesis of a sample from a Log Normal distribution. Default value is 0.
#' @param sdlog numeric, represents the natural logarithm of scale parameter of a Log normal distribution. It is assumed known and its default value is 1.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @return A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Error type I when we use a random sample of size 50 from an Log Normal distribution,
#' # a critical value c = 0.5 and R = 500 to test H_0: theta = 0 vs H_1: theta != 0
#' lognorm_errorI(c=0.5,n=50,theta0=0,sdlog=1,R=500)
#' @importFrom stats rlnorm dlnorm
#' @export lognorm_errorI

lognorm_errorI <- function(c,n=150,theta0=0,sdlog=1,R=15000){
  counter <- 0
  for(j in 1:R){
    muestra <- rlnorm(n,meanlog=theta0,sdlog)
    p_est <- mean(log(muestra))
    ln_ratio <- sum(dlnorm(muestra,meanlog=theta0,sdlog,log=TRUE) - dlnorm(muestra,meanlog=p_est,sdlog,log=TRUE))
    if(ln_ratio < log(c)){
      counter <- counter + 1
    }
  }
  return(list(Rep=R,n=n,c=c,errorI_sim = counter/R))
}
