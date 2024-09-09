#' Empirical Error Type I Associated with a Gamma Distribution
#'
#'\code{gamma_errorI} is used to obtain an empirical error type I when we use a random sample from a Gamma distribution.
#' @param c numeric, represents a positive value that defines a critical region. Default value is 1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the shape parameter under the null hypothesis of a sample from a Gamma distribution. Default value is 1.
#' @param beta numeric, represents the scale parameter of a Gamma distribution. It is assumed known and its default value is 1.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @return A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Error type I when we use a random sample of size 120 from a Gamma distribution,
#' # a critical value c = 0.5 and R = 200 to test H_0: theta = 1.5 vs H_1: theta != 1.5
#' gamma_errorI(0.5,n=120,theta0=1.5,R=200)
#' @importFrom stats rgamma dgamma glm Gamma coefficients
#' @export gamma_errorI

gamma_errorI <- function(c=1,n=150,theta0=1,beta=1,R=15000){
  counter <- 0
  for(j in 1:R){
    muestra <- rgamma(n,shape=theta0,scale=beta)
    shape_est <- coefficients(glm(muestra ~ 1,family = Gamma(link='identity')))
    ln_ratio <- sum(dgamma(muestra,shape=theta0,scale=beta,log=TRUE) - dgamma(muestra,shape=shape_est,scale=beta,log=TRUE))
    if(ln_ratio < log(c)){
      counter <- counter + 1
    }
  }
  return(list(Rep=R,n=n,c=c,errorI_sim = counter/R))
}
