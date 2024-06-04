#' Empirical Error Type I Associated with a Bernoulli Distribution
#'
#'\code{ber_errorI} is used to obtain an empirical error type I when we use a random sample from a Bernoulli distribution.
#' @param c numeric, represents a positive value that defines a critical region. Default value is 1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the location parameter under the null hypothesis of a sample from a Bernoulli distribution. Default value is 0.5.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @return A list with number of replicates, sample size, and critical value that were used in the calculation of error type I
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Error type I when we use a random sample of size 200 from a Bernoulli distribution,
#' # a critical value c = 0.45 and R = 20000 to test H_0: theta = 0.7 vs H_1: theta != 0.7
#' ber_errorI(0.45,n=100,theta0=0.7,R=20000)
#' @importFrom stats rbinom dbinom
#' @export ber_errorI

ber_errorI <- function(c=1,n=150,theta0=0.5,R=12000){
  counter <- 0
  for(j in 1:R){
    muestra <- rbinom(n,size=1,prob=theta0)
    ln_ratio <- sum(dbinom(muestra,size=1,prob=theta0,log=TRUE) - dbinom(muestra,size=1,prob=mean(muestra),log=TRUE))
    if(ln_ratio < log(c)){
      counter <- counter + 1
    }
  }
  return(list(Rep=R,n=n,c=c,errorI_sim = counter/R))
}
