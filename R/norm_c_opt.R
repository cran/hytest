#' Critical Value Given a Nominal Error Type I Associated with a Normal Distribution
#'
#'\code{norm_c_opt} is used to obtain a critical value to achieve a nominal error type I when we use a random sample from a Normal distribution.
#' @param alpha numeric, represents a nominal error type I. Default value is 0.1.
#' @param n numeric, represents the size of the sample. Default value is 100.
#' @param theta0 numeric, represents the probability parameter under the null hypothesis of a sample from a Normal distribution. Default value is 0.5.
#' @param sd numeric, represents the scale parameter of a ]Normal distribution. It is assumed known and its default value is 1.
#' @param c1 numeric, represents a lower bound to the critical value. Default value is 1e-03.
#' @param c2 numeric, represents an upper bound to the critical value. Default value is 0.99.
#' @param R numeric, represents the number of replicates. Default value is 15000.
#' @param delta numeric, represents a precision parameter. Default value is 0.005.
#' @param tolerance numeric, represents a relative precision with respect a given alpha. Default value is 0.01.
#' @param max_iter integer, represents the maximum number of iterations. Default value is 100.
#' @return A list with number of replicates, sample size, nominal error type I, and empirical critical value obtained
#' associated with a likelihood ratio statistic.
#' @references Casella, G. and Berger, R. (2003). Statistical Inference, Second Edition. Duxbury Press.
#' @references Hogg, R., McKean, J., and Craig, A. (2019) Introduction to Mathematical Statistic.  Eighth edition. Pearson.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # Critical value when we use a random sample of size 100 from a Normal distribution
#' # given a nominal error type I equals to 0.1 and R = 10000
#' # to test H_0: theta = 0 vs H_1: theta != 0
#' norm_c_opt(alpha=0.1,n=100,theta0=0,sd=1,R=10000)
#' @export norm_c_opt

norm_c_opt <- function(alpha=0.1,n=100,theta0=0,sd=1,c1=1e-03,c2=0.999,R=15000,delta=0.005,tolerance=0.01,max_iter=100){
  k <- 1
  cond2 <- 1
  while(cond2 > tolerance & k <= max_iter){
    c_aux <- mean(c(c1,c2))
    alpha_c_aux <- norm_errorI(c_aux,n,theta0,sd,R)$errorI_sim
    cond1 <- sign(alpha - alpha_c_aux)
    if(cond1 > 0){
      c1 <- c_aux
    }
    else{
      c2 <- c_aux
    }
    cond2 <- round(abs(alpha_c_aux - alpha)/alpha,5)
    #print(paste('Iter ',k,': c value is ',round(c_aux,5),' and ','its simulated alpha is ',round(alpha_c_aux,5),sep=''))
    k <- k + 1
  }
  if(k <= max_iter){
    return(list(Rep = R, n = n, c_opt = c_aux, errorI_nom = alpha, errorI_sim = alpha_c_aux, cond = cond2, n_iter = k-1))
  }
  else{
    return('The goal was not achieved!. Please, try with other vector of parameters!')
  }
}

