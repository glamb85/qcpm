#' @title qc
#' 
#' @details
#' Internal function. \code{qc} is called by \code{get_internal_parameter_estimation} and 
#'
#' @param x the covariate variables 
#' @param y the dependent variable
#' @param tau  the quantile(s) to be estimated 
#' @return the qc correlation
#' @keywords internal
#' @export
#' 
qc <- function (x, y, tau) 
{
  n <- length(y)
  k <- length(tau)
  rho <- NULL
  for (i in 1:k) {
    ytau <- as.numeric(stats::quantile(y, probs = tau[i]))
    psiy <- rep(tau[i], n) - I(y - ytau < 0)
    rho[i] <- (x - mean(x)) %*% psiy/(n * sqrt(tau[i] - tau[i]^2) * 
                                        stats::sd(x))
  }
  return(list(tau = tau, rho = rho))
}