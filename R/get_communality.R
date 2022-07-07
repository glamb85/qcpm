#' @title get_communality
#' 
#' @details
#' Internal function. \code{get_communality} is called by \code{assessment} and 
#'
#' @param mv_name manifest variable label (string) 
#' @param lv_name latent variable label (string) 
#' @param tau  the quantile(s) to be estimated 
#' @param data  dataset. It includes manifest variable and latent score 
#' for a specific quantile 
#' @return the communality for each manifest varaible
#' @keywords internal
#' @export
#' 
get_communality <- function(mv_name, lv_name, tau, data){
  
  formula_mod <- paste(mv_name, "~", lv_name)
  formula_null <- paste(mv_name, "~", 1)
  
  mod_qr <- suppressWarnings(quantreg::rq(data = data, formula = formula_mod, 
                         tau = tau, method = "fn"))
  mod_null <- suppressWarnings(quantreg::rq(data = data, formula = formula_null, 
                           tau = tau, method = "fn"))
  
  rasw <- mod_qr$rho
  tasw <- mod_null$rho
  
  return(1 - rasw / tasw)
}

