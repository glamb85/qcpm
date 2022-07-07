#' @title get_loadings
#' 
#' @details
#' Internal function. \code{get_loadings} is called by \code{qcpm} and 
#'
#' @param data the matrix of data (manifest variables)
#' @param sets outer model
#' @param mvs number of manifest variables
#' @param lvs number of latent variables
#' @param IDM the path matrix
#' @param tau  the quantile(s) to be estimated 
#' @param LV  the  estimated latent variables
#' @param qcorr  boolean. If it si equal to \code{TRUE},
#' loadings are standardized.
#' @return the loadings estimated for each latent variables
#' @keywords internal
#' @export
#' 
get_loadings <- function(data, sets, mvs,lvs, IDM, tau, LV, qcorr,... ){
  
  Beta0 = matrix(0,mvs,lvs)  
  mvs= length(unlist(sets))
  lvs=ncol(IDM)
  
  blocks = unlist(lapply(sets, length))
  ODMload = matrix(0,mvs,lvs)  
  aux = 0
  if (qcorr == TRUE){
    for (k in 1:lvs){
      for (m in (aux+1):sum(blocks[1:k])){
        
        ODMload[m,k]=suppressWarnings(qc(x = LV[, k], y = data[,m], tau = tau)$rho)
      }
      aux = sum(blocks[1:k])
    }
  }
  if (qcorr == FALSE){
    for (k in 1:lvs){
      for (m in (aux+1):sum(blocks[1:k])){
        Beta0[m,k]=suppressWarnings(quantreg::rq(data[,m]  ~ LV[, k], tau)$coefficients[1])
        ODMload[m,k]=suppressWarnings(quantreg::rq(data[,m]  ~  LV[, k], tau)$coefficients[2])
        
        #Beta0[m,k]=quantreg::rq(scale(data[,m])  ~ LV[, k], tau)$coefficients[1]
        #ODMload[m,k]=quantreg::rq(scale(data[,m])  ~  LV[, k], tau)$coefficients[2]
      }
      aux = sum(blocks[1:k])
    }
  }
  colnames(ODMload)=colnames(IDM)
  # output
  ODMload
}