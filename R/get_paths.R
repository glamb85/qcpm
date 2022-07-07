#' @title get_paths
#' 
#' @details
#' Internal function. \code{get_paths} is called by \code{qcpm} and 
#'
#' @param path_matrix the matrix of path coefficients 
#' @param Y_lvs the matrix of latent variables  
#' @param tau  the quantile(s) to be estimated 
#' @return the path coefficients
#' @keywords internal
#' @export
#' 
get_paths <-  function(path_matrix, Y_lvs, tau, full=TRUE,...)
{
  
  lvs_names = colnames(path_matrix)
  endogenous = as.logical(rowSums(path_matrix))
  num_endo = sum(endogenous)
  Path = path_matrix
  
  for (aux in 1:num_endo) 
  {
    
    # index for endo LV
    k1 <- which(endogenous)[aux]
    # index for indep LVs
    k2 = which(path_matrix[k1,] == 1)
    # elimnate warnings
    path_rq= suppressWarnings(summary(quantreg::rq(Y_lvs[,k1] ~ Y_lvs[,k2],tau)))
    Path[k1,k2] = path_rq$coef[-1,1]
    
  }
  
  paths = as.matrix(Path[which(path_matrix==1)])
  rownames(paths) = get_element(Path)
  
  # output
  paths
}

