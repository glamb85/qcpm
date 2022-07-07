#' @title get_info
#' @details
#' Internal function. \code{get_info} is called by \code{qcpm}.
#' @param data matrix or data frame containing the manifest variables.
#' @param inner A square (lower triangular) boolean matrix representing 
#' the inner model (i.e. the path relationships between latent variables).
#' @param outer list of vectors with column indices or column names
#' from \code{data} indicating the sets of manifest variables forming 
#' each block (i.e. which manifest variables correspond to each block).
#' @param modes character vector indicating the type of measurement for each
#' block. Possible values are: \code{'A', 'B'}. 
#' The length of \code{modes} must be equal to the length of \code{blocks}.
#' @param scheme string indicating the type of inner weighting
#' scheme. It is equal to \code{"factorial"} by default.
#'  Possible values are \code{"centroid"}, \code{"factorial"}.
#' @param tau if sepcifed indicate the specific quantile to be considered
#' @param tau_Alg is the vector of quantile specified by default. It is 
#' equal to (0.25,0.50,0.75).
#' @param fix.quantile is boolean equal to \code{FALSE} or \code{TRUE}  
#' indiciating whether invariance is holded.
#' @param \dots Further arguments passed on to \code{\link{get_info}}.
#' 
#' @return a string containing generalinformations of the inpunt and output 
#' parameters of the qcpm algorithm.
#' @keywords internal
#' @export
#'
get_info	<-	function(data, inner, outer,modes, scheme, tau,tau_Alg,fix.quantile,...)
{
  cat("\n")
  cat("\n")
  cat("---------------------------------------------------------------------------")
  cat("\n")
  cat("QC-PM: Quantile Composite-based Path Modeling","\n")
  cat("---------------------------------------------------------------------------")
  cat("\n")
  if(fix.quantile == TRUE) {
    cat("'fix.quantile == TRUE': tau is fixed to the median in the parameters 
    iterative estimation. Loadings, weights, and communality are admisible
    only for tau=0.5 ","\n")
  }    
  else{cat("Quantile estimation: complete","\n")}
  cat("\n")
  cat("Info models:" ,"\n")
  cat(paste("- Number of LV:", ncol(inner),"\n"))
  cat(paste("- Number of MV:", length(unlist(outer)),"\n"))
  cat("- LV modes:")
  cat(modes,"\n")
  cat(paste("- LV scheme: ", scheme,"\n"))
  cat("\n")
  cat("---------------------------------------------------------------------------")
  cat("\n")
  cat("Info quantile:","\n")
  if(is.null(tau) == TRUE ) {
    cat(paste("Quantile not specided and fixed by default:","\n"))
    cat(paste(tau_Alg),"\n")}
  else{cat(paste("Quantile selected:",tau_Alg[-length(tau_Alg)],"\n"))}
  cat("\n")
  cat("---------------------------------------------------------------------------")
  cat("\n")
  cat("Use also:","\n")
  cat("- for general results datails use function:             summary()","\n")
  cat("- for assessment results datails use function:          assessment()","\n")
  cat("- for boostrap results datails use function:            boot()","\n") 
  cat("- for reliability results of LV use function:           reliability()","\n") 
  cat("---------------------------------------------------------------------------")
  cat("\n")
  }


