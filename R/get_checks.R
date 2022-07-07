#' @title intternal checks 
#' 
#' @details
#' Internal function. \code{get_checks} is called by \code{qcpm}.
#' @param data matrix or data frame containing the manifest variables.
#' @param inner A square (lower triangular) boolean matrix representing 
#' the inner model (i.e. the path relationships between latent variables).
#' @param outer list of vectors with column indices or column names
#' from \code{data} indicating the sets of manifest variables forming 
#' each block (i.e. which manifest variables correspond to each block).
#' @param scheme string indicating the type of inner weighting
#' scheme. It is equal to \code{"factorial"} by default.
#'  Possible values are \code{"centroid"} or  \code{"factorial"}.
#' @param tau if sepcifed indicates the specific quantile to be considered
#' @param \dots Further arguments passed on to \code{\link{get_checks}}.
#' 
#' @return A list containing checked  parameters for internal estimation of 
#' the qcpm algorithm.
#' @keywords internal
#' @export
#' 
get_checks  <- function (data,inner,outer,scheme,tau,... ){
  
  if(sum(is.na(data))>0)
    stop("\n'data' contains 'NA' please impute NA values before running 'qcpm'")
  
  if (class(data)[1]!= "matrix") {data = as.matrix(data) }
  
  if (is.matrix(data) && !is.numeric(data))
    stop("\nInvalid 'data' matrix. Must be a numeric matrix.")
  
  if (nrow(data) == 1)
    stop("\nCannot work with only one row in 'data'")
  
  if (ncol(data) == 1)
    stop("\nCannot work with only one column in 'data'")
  
  if (is.null(rownames(data))==TRUE)
    rownames(data) = 1:nrow(data)
  
  if (is.null(colnames(data))==TRUE) 
    colnames(data) = paste("MV", 1:ncol(data), sep="")
  
  ##########################
  
  #####################################################
  is_square_matrix <- function(x) {
    if (is.matrix(x)) {
      if (nrow(x) == ncol(x)) TRUE else FALSE      
    } else FALSE
  }
  
  is_lower_triangular <- function(x, diag = FALSE) {
    if (is.matrix(x)) {
      all(x[upper.tri(x, diag = diag)] == 0)
    } else FALSE
  }
  
  #####################################################
  
  if (class(inner)[1]!= "matrix" )
    stop("\n'inner' must be a matrix.")
  
  if (!is_square_matrix(inner))
    stop("\n'inner' must be a square matrix.")
  
  if (nrow(inner) == 1)
    stop("\n'inner' must have more than one row")
  
  if (!is_lower_triangular(inner))
    stop("\n'inner' must be a lower triangular matrix")
  
  
  for (j in 1:ncol(inner)) 
  {
    for (i in 1:nrow(inner)) 
    {
      if (length(intersect(inner[i,j], c(1,0))) == 0)
        stop("\nElements in 'inner' must be '1' or '0'")
    }      
  }
  
  if (is.null(rownames(inner))==TRUE) {
    LV_names = paste("LV", 1:ncol(inner), sep = "")
    dimnames(inner) = list(LV_names, LV_names)
  }
  if (is.null(rownames(inner))==FALSE && is.null(colnames(inner))==TRUE) {
    colnames(inner) = rownames(inner)
  }
  if (is.null(colnames(inner))==FALSE && is.null(rownames(inner))==TRUE) {
    rownames(inner) = colnames(inner)
  }
  if (!is.list(outer))
    stop("\n'outer' must be a list.")
  
  if (scheme!="factorial" && scheme!="centroid"){
    warning("Invalid argument 'scheme'. Default 'scheme=factorial' is used.")
    scheme <- "factorial"
  }
  
  if (length(outer) != nrow(inner))
    stop("\nNumber of rows in 'inner' different from length of 'outer'.")
  
  
  if (!is.numeric(tau)==TRUE && !is.null(tau)==TRUE ) 
    stop("\nArgument 'tau' must be a numeric value indicating the quantile")
  
  check.res=list(data=data, inner=inner, outer=outer, scheme=scheme, tau=tau)
  check.res
  
}
