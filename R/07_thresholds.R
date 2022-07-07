#' @title Admissible quantile thresholds
#'
#' @description
#' \code{thresholds} thresholds provides the maximum and minimum admissible quantile threshold.
#'
#' @details
#' The argument \code{x} is data frame that contains the manifest variables used to
#' estimate the qcpm models 
#' 
#'
#' @param x  is a data frame or a data matrix (statistical units x manifest variables).
#' 
#' @return A vector containing the maximum and minimum admisible quantile threshold values.
#' 
#' @author Cristina Davino, Pasquale Dolce, Giuseppe Lamberti, Domenico Vistocco
#'
#' 
#' 
#' @references Davino, C., Dolce, P., Taralli, S. and Vistocco, D. (2020) Composite-based 
#' path modeling for conditional quantiles prediction. An application to assess 
#' health differences at local level in a well-being perspective.
#' \emph{Social Indicators Research}, doi:10.1007/s11205-020-02425-5.
#' 
#' @references Davino, C. and Esposito Vinzi, V. (2016) Quantile composite-based path modeling. 
#' \emph{Advances in Data Analysis and Classification}, \bold{10 (4)}, pp. 
#' 491--520, doi:10.1007/s11634-015-0231-9.
#' 
#' @references Dolce, P., Davino, C. and Vistocco, D. (2021) Quantile composite-based path modeling: 
#' algorithms, properties and applications. \emph{Advances in Data Analysis and Classification},
#' doi:10.1007/s11634-021-00469-0.
#' 
#' 
#' @export
#' @examples
#' 
#' # Example of QC-PM in Well-Being analysis
#' # model with three LVs and reflective indicators
#' 
#' # load library and dataset province
#' library(qcpm)
#' data(province)
#' 
#' thresholds(province)
#' 
#'
thresholds <- function(x){
  
  max_threshold <- function(x){
    vec_max <- apply(x, 2, function(mv){
      tb <- prop.table(table(mv))
      return(tb[length(tb)])
    })
    return(round(1 - max(vec_max), 2))
  }

  min_threshold <- function(x){
    vec_min <- apply(x, 2, function(mv){
      tb <- prop.table(table(mv))
      return(tb[1])
    })
    return(as.numeric(substr(as.character(max(vec_min)), 1, 4)))
  }
  
  return(
  c(tau_min = min_threshold(x), 
    tau_max = max_threshold(x))
  )
}  
  
  
  
  