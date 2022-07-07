#' @exportS3Method summary qcpm
summary.qcpm <- function(object, ...){
  
  y = object
  cat("\n")
  cat("\n")
  cat("----------------------------------------------------------------------")
  cat("\n")
  cat("QC-PM: general results","\n")
  cat("----------------------------------------------------------------------")
  cat("\n")
  if(y$model$fix.quantile == TRUE) {
    cat("'fix.quantile == TRUE' tau is fixed to the median in the parameters 
    iterative estimation. Loadings and weights are admisible only for 
    tau=0.5 ","\n")
  }
  else{cat("quantile estimation: complete","\n")}
  cat("----------------------------------------------------------------------")
  cat("\n")
  cat("Results:","\n")
  cat(" - Outer Loading:", "\n")
  print(y$outer.loadings)
  cat("\n")
  cat("----------------------------------------------------------------------")
  cat("\n")
  cat(" - Outer weights:", "\n")
  print(y$outer.weights)
  cat("\n")
  cat("----------------------------------------------------------------------")
  cat("\n")
  cat(" - Path coefficients:", "\n")
  print(y$path.coefficients)
  cat("\n")
  cat("----------------------------------------------------------------------")
  
  }