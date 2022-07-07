#'@exportS3Method print qcpm
print.qcpm <- function(x, ...)
{
  cat("---------------------------------------------------------------------------")
  cat("\n")
  cat("QC-PM: Quantile  Composite-based Path Modeling", "\n")
  cat("---------------------------------------------------------------------------")
  cat("\n   NAME                ", "DESCRIPTION")
  cat("\n1  $outer.loadings     ", "outer loadings")
  cat("\n2  $outer.weights      ", "outer weights")
  cat("\n3  $path.coefficients  ", "path coefficients")
  cat("\n4  $latent.scores      ", "latent variable scores")
  cat("\n5  $data               ", "original datasets")
  cat("\n6  $model              ", "other model parameters")
  
  cat("\n")
  cat("---------------------------------------------------------------------------")
  cat("\nYou can also use the function 'summary'", "\n\n")    
  invisible(x)
}
