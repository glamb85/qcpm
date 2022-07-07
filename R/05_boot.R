#' @title Inference on QC-PM model parameters (i.e., loadings and path coefficients)
#'
#' @description
#' \code{boot} returns in order the estimates, std. errors, t-values, 
#' p-values, and confidence interval at the specified confidence level 
#' for loadings and path coefficients for each quantile.  
#'
#' @details
#' 
#' The argument \code{qcpm} is an object of class qcpm returned by \code{qcpm} function. 
#' Std. errors are calculated by using the bootstrap method implemented in the 
#' \code{tidy.rq} function of the broom package  (Robinson, 2014). When \code{fix.quantile=TRUE}, 
#' the function boot returns only loading results for the quantile 0.5.
#' 
#'
#' @param qcpm is an object of class \code{qcpm}.
#' @param conf.level is the value used to fix the confidence level to use for the 
#' confidence interval. It is equal to 0.95 by default.
#' @param br specifies the number of bootstrap replications. It is fixed to 
#' \code{200} by default.
#' 
#' @return \item{boot.loadings}{the outer loading results for each considered quantile.}
#' @return \item{boot.path}{the path coefficient results for each considered quantile.}
#' 
#' @author Cristina Davino, Pasquale Dolce, Giuseppe Lamberti, Domenico Vistocco
#'
#' 
#' 
#' @references Davino, C., Dolce, P., Taralli, S. and Vistocco, D. (2020). Composite-based 
#' path modeling for conditional quantiles prediction. An application to assess 
#' health differences at local level in a well-being perspective.
#' \emph{Social Indicators Research}, doi:10.1007/s11205-020-02425-5.
#' 
#' @references Davino, C. and Esposito Vinzi, V. (2016). Quantile composite-based path modeling. 
#' \emph{Advances in Data Analysis and Classification}, \bold{10 (4)}, pp. 
#' 491--520, doi:10.1007/s11634-015-0231-9.
#' 
#' @references Dolce, P., Davino, C. and Vistocco, D. (2021). Quantile composite-based path modeling: 
#' algorithms, properties and applications. \emph{Advances in Data Analysis and Classification},
#' doi:10.1007/s11634-021-00469-0.
#' 
#' @references Robinson, D. (2014). broom: An R package for converting statistical analysis 
#' objects into tidy data frames. Available at 
#' \url{https://CRAN.R-project.org/package=broom}.
#' 
#' @seealso \code{\link{qcpm}}, \code{\link{assessment}}, \code{\link{summary}}, and 
#' \code{\link{reliability}}
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
#' # Define the model using laavan sintax. Use a set of regression formulas defining 
#' # firstly the structural model and then the measurement model
#' model <- "
# Structural model
#' ECOW ~ EDU
#' HEALTH ~ EDU + ECOW
#'
#' # Reflective measurement model
#' EDU =~ EDU1 + EDU2 + EDU3 + EDU4 + EDU5 + EDU6 + EDU7
#' ECOW =~ ECOW1 + ECOW2 + ECOW3 + ECOW4 + ECOW5 + ECOW6
#' HEALTH =~  HEALTH1 + HEALTH2 + HEALTH3
#' "
#'
#' 
#' # Apply qcpm
#' well.qcpm = qcpm(model,province)
#' well.boot = boot(well.qcpm)   
#' well.boot 
#'
boot <-  function(qcpm,conf.level=0.95,br=200)
{
  
  if (class(qcpm) != "qcpm") 
    stop("Argument 'qcpm' must be an object of class 'qcpm'")
  
  if(is.null(qcpm$model$tau)==TRUE) {tau=qcpm$model$tau_Alg} 
  else {tau =c(qcpm$model$tau,0.5)}
  
  Inference.path.summary=list()
  Inference.load.summary=list()
  
  for (i in 1: length(tau)){
    
    scores= qcpm$latent.scores[names(qcpm$latent.scores)==as.character(tau[i])]
    
    path_matrix = qcpm$model$qcpm.inner
    Y_lvs = as.matrix(scores[[1]])
    
    lvs_names = colnames(path_matrix)
    endogenous = as.logical(rowSums(path_matrix))
    num_endo = sum(endogenous)
    Path = low = upper = std = tv = pv = path_matrix
    
    for (aux in 1:num_endo) 
    {
      
      k1 <- which(endogenous)[aux]
      # index for indep LVs
      k2 = which(path_matrix[k1,] == 1)

      
      fit_score = suppressWarnings(quantreg::rq(Y_lvs[,k1] ~ Y_lvs[,k2],tau[i]))
      
      score_inf = as.data.frame(broom::tidy(fit_score,se.type = "boot",
                                            R=br,conf.int = TRUE,
                                            conf.level = conf.level ))
      
      Path[k1,k2] = score_inf[-1,2]
      std[k1,k2] = score_inf[-1,3]
      tv[k1,k2] = score_inf[-1,4]
      pv[k1,k2] = score_inf[-1,5]
      low[k1,k2] = score_inf[-1,6]
      upper[k1,k2] = score_inf[-1,7]
    }
    
    paths = as.matrix(Path[which(path_matrix==1)])
    low = as.matrix(low[which(path_matrix==1)])
    upper = as.matrix(upper[which(path_matrix==1)])
    
    std = as.matrix(std[which(path_matrix==1)])
    tv = as.matrix(tv[which(path_matrix==1)])
    pv = as.matrix(pv[which(path_matrix==1)])
    
    res_path = cbind(paths,std,tv,pv,low,upper)
    
    rownames(res_path) = get_element(Path)
    colnames(res_path) = c("Estimate", "Std. Error", "t value",  "Pr(>|t|)",
                           paste("low ",conf.level,"%",sep=""), paste("upper ",conf.level,"%",sep=""))
    
    Inference.path.summary[[length(Inference.path.summary)+1]]=
      path.coefficients = round(res_path,4)
    
    ################################################# loadings
    
    IDM = qcpm$model$qcpm.inner
    sets = qcpm$model$qcpm.outer.list
    Nom.MVs = get_names_MV( qcpm$model$qcpm.outer.matrix)
    
    lvs = ncol(IDM)
    mvs = length(unlist(sets))
    
    blocks = unlist(lapply(sets, length))
    load = load_low = load_up = load_std = load_tv = load_pv =  matrix(0,mvs,lvs)  
    
    aux = 0
    
    for (k in 1:lvs){
      for (m in (aux+1):sum(blocks[1:k])){
        
        fit_load= suppressWarnings(quantreg::rq(qcpm$data[,m]  ~  Y_lvs[, k],tau[i]))
        
        load_inf = as.data.frame(broom::tidy(fit_load,se.type = "boot",
                                             R=br,conf.int = TRUE,
                                             conf.level = conf.level))
        
        load[m,k]=load_inf[2,2]
        load_std[m,k]=load_inf[2,3]
        load_tv[m,k]=load_inf[2,4]
        load_pv[m,k]=load_inf[2,5]
        load_low[m,k]=load_inf[2,6]
        load_up[m,k]=load_inf[2,7]
        
      }
      aux = sum(blocks[1:k])
    }

    res_load = cbind(
      round(as.matrix(rowSums(load)),4),
      round(as.matrix(rowSums(load_std)),4),
      round(as.matrix(rowSums(load_tv)),4),
      round(as.matrix(rowSums(load_pv)),4),
      round(as.matrix(rowSums(load_low)),4),
      round(as.matrix(rowSums(load_up)),4)
    )
    
    rownames(res_load) = Nom.MVs
    colnames(res_load) = c("Estimate", "Std. Error", "t value",  "Pr(>|t|)",
                           paste("low ",conf.level ,"%",sep=""), paste("upper ",conf.level ,"%",sep=""))
    
    Inference.load.summary[[length(Inference.load.summary)+1]]=
      loadings = round(res_load,4)
    
  }
  
  names(Inference.path.summary) = tau
  names(Inference.load.summary) = tau
  
  cat("\n")
  cat("-----------------------------------------------------------------------")
  cat("\n")
  cat("QC-PM model inference: loadings and path coefficients significance","\n")
  cat("-----------------------------------------------------------------------")
  cat("\n")
  if(qcpm$model$fix.quantile == TRUE) {
    cat("'fix.quantile == TRUE' tau is fixed to the median in the parameters 
    iterative estimation. Loadings  are admisible only for tau=0.5 ","\n")
  }
  cat("\n")
  
  if(qcpm$model$fix.quantile == TRUE && is.null(qcpm$model$tau)==TRUE){
    
    list(boot.loadings=Inference.load.summary[2], boot.path = Inference.path.summary)
  }
  else if(qcpm$model$fix.quantile == TRUE && is.null(qcpm$model$tau)==FALSE){
    
    list(boot.loadings=Inference.load.summary[length(tau)], boot.path = Inference.path.summary[-length(tau)])
  }
  else if(qcpm$model$fix.quantile == FALSE && is.null(qcpm$model$tau)==FALSE){
    
    list(boot.loadings=Inference.load.summary[-length(tau)], boot.path = Inference.path.summary[-length(tau)])
  }
  else if(qcpm$model$fix.quantile == FALSE && is.null(qcpm$model$tau)==TRUE){
    
    list(boot.loadings=Inference.load.summary, boot.path = Inference.path.summary)
  }
  
  }
