#' @title Assessment measures of quantile composite-based path modeling
#'
#' @description
#' \code{assessment} returns the following measures for assessing both the inner 
#' and the outer model: communality of each manifest variable, communality of 
#' each block,redundancy of each manifest variable of endogenous blocks, redundancy 
#' of the endogenous blocks, and \eqn{pseudo-R^2} for each inner equation.
#'
#' @details
#' All the assessment measures discussed in Davino et al. (2016) and Dolce et al. (2021) 
#' are based on \eqn{pseudo-R^2}, proposed by Koenker and Machado (1999), which simulates the 
#' role and interpretation of the \eqn{R^2} in classical regression analysis. The \eqn{pseudo-R^2} is 
#' considered as a local measure of goodness of fit for a particular quantile as it measures 
#' the contribute of the selected regressors to the explanation of the dependent variable, 
#' with respect to the trivial model without regressors. In more technical way, \eqn{pseudo-R^2}
#' compares the residual absolute sum of weighted differences using the selected model with 
#' the total absolute sum of weighted differences using a model with the only intercept. 
#' The \eqn{pseudo-R^2} can be used to assess the inner model measuring the amount of variability of a 
#' given endogenous construct explained by its explanatory constructs. A synthesis of the 
#' evaluations regarding the whole inner model can be obtained by the average of all the \eqn{pseudo-R^2}. 
#' Communality indicates how much of the MV variance is explained by the corresponding construct. 
#' It can be calculated for each MV, and for each block, using the average of MV communalities. 
#' Redundancy measures the percent of the variance of MVs in an endogenous block that is predicted 
#' from the explanatory constructs related to the endogenous construct. Redundancy can be computed 
#' only for each MVs of endogenous blocks and for the whole endogenous blocks, using the average of 
#' MV redundancies. Results are provided for each quantile of interest. When \code{fix.quantile=TRUE}, the 
#' function returns communalities and redundancies only for the quantile 0.5.
#' 
#' 
#' 
#' @param qcpm is an object of class \code{qcpm}
#' 
#' @return \item{Communality}{Communality of each MV. It is the proportion of the MV 
#' variance explained by the corresponding construct.}
#' @return \item{Block_Communality}{Communality of a whole block. It is computed 
#' as average of the MV communalities belonging to that block.}
#' @return \item{Redundancy}{Redundancy of each MV of the endogenous blocks. It measures 
#' the percent of the variance of MVs in endogenous blocks that is predicted from the 
#' explanatory constructs related to the endogenous construct.}
#' @return \item{Block_Redundancy}{Redundancy of a block. It is computed as average of 
#' MV redundancies belonging to that block.}
#' @return \item{pseudo.R2}{The \eqn{pseudo-R^2}. It assesses the goodness of fit of the inner 
#' model.}
#' 
#' @author Cristina Davino, Pasquale Dolce, Giuseppe Lamberti, Domenico Vistocco
#'
#' @references Davino, C., Dolce, P., Taralli, S. and Vistocco, D. (2020). Composite-based 
#' path modeling for conditional quantiles prediction. An application to assess 
#' health differences at local level in a well-being perspective.
#' \emph{Social Indicators Research}, doi:10.1007/s11205-020-02425-5..
#' 
#' @references Davino, C. and Esposito Vinzi, V. (2016). Quantile composite-based path modeling. 
#' \emph{Advances in Data Analysis and Classification}, \bold{10 (4)}, pp. 
#' 491--520, doi:10.1007/s11634-015-0231-9.
#' 
#' @references Davino, C., Esposito Vinzi, V. and Dolce, P. (2016). Assessment and validation in 
#' quantile composite-based path modeling. In: Abdi H.,  Esposito Vinzi, V., Russolillo, G., 
#' Saporta, G., Trinchera, L. (eds.). \emph{The Multiple Facets of Partial Least Squares Methods}, 
#' chapter 13. Springer proceedings in mathematics and statistics. Springer, Berlin
#' 
#' @references Dolce, P., Davino, C. and Vistocco, D. (2021). Quantile composite-based path modeling: 
#' algorithms, properties and applications. \emph{Advances in Data Analysis and Classification},
#' doi:10.1007/s11634-021-00469-0.
#' 
#' @references Koenker, R. and Machado, J.A. (1999). Goodness of fit and related inference processes 
#' for quantile regression. \emph{Journal of the American Statistical Association}, \bold{94 (448)} 
#' pp. 1296--1310, doi: 10.1080/01621459.1999.10473882
#' 
#' @references He, X.M. and Zhu, L.X. (2003). A lack-of-fit test for quantile regression. 
#' \emph{Journal of the American Statistical Association} \bold{98} pp. 1013--1022, 
#' doi: 10.1198/016214503000000963
#' 
#' 
#' 
#' @seealso \code{\link{summary}}, \code{\link{qcpm}}, \code{\link{boot}}, and 
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
#' well.assessment = assessment(well.qcpm)   
#' well.assessment 
#'
assessment = function (qcpm){
  
  if (class(qcpm) != "qcpm") 
    stop("Argument 'qcpm' must be an object of class 'qcpm'")
  
  if(is.null(qcpm$model$tau)==TRUE) {tau_outer=qcpm$model$tau_Alg} 
  else {tau_outer =c(qcpm$model$tau,0.5)}
  
  inner = qcpm$model$qcpm.inner
  outer = qcpm$model$qcpm.outer.matrix
  
  Communality = matrix(NA,dim(outer)[2],length(tau_outer))
  Block_Communality = matrix(NA,length(tau_outer),dim(inner)[2])
  Redundancy = matrix(NA,sum(rowSums(outer)[as.logical(rowSums(inner))]), length(tau_outer))
  Block_Redundancy =matrix(NA,length(tau_outer),sum(as.logical(rowSums(inner))))
  pseudo.r2 = matrix(NA,length(tau_outer),sum(as.logical(rowSums(inner))))
  
  for (k in 1: length(tau_outer)){
    comm_mv=NULL  
    
    scores= qcpm$latent.scores[names(qcpm$latent.scores)==as.character(tau_outer[k])][[1]]
    data.ass = cbind(qcpm$data, scores)
    out.mv.names=list()
    ################################################################ communalities
    for (i in 1: dim(inner)[1]){
      out.mv.names[[length(out.mv.names)+1]] =  names(which(outer[i,]==1))
    }
    
    for (j in 1: length(out.mv.names)){
      
      comm_xi=NULL
      
      for (i in 1: length(out.mv.names[[j]])){
        
        
        
        comm_xi[i]= get_communality(mv_name =out.mv.names[[j]][i], 
                                    lv_name =rownames(outer)[j], 
                                    data = data.ass, tau=tau_outer[k])
      }
      comm_mv=rbind(comm_mv,as.matrix(comm_xi)) 
      
      Block_Communality[k,j] = mean(comm_xi) 
    }
    Communality[,k]=comm_mv  
  }
  
  ################################################################ redundancy
  for (k in 1: length(tau_outer)){
    
    scores= qcpm$latent.scores[names(qcpm$latent.scores)==as.character(tau_outer[k])][[1]]
    data.ass = cbind(qcpm$data, scores)
    
    Y_lvs = as.matrix(scores)
    lvs_names = colnames(inner)
    endogenous = as.logical(rowSums(inner))
    num_endo = sum(endogenous)
    pseudo.r2.names = NULL
    red.names = NULL
    ra.r2=NULL
    Red_MV = NULL
    
    for (aux in 1:num_endo) 
    {
      # index for endo LV
      k1 <- which(endogenous)[aux]
      # index for indep LVs
      k2 = which(inner[k1,] == 1)
      lv_pred = names(which(inner[k1,] == 1))
      
      rasw = suppressWarnings(quantreg::rq(Y_lvs[,k1] ~ Y_lvs[,k2],tau=tau_outer[k])$rho)
      tasw = suppressWarnings(quantreg::rq(Y_lvs[,k1] ~ 1,tau=tau_outer[k])$rho)
      
      ra.r2 = (1 - rasw / tasw)
      
      pseudo.r2.names[aux] = colnames(Y_lvs)[k1]
      pseudo.r2[k,aux] = ra.r2
      
      mv_names = names(which(outer[k1,]==1))
      lv_name = rownames(outer)[k1]
      
      red_xi=NULL
      for(j in 1: length(mv_names)){
        
        red_xi[j]= get_communality(mv_name =mv_names[j], 
                                   lv_name =lv_name,
                                   data = data.ass, tau=tau_outer[k])*ra.r2
      }
      red.names[aux] = colnames(Y_lvs)[k1]
      Block_Redundancy[k,aux] = mean(red_xi)
      
      red.mat = as.matrix(red_xi)
      rownames(red.mat) = paste(lv_name,mv_names,sep="-")
      
      Red_MV=rbind(Red_MV,red.mat)
      Red_MV
      
    }
    
    Redundancy[,k]= Red_MV
    rownames(Redundancy) = rownames(Red_MV)
    
  }
  
  
  cat("\n")
  cat("----------------------------------------------------------------------")
  cat("\n")
  cat("QC-PM model assesement: Cummunality, ridondance, and pseudo-R2","\n")
  cat("----------------------------------------------------------------------")
  cat("\n")
  if(qcpm$model$fix.quantile == TRUE) {
    cat("'fix.quantile == TRUE' tau is fixed to the median in the parameters 
        iterative estimation. Cummunality are admisible only for tau=0.5 ","\n")
  }
  cat("\n")
  
  if(qcpm$model$fix.quantile == TRUE && is.null(qcpm$model$tau)==TRUE){
    
    Communality = as.matrix(Communality[,2])
    Block_Communality = as.matrix(Block_Communality[2,])
    Redundancy = as.matrix(Redundancy)
    Block_Redundancy = as.matrix(Block_Redundancy)
    
    colnames(Redundancy) = tau_outer
    colnames(Block_Redundancy) = pseudo.r2.names
    rownames(Block_Redundancy) =tau_outer
    
    colnames(Communality) = colnames(Block_Communality) = 0.5
    rownames(Communality) = colnames(outer)
    rownames(Block_Communality) = colnames(inner)
    
    pseudo.R2 = pseudo.r2
    colnames(pseudo.R2) =  pseudo.r2.names
    rownames(pseudo.R2) = tau_outer
    
  }
  
  else if(qcpm$model$fix.quantile == TRUE && is.null(qcpm$model$tau)==FALSE){
    
    
    Communality = as.matrix(Communality[,length(tau_outer)])
    Block_Communality = as.matrix(Block_Communality[length(tau_outer),])
    Redundancy = as.matrix(Redundancy[,-length(tau_outer)])
    if(length(qcpm$model$tau)>1){
    Block_Redundancy = as.matrix(Block_Redundancy[-length(tau_outer),])}
    else{
      Block_Redundancy = t(as.matrix(Block_Redundancy[-length(tau_outer),]))}
    
    
    colnames(Redundancy) = rownames(Block_Redundancy) = qcpm$model$tau
    colnames(Block_Redundancy) =pseudo.r2.names
    
    colnames(Communality) = colnames(Block_Communality) = 0.5
    rownames(Communality) = colnames(outer)
    rownames(Block_Communality) = colnames(inner)
    
    
    if(length(qcpm$model$tau)>1){
      pseudo.R2 = as.matrix(pseudo.r2[-length(tau_outer),])}
    else{
      pseudo.R2 = t(as.matrix(pseudo.r2[-length(tau_outer),]))}
    
    rownames(pseudo.R2) =  qcpm$model$tau
    colnames(pseudo.R2) = pseudo.r2.names
    
  }
  
  else if(qcpm$model$fix.quantile == FALSE && is.null(qcpm$model$tau)==FALSE){
    
    Communality = as.matrix(Communality[,-length(tau_outer)])
    
    if(length(qcpm$model$tau)>1){
      Block_Communality = t(as.matrix(Block_Communality[-length(tau_outer),]))}
    else{
      Block_Communality = as.matrix(Block_Communality[-length(tau_outer),])}
    
    Redundancy = as.matrix(Redundancy[,-length(tau_outer)])
    if(length(qcpm$model$tau)>1){
      Block_Redundancy = as.matrix(Block_Redundancy[-length(tau_outer),])}
    else{
      Block_Redundancy = t(as.matrix(Block_Redundancy[-length(tau_outer),]))}
    
    colnames(Redundancy) = rownames(Block_Redundancy) = qcpm$model$tau
    colnames(Block_Redundancy) =pseudo.r2.names
    
    colnames(Communality) = colnames(Block_Communality) = qcpm$model$tau
    rownames(Communality) = colnames(outer)
    rownames(Block_Communality) = colnames(inner)
    
    if(length(qcpm$model$tau)>1){
      pseudo.R2 = as.matrix(pseudo.r2[-length(tau_outer),])}
    else{
      pseudo.R2 = t(as.matrix(pseudo.r2[-length(tau_outer),]))}
    
    rownames(pseudo.R2) =  qcpm$model$tau
    colnames(pseudo.R2) = pseudo.r2.names
    
    
  }
  
  else if(qcpm$model$fix.quantile == FALSE && is.null(qcpm$model$tau)==TRUE){
    
    colnames(Communality) =  tau_outer
    rownames(Communality) = get_names_MV(outer)
    rownames(Block_Communality) = tau_outer
    colnames(Block_Communality) = colnames(inner)
    colnames(Redundancy) = tau_outer
    colnames(Block_Redundancy) = colnames(pseudo.r2) =  pseudo.r2.names
    rownames(Block_Redundancy) = rownames(pseudo.r2) = tau_outer
    pseudo.R2 = pseudo.r2
    
  }
  
  
  assessment.results = list(Communality=round(Communality,4),
                            Block_Communality=round(Block_Communality,4),
                            Redundancy = round(Redundancy,4),
                            Block_Redundancy=round(Block_Redundancy,4),
                            pseudo.R2=round(pseudo.R2,4))
  
  
  return(assessment.results)
}

