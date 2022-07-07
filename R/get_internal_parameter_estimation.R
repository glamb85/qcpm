#' @title get_internal_parameter_estimation
#' 
#' @details
#' Internal function. \code{get_internal_parameter_estimation} is called by \code{qcpm} and 
#'
#' @param X matrix or data frame containing the manifest variables.
#' @param modes character vector indicating the type of measurement for each
#' block. Possible values are: \code{"A", "B"}. 
#' The length of \code{modes} must be equal to the length of \code{blocks}.
#' @param tau if sepcifed indicate the specific quantile to be considered
#' @param lvs  the number of latent variables
#' @param lvs.names  the label of latent variables
#' @param IDM  the path matrix
#' @param sets  the outer model
#' @param scheme  the internal scheme
#' @param method  rq method. It is equal to \code{"br"} by default
#' @param fix.quantile is boolean equal to \code{FALSE} or \code{TRUE} indicating 
#' whether invariance is holded.
#' @param tol decimal value indicating the tolerance criterion for the iterations (tol=0.00001).
#' @param maxiter integer indicating the maximum number of iterations (maxiter=100 by default).
#'
#' @return the outer weights
#' 
#' @keywords internal
#' @export
#'
get_internal_parameter_estimation <-function(X, modes, tau, lvs, lvs.names, IDM, 
                                             sets,scheme,method,fix.quantile, tol,
                                             maxiter)
{
  X = as.matrix(X)
  
  n<-nrow(X)
  
  # standardized Mvs
  X=scale(X)*sqrt((n)/(n-1))
  
  # number of manifest variables
  mvs = ncol(X)
  
  dimnames(IDM) = list(lvs.names, lvs.names)
  
  # create an object "blocks" which indicates the number
  # of variables in each block.
  blocks = unlist(lapply(sets, length))
  
  
  # MODE A E B
  # Outer Design Matrix.
  ODM = matrix(0,mvs,lvs)
  ODM_pvalue = matrix(0,mvs,lvs)
  ODM_pvaluebeta0 = matrix(0,mvs,lvs)
  ODMbeta0 = matrix(0,mvs,lvs)
  ODMbeta0_ModeB=matrix(0,1,lvs)
  sumODMbeta0=matrix(0,1,lvs)
  colnames(ODMbeta0_ModeB)=lvs.names
  colnames(ODM)=lvs.names
  colnames(ODM_pvalue)=lvs.names
  colnames(ODM_pvaluebeta0)=lvs.names
  Intercept=matrix(0,1,lvs)
  colnames(Intercept)=lvs.names
  
  aux = 1
  for (k in 1:lvs){
    ODM[aux:sum(blocks[1:k]),k] = rep(1,blocks[k]) # values of 1
    #CRI non mettiamo 1 alle intercette iniziali, le lasciamo a zero
    #ODMbeta0[aux:sum(blocks[1:k]),k] = rep(1,blocks[k]) # CRI
    ODM_pvaluebeta0[aux:sum(blocks[1:k]),k] = rep(100,blocks[k])
    ODM_pvalue[aux:sum(blocks[1:k]),k] = rep(100,blocks[k])
    aux = sum(blocks[1:k]) + 1
  }
  ODM
  # initialization of the weight matrix W
  Wbeta1 = ODM %*% diag(1/(apply(X %*% ODM,2,stats::sd)),lvs,lvs)
  wbeta1.old = rowSums(Wbeta1) 
  w.dif = 1 # scalar for checking convergence
  # inizializzo la matrice delle intercette outer
  # per Mode A (serve una intercetta per ciascuna MV)
  Wbeta0 = ODMbeta0 %*% diag(1/(apply(X %*% ODM,2,stats::sd)),lvs,lvs)#CRI
  wbeta0.old = rowSums(Wbeta0)
  #w.old=c(wbeta1.old,wbeta0.old)
  w.old=c(wbeta1.old)
  #senza intercetta
  #w.old=c(wbeta1.old)
  
  IDM2=IDM
  itermax = 0 # index of number of iterations
  W.DIF=NULL
  ITER=NULL
  
  while (w.dif > tol && !is.na(w.dif) &&  itermax < maxiter){
    # SOLO MODE QA E QB
    #Y = X %*% W  # external estimation of latent variables 'Y' 
    # FINE SOLO MODE QA E QB
    
    # external estimation of latent variables 'Y'
    Y=X %*% Wbeta1
    
    for (k in 1:ncol(IDM)){ 
      follow <- IDM[, k] == 1
      if (sum(follow) > 0) {
        qcov<-matrix(0,nrow=sum(follow), ncol=1)
        Y_temp<-as.matrix(Y[, follow])
        for (i in 1:sum(follow)){     
          qcov[i]<-qc(x = Y[, k], y = Y_temp[, i], tau = tau)$rho
        }
        IDM2[follow, k] <-qcov
      }
    }
    
    E=IDM2+t(IDM2)
    if ((scheme="Factorial")==FALSE) {
      E=sign(E)
    }
    
    #print(E)
    Z = Y %*% E   # internal estimation of latent variables 'Z'
    
    #MODE A E B
    # computin outer weights
    aux = 0
    
    
    if(fix.quantile==TRUE){tau_outer=0.5}
    else {tau_outer = tau}
    
    for (k in 1:lvs){
      if (modes[k]=="A"){ # reflective way
        for (m in (aux+1):sum(blocks[1:k])){
          ODMbeta0[m,k]=suppressWarnings(quantreg::rq(X[,m]  ~  Z[, k], tau=tau_outer,method=method)$coefficients[1])
          ODM[m,k]=suppressWarnings(quantreg::rq(X[,m]  ~  Z[, k], tau=tau_outer,method=method)$coefficients[2])
        }
      }
      else if (modes[k]=="B"){ # formative way
        X.blok = matrix(X[,(aux+1):sum(blocks[1:k])], nrow=nrow(X), ncol=length(sets[[k]]))
        ODMbeta0[(aux+1):sum(blocks[1:k]),k] = suppressWarnings(quantreg::rq(Z[, k] ~X.blok, tau=tau_outer,method=method)$coefficients[1])
        ODMbeta0_ModeB[1,k] = suppressWarnings(quantreg::rq(Z[, k] ~X.blok, tau=tau_outer,method=method)$coefficients[1])
        ODM[(aux+1):sum(blocks[1:k]),k] = suppressWarnings(quantreg::rq(Z[, k] ~X.blok, tau=tau_outer,method=method)$coefficients[2:(ncol(X.blok)+1)])
      }
      aux = sum(blocks[1:k]) 
    }
    
    LV=X %*% ODM  
    
    aux = 0
    for (k in 1:lvs){
      X.blok = matrix(X[,(aux+1):sum(blocks[1:k])], nrow=nrow(X), ncol=length(sets[[k]]))
      Intercept[1,k] = suppressWarnings(quantreg::rq(Z[, k] ~X.blok, tau=tau_outer,method=method)$coefficients[1])
      aux = sum(blocks[1:k])
    }
    sweep(LV, 2, Intercept, "+")
    
    # standardizzo le LV e determino i nuovi pesi
    
    Wbeta1 = ODM %*% diag(1/(apply(LV,2,stats::sd)),lvs,lvs)
    #Wbeta1 = ODM
    colnames(Wbeta1)=colnames(IDM)
    rownames(Wbeta1)=colnames(X)
    wbeta1.new = rowSums(Wbeta1)
    
    w.new=wbeta1.new
    
    # compute outer weights difference
    w.dif = sum((w.old - w.new)^2)
    w.old = w.new
    itermax = itermax + 1
    ITER=rbind(ITER,itermax)
    W.DIF=rbind(W.DIF,w.dif) 
  } 
  
  
  list(outer.weights=wbeta1.new, Z=Z, 
       Wbeta1=Wbeta1, tau=tau,W.DIF=W.DIF,ITER=ITER, IDM=IDM, Intercept=Intercept)
  
}

