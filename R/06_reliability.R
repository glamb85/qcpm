#' @title Measurement model reliability and internal consistence
#'
#' @description
#' \code{reliability} returns the classical indices used in PLS-PM to assess 
#' the reliability and internal consistence of the measurement model (Hair et al., 2019). 
#' In order it provides: Cronbach's alpha, Dillon-Goldstein's rho, the Dijkstra-Henseler rho, and 
#' first and second eigenvalue of the correlation matrix of the manifest variables. The function 
#' also returns the outer mode (A or B) and the number of manifest variables for each block.
#' 
#' @details
#' The function only returns Dijkstra-Henseler rho values for quantile 0.5. When mode 
#' B is selected, or there are some intra-block inverse correlations, the Dijkstra-Henseler rho, 
#' Cronbach's alpha, and Dillon-Goldstein's rho are not calculated.
#' 
#'
#' @param qcpm is an object of class \code{qcpm}
#' 
#' @return A table containing, for each block, the outer mode (A or B), 
#' the number of manifest variables, Cronbach's alpha, Dillon-Goldstein's rho, 
#' Dijkstra-Henseler rho, and first and second eigenvalue of the manifest variable 
#' correlation matrix.
#' 
#' 
#' @author Cristina Davino, Pasquale Dolce, Giuseppe Lamberti, Domenico Vistocco
#'
#' 
#' 
#' @references Hair, J.F., Risher, J.J., Sarstedt, M. and Ringle, C.M. (2019). When to 
#' use and how to report the results of PLS-SEM. \emph{European Business Review}, \bold{31 (1)},
#'  pp. 2--24, doi: 10.1108/EBR-11-2018-0203
#' 
#' @references Sanchez, G. (2013). PLS Path Modeling with R Trowchez Editions. Berkeley, 2013. 
#' Available at \url{https://www.gastonsanchez.com/PLS_Path_Modeling_with_R.pdf}.
#' 
#' @seealso \code{\link{qcpm}}, \code{\link{assessment}}, \code{\link{boot}}, and 
#' \code{\link{summary}}
#' 
#' @importFrom stats cor princomp sd var
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
#' reliability(well.qcpm)   
#'
reliability <- function(qcpm)
{
  if (class(qcpm) != "qcpm") 
    stop("Argument 'qcpm' must be an object of class 'qcpm'")
  
  DM = qcpm$data
  blocks = qcpm$model$qcpm.outer.list
  modes =  qcpm$model$qcpm.modes
  if(is.null(qcpm$model$tau) == FALSE && qcpm$model$tau==0.5){
    ow = qcpm$model$outer.weights2[,1]
  }
  else{
    ow = qcpm$model$outer.weights2[,colnames(qcpm$model$outer.weights2)=="0.5"]
  }
  # inputs setting
  lvs = length(blocks) 
  lvs_names = rownames(qcpm$model$qcpm.inner)
  
  blockinds = NULL
  for (i in 1: lvs) blockinds=c(blockinds,rep(i,length(blocks[[i]])))
  
  block_sizes = lengths(blocks)
  #  blocklist = unlist(lapply(block_sizes, function(x) rep(x, x)))
  obs = nrow(DM)
  sdvf = sqrt((nrow(DM)-1) / nrow(DM)) 
  
  # Unidimensionality
  Alpha = rep(1, lvs)    # Cronbach's Alpha for each block
  Rho = rep(1, lvs)      # D.G. Rho for each block
  RhoA = rep(1,lvs)
  eig.1st = rep(1, lvs)  # first eigenvalue
  eig.2nd = rep(0, lvs)  # second eigenvalue
  
  # calculate indices
  for (aux in 1:lvs)
  {
    if (block_sizes[aux] != 1) 
    { 
      # scaling data
      DM.block = DM[,blockinds==aux]
      weights.block = ow[blockinds==aux]
      stdev.X = apply(DM.block, 2, sd) * sdvf 
      X_uni = scale(DM.block, scale=stdev.X)
      if (nrow(X_uni) < ncol(X_uni)) {   # more columns than rows
        acp = princomp(t(X_uni)) 
        X.rho = t(X_uni)
      } else {   # more rows than columns
        acp = princomp(X_uni)
        X.rho = X_uni
      }
      if (modes[aux] == "A") 
      {
        p = ncol(X_uni)
        
        if(sum(sign(cor(X_uni))==-1)>0){
          Alpha[aux] = NA
          Rho[aux] = NA
          RhoA[aux] = NA  
          warning(paste("coefficients of block",lvs_names[aux], "are not calculated because of intra-block inverse correlations."))
        }
        else{
          # cronbach's alpha
          a.denom = var(rowSums(X_uni)) * sdvf^2
          a.numer = 2 * sum(cor(X_uni)[lower.tri(cor(X_uni))])
          alpha = (a.numer / a.denom) * (p / (p - 1))
          Alpha[aux] <- ifelse(alpha < 0, NA, alpha)
          # dillon-goldstein rho
          
          numer_rho <- colSums(cor(X.rho, acp$scores[,1]))^2
          denom_rho <- numer_rho + (p - colSums(cor(X.rho, acp$scores[,1])^2) )
          Rho[aux] <- numer_rho / denom_rho
          
          # rho_DH
          S = stats::cov(scale(DM.block))
          w = weights.block
          ww <- w %*% t(w)
          numer <- t(w) %*% (S-diag(diag(S))) %*% w
          denom <- t(w) %*% (ww-diag(diag(ww)) ) %*% w
          rhoa  <- (t(w) %*% w)^2 * numer / denom
          RhoA[aux] <- ifelse(rhoa < 0, NA, rhoa)
        }}
      else{
        Alpha[aux] = NA
        Rho[aux] = NA
        RhoA[aux] = NA
        warning(paste("coefficients of block",lvs_names[aux], "are not calculated because the construct is formative"))
        
      }
      eig.1st[aux] = acp$sdev[1]^2
      eig.2nd[aux] = acp$sdev[2]^2
    }
  }
  reli = data.frame(Mode = modes, 
                      MVs = round(block_sizes,4),
                      C.alpha = round(Alpha,4), 
                      DG.rho = round(Rho,4),
                      rho.A  = round(RhoA,4),
                      eig.1st = round(eig.1st,4), 
                      eig.2nd = round(eig.2nd,4))
  rownames(reli) = lvs_names
  return(reli)
}