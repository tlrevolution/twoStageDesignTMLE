#' Utilities
#' setV
#'  Set the number of cross-validation 
#'  folds as a function of effective sample size
#'  See Phillips 2023 doi.org/10.1093/ije/dyad023
#' @param n.effective the effective sample size
#'
#' @return the number of cross-validation folds
#' @export
#'
setV <- function(n.effective){
  
  if(n.effective <= 30){
    V <- n.effective 
  } else if (n.effective <= 500) {
    V <- 20
  } else if (n.effective <= 1000) {
    V <- 10
  } else if (n.effective <= 10000){
    V <- 5
  } else {
    V <- 2
  }
  return(V)
}

#' twoStageDesignTMLENews
#' Get news about recent updates and bug fixes
#' @param ... ignored 
#' @export
twoStageDesignTMLENews <- function(...){
  utils::RShowDoc("NEWS", package="twoStageTMLE",...)
}


#' .getColNames
#'  get all the column names for variables
#'  to condition on when evaluating Stage 1
#'  missingness probabilities
#'  valid entries are Y, A, (both vectors)
#'  where W is a matrix
#' @param condSetNames which sets of variables to include in outcome regression model
#' @param Wnames names of columns in baseline covariates \code{W}
#' @param Vnames names of columns in MSM strata-definition \code{V}
#' @noRd
#' @return column names for variables to condition on when
#'  modeling missingness probabilites
#'
.getColNames <- function(condSetNames, Wnames, Vnames = NULL){
  orig.colnames <- NULL
  for (i in 1:length(condSetNames)){
    if (condSetNames[i] %in% c("A", "Y")) {
      orig.colnames <- c(orig.colnames, condSetNames[i]) 
    } else if (condSetNames[i] == "W"){
      orig.colnames <- c(orig.colnames, Wnames)
    } else if (condSetNames[i] == "V"){
      orig.colnames <- c(orig.colnames, Vnames)
    }
  }
  return(orig.colnames)
}

#' .evalAugW
#'  calls TMLE to use super learner to evalute preliminary predictions for 
#'  Q(0,W) and Q(1,W) conditioning on stage 1 covariates
#' @param Y outcome vector
#' @param A binary treatment indicator
#' @param W covariate matrix
#' @param Delta outcome missingness indicator
#' @param id identifier of i.i.d. unit
#' @param family outcome regression family
#' @param SL.library super learner library for outcome regression modeling 
#' @return \code{W.Q}, nx2 matrix of outcome predictions based on stage 1 
#' covariates
#' @export
evalAugW <- function(Y, A, W, Delta, id, family, SL.library){
  n <- length(Y)
  n.effective <- sum(Delta)
  if(family == "binomial"){
    n.effective <- min(c(table(Y[Delta == 1]) * 5, n.effective))
  }
  V.Q <- setV(n.effective)
  W.Q <- tmle::tmle(Y, A, W, Delta = Delta, id = id, family = family, 
              g1W = rep(0.5, n), pDelta1 = cbind(A0 = rep(1,n), A1 = rep(1,n)),
              V.Q = V.Q, 
              Q.SL.library = SL.library)$Qinit$Q
  colnames(W.Q) <- c("W.Q0", "W.Q1")
  return(W.Q)
}

#' .bound
#' truncate input values to lie between upper and lower bounds
#' @noRd
#' @param x, any numeric scalar, vector, matrix, data.frame
#' @param bounds vector whose min and max values are the bounds
#' @returns \code{x} with outliers truncated
.bound <- function (x, bounds) 
{
  x[x > max(bounds)] <- max(bounds)
  x[x < min(bounds)] <- min(bounds)
  return(x)
}
