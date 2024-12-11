#' estimatePi
#'
#' Typically not called directly by the user. Function for modeling
#' the two-stage missingness mechanism and evaluating conditional 
#' probabilities for each observation
#' 
#' @param Y outcome
#' @param A binary treatment indicator
#' @param W covariate matrix observed on everyone
#' @param condSetNames Variables to include as predictors of missingness
#'  in \code{W.stage2}, any combination of \code{Y, A}, and either \code{W} 
#'  (for all covariates in \code{W}) or individual covariate names in \code{W}
#' @param W.Q additional covariates based on preliminary outcome regression
#' @param Delta.W binary indicator of missing second stage covariates
#' @param piform parametric regression formula for estimating \code{pi} 
#' @param pi.SL.library super learner library for estimating \code{pi} 
#' @param V number of cross validation folds for estimating \code{pi} 
#'  using super learner
#' @param discreteSL Use discrete super learning when \code{TRUE}, otherwise
#'  ensemble super learning
#' @param id Identifier of independent units of observation, e.g., clusters
#' @param verbose When \code{TRUE} prints informational messages 
#' @param pi optional vector of user-specified probabilities
#' @param obsWeights optional weights for evaluating pi
#' @return list containing the predicted probabilities, estimation method
#' coefficients in parametric regression model (if piform supplied),
#' indicator of whether discrete or ensemble SL was used.
#' 
#' @importFrom stats as.formula terms
#' @export
#'
estimatePi <- function(Y, A, W, condSetNames, W.Q, Delta.W, piform,
                       pi.SL.library, id, V, discreteSL, verbose, pi=NULL,
                       obsWeights = rep(1, nrow(W))){
  if (is.null(pi)){		
    validCondSetNames <- c("A", "W","Y", colnames(W))
    if (!(all(condSetNames %in% validCondSetNames))) {
      stop("condSetNames must be any combination of 'A', 'Y', specific column
      names in W, or 'W' to include all columns.")
    }		
    if (any(condSetNames == "Y")){
      if (any(is.na(Y))){
        stop("Cannot condition on the outcome to evaluate sampling probabilities
             when some outcome values are missing")
      }
    }
    
    if(is.null(W.Q)){
      d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames))
      colnames(d.pi)[-1] <- temp <- .getColNames(condSetNames, 
                                                 c(colnames(W)), colnames(V))	
    } else {
      d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames), W.Q)
      colnames(d.pi)[-1] <- .getColNames(condSetNames, 
                                         c(colnames(W), colnames(W.Q)), colnames(V))	
    }
    # stop if  piform includes Y on the right hand side
    # (don't automatically replace because too difficult if interaction terms are in 
    # the model) 
    if(!is.null(piform)){
      varNames <- as.character(attr(terms(as.formula(piform)), "variables"))[-(1:2)]
      if ("Y" %in% varNames) {
        stop("Error in specified parametric formula for evaluating sampling
             probabilities. To condition on Y use 'Y.orig' on the right hand side
             of the regression formula")
      }
    }
    res.twoStage <- tmle::estimateG(d = d.pi, g1W=NULL, gform = piform,
                      SL.library = pi.SL.library, id=id, V=V, 
                      message = "sampling weights", outcome = "A",
                      discreteSL = discreteSL, obsWeights = obsWeights,
                      verbose=verbose)
    names(res.twoStage)[1] <- "pi"
  } else {
    res.twoStage <- list()
    res.twoStage$pi <- pi
    res.twoStage$type <- "User supplied values"
    res.twoStage$coef <- NA
    res.twoStage$discreteSL <- NULL
  }
  return(res.twoStage)
}

