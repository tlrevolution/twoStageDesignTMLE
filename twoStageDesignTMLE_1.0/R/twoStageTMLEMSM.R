#' twoStageTMLEmsm
#'
#' Inverse probability of censoring weighted TMLE for evaluating MSM
#' parameters when the full set of covariates is available on only
#' a subset of observations, as in a 2-stage design.
#'
#' @param Y  outcome of interest (missingness allowed)
#' @param A  binary treatment indicator 
#' @param W  matrix or data.frame of covariates measured on entire population
#' @param V vector, matrix, or dataframe of covariates used to define MSM strata
#' @param Delta.W Indicator of inclusion in subset with additional information
#' @param W.stage2 matrix or data.frame of covariates measured in subset population
#			  (in same order as the corresponding rows in \code{W}. DO NOT include
#			   rows for subjects not included in the subset)
#' @param Delta binary indicator that outcome Y is observed
#' @param	pi optional vector of sampling probabilities
#' @param piform parametric regression formula for estimating \code{pi} (see Details)
#' @param pi.SL.library super learner library for estimating \code{pi} (see Details)
#' @param V.pi optional number of cross-validation folds for super learning 
#'    (ignored when piform or pi is provided)
#' @param	pi.discreteSL flag to indicate whether to use ensemble or discrete 
#'    super learning (ignored when piform or pi is provided)
#' @param condSetNames Variables to include as predictors of missingness
#'  in \code{W.stage2}, any combination of \code{Y, A}, and either \code{W}  
#'  (for all covariates in \code{W}), or individual covariate names in \code{W}
#' @param id optional indicator of independent units of observation
#' @param	Q.family outcome regression family, "gaussian" or "binomial"
#' @param augmentW set to \code{TRUE} to augment \code{W} with predicted outcome 
#'    values when \code{A = 0} and \code{A = 1}
#' @param augW.SL.library super learner library for preliminary outcome
#' regression model (ignored when \code{augmentW} is \code{FALSE})
#' @param rareOutcome when \code{TRUE} sets \code{V.Q = 20, Q.discreteSL = TRUE}, 
#'  \code{Q.SL.library} includes glm, glmnet, bart
#' @param	verbose when \code{TRUE} prints informative messages 
#' @param	...	other arguments passed to the \code{tmleMSM} function 
#' 
#' @return 
#' \describe{
#' Object of class "twoStageTMLE"
#'   \item{tmle}{Treatment effect estimates and summary information from 
#'   call to \code{tmleMSM} function}
#'   \item{twoStage}{IPCW weight estimation summary, 
#'   \code{pi} are the probabilities,\code{coef} are SL weights or coefficients 
#'   from glm fit, \code{type} of estimation procedure, 
#'   \code{discreteSL} flag indicating whether discrete super learning was used}
#'   \item{augW}{Matrix of predicted outcomes based on stage 1 covariates only}
#' }
#' 
#' @details
#' When using \code{piform} to specify a parametric model for pi that conditions 
#' on the outcome use \code{Delta.W} as the dependent variable and \code{Y.orig}
#' on the right hand side of the formula instead of \code{Y}. When writing a 
#' user-defined SL wrapper for inclusion in \code{pi.SL.library} use \code{Y}
#' on the left hand side of the formula. If specific covariate names are
#' used on the right hand side use \code{Y.orig} to condition
#' on the outcome. 
#' 
#' @examples
#' n <- 1000
#' set.seed(10)
#' W1 <- rnorm(n)
#' W2 <- rnorm(n)
#' W3 <- rnorm(n)
#' A <- rbinom(n, 1, plogis(-1 + .2*W1 + .3*W2 + .1*W3))
#' Y <- 10 + A + W1 + W2 + A*W1 + W3 + rnorm(n)
#' Y.bin <- rbinom(n, 1, plogis(-4.6 - 1.8* A + W1 + W2 -.3 *A*W1 + W3))
#' # Set 400 obs with data on W3, more likely if W1 > 1
#' n.sample <- 400
#' p.sample <- 0.5 + .2*(W1 > 1)
#' rows.sample <- sample(1:n, size = n.sample, p = p.sample)
#' Delta.W <- rep(0,n)
#' Delta.W[rows.sample] <- 1
#' W3.stage2 <- cbind(W3 = W3[Delta.W==1])
#'
#' # 1. specify parametric models, misspecified outcome model (not recommended)
#' result1.MSM <- twoStageTMLEmsm(Y=Y, A=A, V= cbind(W1), W=cbind(W2), 
#' Delta.W = Delta.W, W.stage2 = W3.stage2, augmentW = FALSE,
#' piform = "Delta.W~ I(W1 > 0)", MSM = "A*W1", augW.SL.library = "SL.glm",
#' Qform = "Y~A+W1",gform="A~W1 + W2 +W3", hAVform = "A~1", verbose=TRUE)

#'summary(result1.MSM)
#'
#' # 2. Call again, passing in previously estimated observation weights, 
#' # note that specifying a correct model for Q improves efficiency
#' result2.MSM <- twoStageTMLEmsm(Y=Y, A=A, V= cbind(W1), W=cbind(W2), 
#' Delta.W = Delta.W, W.stage2 = W3.stage2, augmentW = FALSE,
#' pi = result1.MSM$twoStage$pi, MSM = "A*W1",
#' Qform = "Y~ A + W1 + W2 + A*W1 + W3",gform="A~W1 + W2 +W3", hAVform = "A~1")
#' cbind(SE.Qmis = result1.MSM$tmle$se, SE.Qcor = result2.MSM$tmle$se)
#' 
#' \donttest{
#' #Binary outcome, augmentW, rareOutcome
#' result3.MSM <- twoStageTMLEmsm(Y=Y.bin, A=A, V= cbind(W1), W=cbind(W2), 
#' Delta.W = Delta.W, W.stage2 = W3.stage2, augmentW = TRUE,
#' piform = "Delta.W~ I(W1 > 0)", MSM = "A*W1", gform="A~W1 + W2 +W3",
#'  Q.family = "binomial", rareOutcome=TRUE)
#'}
#' 
#' @seealso
#' * [tmle::tmleMSM()] for details on customizing the estimation procedure
#' * [twoStageTMLE()] for estimating marginal effects
#' @export
twoStageTMLEmsm <- function(Y, A, W, V, Delta.W, W.stage2, 
      Delta = rep(1, length(Y)),pi=NULL, piform=NULL,
      pi.SL.library = c("SL.glm", "SL.gam", "SL.glmnet", "tmle.SL.dbarts.k.5"),
      V.pi=10, pi.discreteSL = TRUE, condSetNames = c("A","V", "W", "Y"), 
      id = NULL, Q.family="gaussian", augmentW = TRUE, 
      augW.SL.library = c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2"),
      rareOutcome = FALSE, 
      verbose=FALSE, ...) {
    
    if(is.null(id)){id <- 1:length(Y)}
    
    if (is.vector(W)){
  		W <- as.matrix(W)
  		colnames(W) <- "W1"
  	}
  
 	 if (is.null(colnames(W))){
  		colnames(W) <- paste0("W", 1:ncol(W))
  	}
  
    
    if(is.vector(W.stage2)){
      W.stage2 <- as.matrix(W.stage2)
      colnames(W.stage2) <- "W.stage2"
    }
    
    if (augmentW){	
      W.Q <- evalAugW(Y, A, W, Delta, id, Q.family, augW.SL.library)			
    } else {
      W.Q <- NULL
    }
    
    V <- as.matrix(V)
    
    # # Evaluate conditional sampling probabilities 
    res.twoStage <- estimatePi(Y=Y, A=A, W=W, condSetNames=condSetNames, W.Q=W.Q, 
                               Delta.W = Delta.W, piform = piform, id = id,
                               pi.SL.library = pi.SL.library,  V = V.pi, 
                               discreteSL = pi.discreteSL, 
                               verbose=verbose, pi = pi)
    # if (is.null(pi)){		
    #   validCondSetNames <- c("A", "V", "W","Y")
    #   if (!(all(condSetNames %in% validCondSetNames))) {
    #     stop("condSetNames must be any combination of 'A', 'V', 'W', 'Y'")
    #   }		
    #   if (any(condSetNames == "Y")){
    #     if (any(is.na(Y))){
    #       stop("Cannot condition on the outcome to evaluate sampling probabilities when some outcome values are missing")
    #     }
    #   }
    #   
    #   if(is.null(W.Q)){
    #     d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames))
    #     colnames(d.pi)[-1] <- .getColNames(condSetNames, 
    #                              c(colnames(W)), colnames(V))	
    #   } else {
    #     d.pi <- data.frame(Delta.W = Delta.W, mget(condSetNames), W.Q)
    #     colnames(d.pi)[-1] <- .getColNames(condSetNames, 
    #                              c(colnames(W), colnames(W.Q)), colnames(V))	
    #   }
    # 
    #   res.twoStage <- tmle::estimateG(d = d.pi, g1W=NULL, gform = piform,
    #         SL.library = pi.SL.library, id=id, V=V.pi,
    #         message = "sampling weights", outcome = "A", 
    #         discreteSL = pi.discreteSL, obsWeights = rep(1, nrow(W)),
    #         verbose=verbose)
    #   names(res.twoStage)[1] <- "pi"		
    # } else {
    #   res.twoStage <- list()
    #   res.twoStage$pi <- pi
    #   res.twoStage$type <- "User supplied values"
    #   res.twoStage$coef <- NA
    #   res.twoStage$discreteSL <- NULL
    # }
    # 
    ub <- sqrt(sum(Delta.W)) * log(sum(Delta.W)) / 5
    obsWeights <- .bound(Delta.W/res.twoStage$pi, c(0, ub))
    
    # Now set call to tmle on full data
    argList <- list(...)
    argList$Y <- Y[Delta.W==1]
    argList$A <- A[Delta.W==1] 
	if (is.null(W.Q)){
	  	  argList$W <- cbind(W[Delta.W==1,, drop=FALSE], W.stage2)
	  } else {
	  	 argList$W <- cbind(W[Delta.W==1,, drop=FALSE], W.Q[Delta.W==1,], W.stage2)
	  }
    argList$V = V[Delta.W==1, , drop = FALSE]
    argList$Delta <- Delta[Delta.W==1]
    argList$obsWeights <- obsWeights[Delta.W==1]
    argList$id <- id[Delta.W==1]
    argList$family <- Q.family
    argList$verbose <- verbose
    if (rareOutcome){
      argList$Q.SL.library <- c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2") 
      argList$Q.discreteSL <- TRUE
      argList$V.Q <- 20
    }		
    
    result <- list()	
    result$tmleMSM <- try(do.call(tmle::tmleMSM, argList))			 			
    
    if(inherits(result$tmleMSM, "try-error")){
      warning("Error calling tmleMSM. Estimated sampling probabilites are being returned")
      result$tmle <- NULL
    }
    result$twoStage <- res.twoStage
    result$augW <- W.Q
    class(result) <- "twoStageTMLE"  # for tmleMSM can use same class
    return(result)					
}