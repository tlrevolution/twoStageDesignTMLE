#' summary.twoStageTMLE
#'
#' Summarizes estimation procedure for missing 2nd stage covariates
#'
#' @param object An object of class \code{twoStageTMLE}
#' @param ... Other arguments passed to the tmle function in the tmle package
#' 
#' @return A list containing the missingness model, terms, coefficients, type,
#	and whether discrete or ensemble Super Learning was sused
#'
#' @export
summary.twoStage <- function(object,...){
	if (!is.null(object$coef)) {
		picoef <- object$coef
        if (inherits(picoef, "matrix")) {
            piterms <- colnames(picoef)
         }else {
            piterms <- names(picoef)
          }
         pimodel <- paste("Delta.W ~ 1")
         if (length(piterms) > 1) {
             pimodel <- paste("Delta.W ~ ", paste(piterms, collapse = " + "))
         }
    } else {
    	pimodel <- piterms <- picoef <- NULL
    }
    return(list(pimodel=pimodel, piterms=piterms, picoef=picoef, pitype=object$type,
						pidiscreteSL=object$discreteSL))   
}
	
	
#' print.summary.twoStageTMLE
#'
#' @param x an object of class summary.twoStageTMLE
#' @param ... additional arguments (i)
#' @importFrom tmle print.tmle
#'
#' @return print object
#' @export
#' 
#' 
#' @method print summary.twoStageTMLE
#'
print.summary.twoStageTMLE <- function(x,...){
	if (inherits(x, "summary.twoStageTMLE")){
		cat("\n Estimation of Pi (subset sampling mechanism)\n")
        cat("\t Procedure:", x$twoStage$pitype)
        if (!(is.null(x$twoStage$pidiscreteSL))) {
            if (x$twoStage$pidiscreteSL) {
                cat(", discrete")
            }
            else {
                cat(", ensemble")
            }
        }
        if (!(is.null(x$twoStage$piAUC))) {
            cat("\t Empirical AUC =", round(x$twoStage$piAUC, 4), "\n")
        }
        cat("\n")
        if (!(is.na(x$twoStage$picoef[1]))) {
            cat("\t Model:\n\t\t", x$twoStage$pimodel, "\n")
            cat("\n\t Coefficients: \n")
            terms <- sprintf("%15s", x$twoStage$piterms)
            extra <- ifelse(x$twoStage$picoef >= 0, "  ", " ")
            for (i in 1:length(x$twoStage$picoef)) {
                cat("\t", terms[i], extra[i], x$twoStage$picoef[i], "\n")
            }
        }
        cat("\n")
		print(x$tmle)
	}
}

#' print.twoStageTMLE
#'
#' @param x an object of class twoStageTMLE
#' @param ... additional arguments (i)
#'
#' @return print tmle results using print.tmle
#' method from tmle package
#' @importFrom tmle print.tmle 
#' @export
#' 
#'
#' @method print twoStageTMLE
#'
print.twoStageTMLE <- function(x,...){
	cat("Subset calibration TMLE\n")
	if (inherits(x, "twoStageTMLE")){
		print(x$tmle)
	}
}

#' summary.twoStageTMLE
#'
#' @param object an object of class twoStageTMLE
#' @param ... additional arguments (ignored)
#'
#' @return list summarizing the two-stage procedure components, 
#'  summary of the twoStage missingness estimation
#'  summary of the tmle for estimating the parameter
#' @export
#' 
#' 
#' @method summary twoStageTMLE
#'
summary.twoStageTMLE <- function(object,...){
	# summary for estimating Pi here
	if (inherits(object, "twoStageTMLE")){
		sum.twoStageTMLE <- list()
		sum.twoStageTMLE$twoStage <- summary.twoStage(object$twoStage)
		sum.twoStageTMLE$tmle <- summary(object$tmle)
		class(sum.twoStageTMLE) <- "summary.twoStageTMLE"
		return(sum.twoStageTMLE)
	}
}
