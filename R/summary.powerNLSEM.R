#' Summary function for powerNLSEM objects
#' @import stats
#' @param object Result of powerNLSEM function estimating the MSPE. \code{object} must be of class \code{"powerNLSEM"}.
#' @param test Should the parameter be tested with a directed hypothesis (onesided) or with an undirected hypothesis (twosided, also equivalent to Wald-Test for single parameter). Default to \code{NULL} (if \code{NULL}, \code{test} of the original MSPE is used).
#' @param alpha Type I-error rate for significance decision. Default to \code{NULL} (if \code{NULL}, \code{alpha} of the original MSPE is used).
#' @param ... Further arguments to use in \code{summary}.
#' @return summary of powerNLSEM object
#' @exportS3Method summary powerNLSEM
#' @export

summary.powerNLSEM <- function(object, test = NULL, alpha = NULL, ...){
     if(is.null(test)) test <- object$args$test
     if(is.null(alpha)) alpha <- object$args$alpha

     if(tolower(test) == "onesided")
     {
          trueMatrix <- matrix(rep(object$truth, each = nrow(object$est)),
                               ncol = ncol(object$est))
          pvalue <- pnorm(sign(trueMatrix)*as.matrix(object$est / object$se),
                          lower.tail = FALSE)
     }else if(tolower(test) == "twosided")
     {
          pvalue <- 2*pnorm(as.matrix(abs(object$est) / object$se),
                            lower.tail = FALSE)
     }
     significanceDecision <- data.frame(pvalue < alpha); significanceDecision$Ns <- object$Ns
     significanceDecision <- significanceDecision[object$fitOK, ]
     names(significanceDecision) <- c(names(object$est), "Ns")
     object$significanceDecision <- significanceDecision
     object$powersPerN <- aggregate(.~Ns, data = object$significanceDecision, FUN = mean)

     class(object) <- "summary.powerNLSEM"
     return(object)
}

