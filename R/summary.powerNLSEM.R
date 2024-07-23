#' Summary function for powerNLSEM objects
#' @import stats
#' @param object Result of powerNLSEM function estimating the MSPE. \code{object} must be of class \code{"powerNLSEM"}.
#' @param test Should the parameter be tested with a directed hypothesis (onesided) or with an undirected hypothesis (twosided, also equivalent to Wald-Test for single parameter). Default to \code{NULL} (if \code{NULL}, \code{test} of the original MSPE is used).
#' @param alpha Type I-error rate for significance decision. Default to \code{NULL} (if \code{NULL}, \code{alpha} of the original MSPE is used).
#' @param ... Further arguments to use in \code{summary}.
#' @return summary of powerNLSEM object
#' @examples
#' \donttest{
#' # write model in lavaan syntax
#' model <- "
#' # measurement models
#'           X =~ 1*x1 + 0.8*x2 + 0.7*x3
#'           Y =~ 1*y1 + 0.85*y2 + 0.78*y3
#'           Z =~ 1*z1 + 0.9*z2 + 0.6*z3
#'
#' # structural models
#'           Y ~ 0.3*X + .2*Z +  .2*X:Z
#'
#' # residual variances
#'          Y~~.7975*Y
#'          X~~1*X
#'          Z~~1*Z
#'
#' # covariances
#'          X~~0.5*Z
#'
#' # measurement error variances
#'          x1~~.1*x1
#'          x2~~.2*x2
#'          x3~~.3*x3
#'          z1~~.2*z1
#'          z2~~.3*z2
#'          z3~~.4*z3
#'          y1~~.5*y1
#'          y2~~.4*y2
#'          y3~~.3*y3
#' "
#' # run model-implied simulation-based power estimation
#' # for the effects: c("Y~X", "Y~Z", "Y~X:Z")
#' Result_Power <- powerNLSEM(model = model, POI = c("Y~X", "Y~Z", "Y~X:Z"),
#'                            method = "UPI", search_method = "adaptive",
#'                            steps = 10, power_modeling_method = "probit",
#'                            R = 1000, power_aim = .8, alpha = .05,
#'                            alpha_power_modeling = .05,
#'                            CORES = 1, seed = 2024)
#'
#' Result_Power
#' summary(Result_Power)
#' }
#' @export
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

