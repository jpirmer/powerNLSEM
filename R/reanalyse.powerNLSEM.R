#' Reanalyse powerNLSEM object
#' @param out object of class powerNLSEM
#' @param power_modeling_method Character indicating the power modeling method used. Default to \code{NULL}, indicating to use the same power modeling method as was used in the \code{powerNLSEM} object.
#' @param powerLevels Power levels for which the desired sample sizes should be computed. Needs to be a vector. Default to \code{NULL} indicating to use same power rate used in \code{powerNLSEM} object.
#' @param alpha Alpha value used for confidence intervals. Default to \code{NULL} indicating to use same alpha used in \code{powerNLSEM} object.
#' @returns Returns list of desired sample sizes per effect for each \code{powerLevel}. \code{Nall} refers to the sample size required per power level for all coefficients. \code{Npower} is a matrix containing the desired sample sizes per effect for every power level.
#' @import stats
#' @import utils
#' @export

reanalyse.powerNLSEM <- function(out, powerLevels = NULL, power_modeling_method = NULL, alpha = NULL, uncertainty_method = NULL)
{
     if(class(out)[1] != "powerNLSEM") stop("powerNLSEM object required.")
     if(is.null(powerLevels)) powerLevels <- out$power
     if(is.null(power_modeling_method)) power_modeling_method <- out$power_modeling_method
     if(is.null(alpha)) alpha <- out$alpha
     if(is.null(uncertainty_method)) uncertainty_method <- ""
     search_method <- out$search_method
     method <- out$method

     if(!is.vector(powerLevels)) stop("powerLevels needs to be a vector.")

     Sigs <- out$SigDecisions
     temp_list <- lapply(1:(ncol(Sigs)-1),
                         FUN = function(i) {fit <- glm(Sigs[,i]~I(sqrt(Ns)), data = Sigs,
                                                                        family = binomial(link = power_modeling_method))
                                            if((nrow(Sigs)-sum(Sigs[,i])<5) | fit$deviance < 10^-8) fit <- NULL # probably not converged
                                            return(fit)})
     Npower <- sapply(temp_list, function(FIT){sapply(powerLevels,
                                                      FUN = function(powerlevel){
                                               if(is.null(FIT)) return(NA)
                                               N_alpha <- find_n_from_glm(fit = FIT, pow = powerlevel, alpha = alpha,
                                               uncertainty_method = uncertainty_method, Nmax = 10^6,
                                               power_modeling_method = "probit")
               if(N_alpha == Inf) return(-Inf)
               return(N_alpha)})})
     if(length(powerLevels) == 1) Npower <- t(Npower)
     colnames(Npower) <- colnames(Sigs[,-ncol(Sigs)])
     rownames(Npower) <- powerLevels

     Nall <- apply(Npower, 1, function(x) max(x, na.rm = TRUE))

     # return ----
     out <- list("Nall" = Nall, "Npower" = Npower)
     out$power <- powerLevels
     out$beta <- 1-out$power
     out$alpha <- alpha
     out$method <- method
     out$search_method <- search_method
     out$power_modeling_method <- power_modeling_method

     class(out) <- c("powerNLSEM.reanalyzed", "list")
     return(out)
}

