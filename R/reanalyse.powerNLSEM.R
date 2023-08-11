#' Reanalyse powerNLSEM object
#' @param out object of class powerNLSEM
#' @param power_modeling_method Character indicating the power modeling method used. This is only relevant when plot = "power_model" is used. Default to NULL, indicating to use the same power modeling method as was used in the powerNLSEM function.
#' @param powerLevels Power levels for which the desired sample sizes should be computed. Needs to be a vector.
#' @param alpha Alpha value used for confidence intervals, when se = TRUE. Default to NULL, indicating to use the same alpha as was used in the powerNLSEM function.
#' @returns Returns list of sample sizes per effect. \code{Nall} refers to the sample size required per power level for all coefficients. \code{Npower} is a matrix containing the desired sample sizes per effect for every power level.
#' @import stats
#' @import utils
#' @export

reanalyse.powerNLSEM <- function(out, power_modeling_method, powerLevels, alpha, uncertainty_method = "")
{
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

     out <- list("Nall" = Nall, "Npower" = Npower)
     return(out)
}

