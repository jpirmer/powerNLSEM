#' Reanalyse powerNLSEM object
#' @param out object of class powerNLSEM
#' @param power_modeling_method Character indicating the power modeling method used. Default to \code{NULL}, indicating to use the same power modeling method as was used in the \code{powerNLSEM} object.
#' @param test Should the parameter be tested with a directed hypothesis (onesided) or with an undirected hypothesis (twosided, also equivalent to Wald-Test for single parameter). Default to \code{NULL}, then the same as in fitted \code{powerNLSEM} object is used.
#' @param powerLevels Power levels for which the desired sample sizes should be computed. Needs to be a vector. Default to \code{NULL} indicating to use same power rate used in \code{powerNLSEM} object.
#' @param alpha Type I-error rate for significance decision. Default to \code{.05}.
#' @param alpha_power_modeling Type I-error rate for confidence band around predicted power rate. Used to ensure that the computed \code{N} keeps the desired power value (with the given Type I-error rate \code{alpha_power_modeling} divided by 2). If set to 1, no confidence band is used. Default to \code{.05}.
#' @returns Returns list of desired sample sizes per effect for each \code{powerLevel}. \code{Nall} refers to the sample size required per power level for all coefficients. \code{Npower} is a matrix containing the desired sample sizes per effect for every power level.
#' @import stats
#' @import utils
#' @export

reanalyse.powerNLSEM <- function(out, test = NULL,
                                 powerLevels = NULL, power_modeling_method = NULL,
                                 alpha = NULL, alpha_power_modeling = NULL)
{
     if(class(out)[1] != "powerNLSEM") stop("powerNLSEM object required.")
     if(is.null(powerLevels)) powerLevels <- out$power
     if(is.null(power_modeling_method)) power_modeling_method <- out$power_modeling_method
     if(is.null(alpha)) alpha <- out$alpha
     if(is.null(alpha_power_modeling)) alpha_power_modeling <- out$alpha_power_modeling
     if(is.null(test)) test <- out$test
     search_method <- out$search_method
     method <- out$method

     if(!is.vector(powerLevels)) stop("powerLevels needs to be a vector.")

     # get significance decisions
     if(tolower(test) == "onesided")
     {
          truth <- out$truth
          trueMatrix <- matrix(rep(truth, each = nrow(out$est)),
                               ncol = ncol(out$est))
          pvalue <- pnorm(sign(trueMatrix)*as.matrix(out$est / out$se),
                          lower.tail = FALSE)
     }else if(tolower(test) == "twosided")
     {
          pvalue <- 2*pnorm(as.matrix(abs(out$est) / out$se),
                            lower.tail = FALSE)
     }
     Sigs <- data.frame(pvalue < alpha); names(Sigs) <- names(out$est)
     Sigs$Ns <- out$Ns
     Sigs <- Sigs[out$fitOK, ] # remove false convergences
     Sigs <- na.omit(Sigs)
     ind_min <- which.min(colMeans(Sigs))

     if(tolower(power_modeling_method) == "wald")
     {
          fit <- suppressWarnings(fitWaldglm(sig = Sigs[,ind_min],
                                             Ns = Sigs$Ns))
          # if Wald-GLM did not converge, retry with probit (one time)
          if(all(is.na(fit$est))) power_modeling_method <- "probit"
     }
     if(tolower(power_modeling_method) != "wald"){
          fit <- glm(Sigs[,ind_min] ~ I(sqrt(Ns)),
                     family = binomial(link = power_modeling_method),
                     data = Sigs)
     }
     if(tolower(power_modeling_method) != "wald")
     {
          Npower <- sapply(powerLevels,
                           FUN = function(powerlevel){
                                if(is.null(fit)) return(NA)
                                N_alpha <- find_n_from_glm(fit = fit, pow = powerlevel,
                                                           uncertainty_method = uncertainty_method, Nmax = 10^6,
                                                           alpha_power_modeling  = alpha_power_modeling)
                                if(N_alpha == Inf) return(-Inf)
                                return(N_alpha)})
     }else{
          fitWald <- fit
          N_temp <- 1:10^6
          df_temp_Wald <- Wald_pred_confint(fitWald, N_interest = N_temp,
                                            alpha = alpha_power_modeling)
          P_LB_temp <- df_temp_Wald$P_lb
          Npower <- sapply(X = powerLevels, function(p) min(N_temp[P_LB_temp >= p]))
     }



#
#      temp_list <- lapply(1:(ncol(Sigs)-1),
#                          FUN = function(i) {if(tolower(power_modeling_method) == "wald"){
#                               fit <- fitWaldglm(sig = Sigs[,i], Ns = Sigs$Ns)
#                               if(all(is.na(fit$est))) return(NULL)
#                               }else{
#                                    fit <- glm(Sigs[,i]~I(sqrt(Ns)), data = Sigs,
#                                          family = binomial(link = power_modeling_method))
#                               if((nrow(Sigs)-sum(Sigs[,i])<5) | fit$deviance < 10^-8) fit <- NULL # probably not converged
#                          }
#                               return(fit)})
#      Npower <- sapply(temp_list, function(FIT){sapply(powerLevels,
#                                                       FUN = function(powerlevel){
#                                                if(is.null(FIT)) return(NA)
#                                                N_alpha <- find_n_from_glm(fit = FIT, pow = powerlevel, alpha = alpha,
#                                                uncertainty_method = uncertainty_method, Nmax = 10^6,
#                                                power_modeling_method = power_modeling_method)
#                if(N_alpha == Inf) return(-Inf)
#                return(N_alpha)})})
#      if(length(powerLevels) == 1) Npower <- t(Npower)
#      colnames(Npower) <- colnames(Sigs[,-ncol(Sigs)])
#      rownames(Npower) <- powerLevels

     # return ----
     out <- list("Npower" = Npower)
     out$power <- powerLevels
     out$beta <- 1-out$power
     out$alpha <- alpha
     out$alpha_power_modeling <- alpha_power_modeling
     out$method <- method
     out$test <- test
     out$search_method <- search_method
     out$power_modeling_method <- power_modeling_method

     class(out) <- c("powerNLSEM.reanalyzed", "list")
     return(out)
}

