#' plot powerNLSEM object
#' @param x object of class powerNLSEM
#' @param test Should the parameter be tested with a directed hypothesis (onesided) or with an undirected hypothesis (twosided, also equivalent to Wald-Test for single parameter). Default to \code{NULL}, then the same as in fitted \code{powerNLSEM} object in \code{x} is used.
#' @param plot Character indicating what type of plot to create. Default to \code{"power_model"}, referencing to the prediction of significant parameters using the model specified in \code{power_modeling_method}.
#' @param power_modeling_method Character indicating the power modeling method used. This is only relevant when \code{plot = "power_model"} is used. Default to \code{NULL}, indicating to use the same power modeling method as was used in the \code{powerNLSEM} function.
#' @param se Logical indicating to use confidence intervals based on normal approximation using the standard errors. Default to \code{FALSE}.
#' @param power_aim Power level to be included into the plot with respective N. If \code{NULL} the same power level as in the \code{powerNLSEM} function will be used. If set to \code{0} no power level and corresponding N will be plotted. Default to \code{NULL}, indicating to use the same power modeling method as was used in the \code{powerNLSEM} function.
#' @param alpha Alpha value used for confidence intervals, when \code{se = TRUE}. Default to \code{NULL}, indicating to use the same alpha as was used in the powerNLSEM function. This does not influence the significance decision, although same alpha is used per default.
#' @param alpha_power_modeling Type I-error rate for confidence band around predicted power rate. Used to ensure that the computed \code{N} keeps the desired power value (with the given Type I-error rate \code{alpha_power_modeling} divided by 2). If set to 1, no confidence band is used. Default to \code{.05}.
#' @param min_num_bins minimal number of bins used for aggregating results. Default to 10.
#' @param defaultgg Logical to return default ggplot object. Default to \code{FALSE}, which returns \code{theme_minimal} and other changes in theme.
#' @param ... Additional arguments passed on to the plot function.
#' @returns Returns \code{ggplot} object of the type specified in plot.
#' @import ggplot2
#' @import stats
#' @import utils
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
#' plot(Result_Power)
#' }
#' @export

plot.powerNLSEM <- function(x, test  = NULL,
                            plot = "power_model",
                            power_modeling_method = NULL, se = FALSE,
                            power_aim = NULL, alpha = NULL,
                            alpha_power_modeling = NULL,
                            min_num_bins = 10,
                            defaultgg = FALSE,
                            ...)
{
     out <- x

     # to get around R checks Ns, Power, Effect, Power_ub, Power_lb are names within a data.frame which are called in ggplot and hence are not defined gloabally:
     Ns <- Power <- Effect <- Power_ub <- Power_lb <- NULL
     if(is.null(power_modeling_method))
     {
          power_modeling_method <- out$power_modeling_method
     }
     if(is.null(alpha))
     {
          alpha <- out$alpha
     }
     if(is.null(test)) test <- out$test
     if(is.null(alpha_power_modeling)) alpha_power_modeling <- out$alpha_power_modeling

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

     if(tolower(plot) == "power_model")
     {
          if(power_modeling_method == "probit")
          {
               temp  <- sapply(1:(ncol(Sigs)-1), function(i) predict(newdata = data.frame("Ns" = c(min(Sigs$Ns, na.rm = TRUE):max(Sigs$Ns, na.rm = TRUE))),
                                                                     glm(Sigs[,i]~I(sqrt(Ns)), data = Sigs,
                                                                         family = binomial(link = "probit")),
                                                                     se.fit = TRUE), simplify = FALSE)
               Probit <- c(); Probit_UB <- c(); Probit_LB <- c()
               for(i in 1:length(temp))
               {
                    Probit <- cbind(Probit, temp[[i]]$fit)
                    Probit_UB <- cbind(Probit_UB, temp[[i]]$fit + qnorm(1-alpha_power_modeling/2)*temp[[i]]$se.fit)
                    Probit_LB <- cbind(Probit_LB, temp[[i]]$fit - qnorm(1-alpha_power_modeling/2)*temp[[i]]$se.fit)
               }
               powers <- pnorm(Probit)
               powers_UB <- pnorm(Probit_UB)
               powers_LB <- pnorm(Probit_LB)
          }else if(tolower(power_modeling_method) == "wald")
          {
               temp  <- lapply(1:(ncol(Sigs)-1), function(i) { fitWald <- suppressWarnings(fitWaldglm(sig = na.omit(Sigs)[, i],
                                                                                     Ns = na.omit(Sigs)$Ns))
               if(all(is.na(fitWald$est))){
                    NonCov <- data.frame(matrix(NA, ncol = 3,
                                     nrow =  length(c(min(Sigs$Ns, na.rm = TRUE):max(Sigs$Ns, na.rm = TRUE)))))
                    names(NonCov) <- c("P", "P_lb", "P_ub")
                    return(NonCov)
               }
                    WaldCI <- suppressWarnings(Wald_pred_confint(out = fitWald,
                                                                 N_interest = c(min(Sigs$Ns, na.rm = TRUE):max(Sigs$Ns, na.rm = TRUE)),
                                                                 alpha = alpha_power_modeling))
                    return(WaldCI)})
               powers <- c(); powers_UB <- c(); powers_LB <- c()
               for(i in 1:length(temp))
               {
                    powers <- cbind(powers, temp[[i]]$P)
                    powers_UB <- cbind(powers_UB, temp[[i]]$P_ub)
                    powers_LB <- cbind(powers_LB, temp[[i]]$P_lb)
               }
               # replace unplausible values
               powers[powers > 1] <- 1
               powers_UB[powers_UB > 1] <- 1
               powers_LB[powers_LB > 1] <- 1
               powers[powers < 0] <- 0
               powers_UB[powers_UB < 0] <- 0
               powers_LB[powers_LB < 0] <- 0
          }else if(power_modeling_method == "logit")
          {
               temp  <- sapply(1:(ncol(Sigs)-1), function(i) predict(newdata = data.frame("Ns" = c(min(Sigs$Ns, na.rm = TRUE):max(Sigs$Ns, na.rm = TRUE))),
                                                                     glm(Sigs[,i]~I(sqrt(Ns)), data = Sigs,
                                                                         family = binomial(link = "logit")),
                                                                     se.fit = TRUE), simplify = FALSE)
               Logit <- c(); Logit_UB <- c(); Logit_LB <- c()
               for(i in 1:length(temp))
               {
                    Logit <- cbind(Logit, temp[[i]]$fit)
                    Logit_UB <- cbind(Logit_UB, temp[[i]]$fit + qnorm(1-alpha_power_modeling/2)*temp[[i]]$se.fit)
                    Logit_LB <- cbind(Logit_LB, temp[[i]]$fit - qnorm(1-alpha_power_modeling/2)*temp[[i]]$se.fit)
               }
               powers <- exp(Logit)/(1 + exp(Logit))
               powers_UB <- exp(Logit_UB)/(1 + exp(Logit_UB))
               powers_LB <- exp(Logit_LB)/(1 + exp(Logit_LB))
          }

          if(se){
               nonconvergence_UB <- apply(powers_UB, 2, function(x) all(x==1, na.rm = TRUE))
               nonconvergence_LB <- apply(powers_LB, 2, function(x) all(x==0, na.rm = TRUE))
               nonconvergence_both <- colMeans(powers_UB - powers_LB, na.rm = TRUE) > .95 # large CIs

               nonconvergence_UB[is.na(nonconvergence_UB)] <- TRUE
               nonconvergence_LB[is.na(nonconvergence_LB)] <- TRUE
               nonconvergence_both[is.na(nonconvergence_both)] <- TRUE

               if(any(nonconvergence_UB)) powers_UB[, nonconvergence_UB] <- powers[, nonconvergence_UB]
               if(any(nonconvergence_LB)) powers_LB[, nonconvergence_LB] <- powers[, nonconvergence_LB]
               if(any(nonconvergence_both)) powers_LB[, nonconvergence_both] <- powers_UB[, nonconvergence_both] <- powers[, nonconvergence_both]

               df_pred <- cbind(powers, powers_UB, powers_LB, c(min(Sigs$Ns, na.rm = TRUE):max(Sigs$Ns, na.rm = TRUE))) |> data.frame()
               names(df_pred) <- c(names(Sigs)[names(Sigs)!="Ns"],
                                   paste0("ub_", names(Sigs)[names(Sigs)!="Ns"]),
                                   paste0("lb_", names(Sigs)[names(Sigs)!="Ns"]),
                                   "Ns")
               df_pred <- df_pred[order(df_pred$Ns),]
               df_long <- reshape(df_pred,
                                  varying = list(names(df_pred)[!(grepl(pattern = "ub_", names(df_pred)) |
                                                                       grepl(pattern = "lb_", names(df_pred))) &
                                                                     (names(df_pred) != "Ns")],
                                                 names(df_pred)[grepl(pattern = "ub_", names(df_pred))],
                                                 names(df_pred)[grepl(pattern = "lb_", names(df_pred))]),
                                  direction = "long", v.names = c("Power", "Power_ub", "Power_lb"),
                                  times = names(df_pred)[!(grepl(pattern = "ub_", names(df_pred)) |
                                                                grepl(pattern = "lb_", names(df_pred)))
                                                         & names(df_pred) != "Ns"],
                                  timevar = "Effect")
               gg <- ggplot(df_long, aes(x = Ns, y = Power, col = Effect, fill = Effect))+
                    ylab("Predicted Power")+xlab("N")+
                    geom_ribbon(aes(x=Ns, y = Power, ymax = Power_ub, ymin = Power_lb), alpha = .3, lwd = 0.1)+
                    geom_line(lwd=1)+
                    ggtitle(paste0("Model implied power with confidence bands for ", out$method),
                            subtitle = paste0("using ", power_modeling_method, " regression for ", test, " test"))
          }else{
               df_pred <- cbind(powers, c(min(Sigs$Ns, na.rm = TRUE):max(Sigs$Ns, na.rm = TRUE))) |> data.frame(); names(df_pred) <- names(Sigs)
               df_pred <- df_pred[order(df_pred$Ns),]
               df_long <- reshape(df_pred, varying = list(names(df_pred)[names(df_pred) != "Ns"]),
                                  direction = "long", v.names = c("Power"),
                                  times = names(df_pred)[names(df_pred) != "Ns"], timevar = "Effect")
               gg <- ggplot(data = df_long, aes(Ns, Power, col = Effect))+
                     ylab("Predicted Power")+xlab("N")+
                     geom_line(lwd = 1)+ggtitle(paste0("Model implied power for ", out$method),
                                               subtitle = paste0("using ", power_modeling_method, " regression for ", test, " test"))
          }

     }else if(tolower(plot) == "empirical")
     {
          SUMMARY <- aggregate(.~ Ns, data = Sigs, FUN = function(x) mean(x = x, na.rm = TRUE))
          SUMMARY$Ns_count <- as.numeric(table(Sigs$Ns))
          NAMES <- names(SUMMARY)
          SUMMARY <- data.frame(SUMMARY)
          names(SUMMARY) <- NAMES

          grouping <- rep(0, length(SUMMARY$Ns_count)); i <- 1; ind_min <- 1
          while(i < length(SUMMARY$Ns_count)+1)
          {
               if(suppressWarnings(min(which(cumsum(SUMMARY$Ns_count[grouping == 0]) >= min_num_bins))) == Inf){
                    grouping[grouping == 0] <- i - 1
                    break
               }
               ind_max <- min(which(cumsum(SUMMARY$Ns_count[grouping == 0]) >= min_num_bins)) + ind_min - 1
               grouping[ind_min:ind_max] <- i
               ind_min <- ind_max + 1; i <- i + 1
          }
          SUMMARY_agg <- c()
          for(i in unique(grouping))
          {
               SUMMARY_agg <- rbind(SUMMARY_agg,
                                    colSums(SUMMARY[grouping == i,, drop = FALSE] *
                                                 SUMMARY$Ns_count[grouping == i])/
                                         sum(SUMMARY$Ns_count[grouping == i]))
          }; SUMMARY_agg <- data.frame(SUMMARY_agg); names(SUMMARY_agg) <- NAMES
          SUMMARY_agg$Ns_count <- NULL



          df_long <- reshape(SUMMARY_agg, varying = list(names(SUMMARY_agg)[names(SUMMARY_agg) != "Ns"]),
                             direction = "long", v.names = c("Power"),
                             times = names(SUMMARY_agg)[names(SUMMARY_agg) != "Ns"], timevar = "Effect")
          gg <- ggplot(data = df_long, aes(Ns, Power, col = Effect, fill = Effect))+geom_point(cex = .1)+
               geom_smooth(method = "loess", formula = "y~x", level = 1-alpha_power_modeling)+
               ylab("Predicted Power")+xlab("N")+
               ggtitle(paste0("Model implied power with confidence bands for ", out$method),
                       subtitle = paste0("using LOESS  for ", test, " test"))
     }

     if(is.null(power_aim) &
        (test == out$test) &
        (power_modeling_method == out$power_modeling_method) &
        (alpha == out$alpha) &
        (alpha_power_modeling == out$alpha_power_modeling)){
          gg <- gg + geom_hline(yintercept = out$power, lwd = .5, lty = 3)+
               geom_vline(xintercept = out$N, lwd = .5, lty = 3)
     }else if(all(power_aim < 1) & all(power_aim > 0))
     {
          if(is.null(power_aim)) power_aim <- out$power
          temp <- reanalyse.powerNLSEM(out, powerLevels = power_aim, test = test,
                               power_modeling_method = power_modeling_method,
                               alpha = alpha, alpha_power_modeling = alpha_power_modeling)
          gg <- gg + geom_hline(yintercept = temp$power, lwd = .5, lty = 3)+
               geom_vline(xintercept = temp$Npower, lwd = .5, lty = 3)

     }
     gg <- gg + theme_minimal(base_size = 16)
     return(gg)
}
