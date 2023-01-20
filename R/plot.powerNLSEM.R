#' plot powerNLSEM object
#' @param out object of class powerNLSEM
#' @param min_num_bins minimal number of bins used for aggregating results. Default to 10.
#' @param plot Character indicating what type of plot to create. Default to "power_model", referencing to the prediction of significant parameters using the model specified in power_modeling_method.
#' @param power_modeling_method Character indicating the power modeling method used. This is only relevant when plot = "power_model" is used. Default to NULL, indicating to use the same power modeling method as was used in the powerNLSEM function.
#' @param se Logical indicating to use confidence intervals based on normal approximation using the standard errors. Default to FALSE.
#' @param alpha Alpha value used for confidence intervals, when se = TRUE. Default to NULL, indicating to use the same alpha as was used in the powerNLSEM function.
#' @returns Returns ggplot object of the type specified in plot.
#' @import ggplot2
#' @import dplyr
#' @export

plot.powerNLSEM <- function(out, min_num_bins = 10, plot = "power_model", power_modeling_method = NULL, se = F, alpha = NULL)
{
     if(is.null(power_modeling_method))
     {
          power_modeling_method <- out$power_modeling_method
     }
     if(is.null(alpha))
     {
          alpha <- out$alpha
     }
     if(tolower(plot) == "power_model")
     {
          Sigs <- out$SigDecisions
          if(power_modeling_method == "logit")
          {
               if(se)
               {
                    temp  <- sapply(1:(ncol(Sigs)-1), function(i) predict(newdata = data.frame("Ns" = c(min(Sigs$Ns, na.rm = T):max(Sigs$Ns, na.rm = T))),
                                                                          glm(Sigs[,i]~I(sqrt(Ns)), data = Sigs,
                                                                              family = binomial(link = "logit")),
                                                                          se.fit = T), simplify = F)
                    Logit <- c(); Logit_UB <- c(); Logit_LB <- c()
                    for(i in 1:length(temp))
                    {
                         Logit <- cbind(Logit, temp[[i]]$fit)
                         Logit_UB <- cbind(Logit_UB, temp[[i]]$fit + qnorm(1-alpha/2)*temp[[i]]$se.fit)
                         Logit_LB <- cbind(Logit_LB, temp[[i]]$fit - qnorm(1-alpha/2)*temp[[i]]$se.fit)
                    }
                    powers <- exp(Logit)/(1 + exp(Logit))
                    powers_UB <- exp(Logit_UB)/(1 + exp(Logit_UB))
                    powers_LB <- exp(Logit_LB)/(1 + exp(Logit_LB))
                    df_pred <- cbind(powers, powers_UB, powers_LB, c(min(Sigs$Ns, na.rm = T):max(Sigs$Ns, na.rm = T))) |> data.frame()
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
                         geom_hline(yintercept = out$power, lwd = .5, lty = 3)+ylab("Predicted Power")+xlab("N")+
                         geom_vline(xintercept = out$N, lwd = .5, lty = 3)+
                         geom_ribbon(aes(x=Ns, y = Power, ymax = Power_ub, ymin = Power_lb), alpha = .3, lwd = 0.1)+
                         geom_line(lwd=1)+theme_minimal(base_size = 16)+
                         ggtitle("Model implied power with confidence bands")

               }else{
                    Logit <- sapply(1:(ncol(Sigs)-1), function(i) predict(newdata = data.frame("Ns" = c(min(Sigs$Ns, na.rm = T):max(Sigs$Ns, na.rm = T))),
                                                                          glm(Sigs[,i]~I(sqrt(Ns)), data = Sigs,
                                                                              family = binomial(link = "logit"))))
                    powers <- exp(Logit)/(1 + exp(Logit))
                    df_pred <- cbind(powers, c(min(Sigs$Ns, na.rm = T):max(Sigs$Ns, na.rm = T))) |> data.frame(); names(df_pred) <- names(Sigs)
                    df_pred <- df_pred[order(df_pred$Ns),]
                    df_long <- reshape(df_pred, varying = list(names(df_pred)[names(df_pred) != "Ns"]),
                                       direction = "long", v.names = c("Power"),
                                       times = names(df_pred)[names(df_pred) != "Ns"], timevar = "Effect")
                    gg <- ggplot(data = df_long, aes(Ns, Power, col = Effect))+
                         geom_hline(yintercept = out$power, lwd = .5, lty = 3)+ylab("Predicted Power")+xlab("N")+
                         geom_vline(xintercept = out$N, lwd = .5, lty = 3)+theme_minimal(base_size = 16)+
                         geom_line(lwd = 1)+ggtitle("Model implied power")

               }

          }
     }else if(tolower(plot) == "empirical")
     {
          SUMMARY <- out$SigDecisions %>% group_by(Ns) %>% summarise_all(.funs = function(x) mean(x = x, na.rm = T)) %>% ungroup()
          SUMMARY$Ns_count <- as.numeric(table(out$SigDecisions$Ns))
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
                                    colSums(SUMMARY[grouping == i,, drop = F] *
                                                 SUMMARY$Ns_count[grouping == i])/
                                         sum(SUMMARY$Ns_count[grouping == i]))
          }; SUMMARY_agg <- data.frame(SUMMARY_agg); names(SUMMARY_agg) <- NAMES
          SUMMARY_agg$Ns_count <- NULL



          df_long <- reshape(SUMMARY_agg, varying = list(names(SUMMARY_agg)[names(SUMMARY_agg) != "Ns"]),
                             direction = "long", v.names = c("Power"),
                             times = names(SUMMARY_agg)[names(SUMMARY_agg) != "Ns"], timevar = "Effect")
          gg <- ggplot(data = df_long, aes(Ns, Power, col = Effect))+geom_point(cex = .1)+
               geom_smooth(method = "loess", formula = "y~x")+
               geom_hline(yintercept = out$power, lwd = .5, lty = 3)+
               geom_vline(xintercept = out$N, lwd = .5, lty = 3)+theme_minimal(base_size = 16)
     }

     return(gg)
}
