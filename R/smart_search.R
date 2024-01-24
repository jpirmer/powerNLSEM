#' @import stats
#' @import parallel
#' @importFrom pbapply pbsapply

smart_search <- function(POI,
                         method, lavModel,
                         lavModel_Analysis,
                         data_transformations,
                         search_method,
                         power_modeling_method,
                         N_start = nrow(lavModel)*10, type = "u",
                         R = 1000, steps = 10,
                         power_aim = .8, alpha = .05,
                         lb = nrow(lavModel),
                         switchStep = round(steps/2),
                         CORES, verbose = TRUE,
                         uncertainty_method = "",
                         FSmethod = "SL",
                         seeds,
                         test = "onesided",
                         ...)
{
     dotdotdot <- list(...)
     sim_seeds <- seeds
     Reps <- get_Reps(type = type, R = R, steps = steps)
     Power_interval <- c(rep(0,switchStep), seq(.1, .01, -(.1-.01)/(steps-switchStep-1)))
     Rel_tol <- seq(.5, 1, (1-.5)/(steps-1))
     Conditions <- data.frame(Reps, Power_interval, Rel_tol)
     Nfinal <- c(); Nnew <- N_start; df <- c()
     df_est <- c(); df_se <- c(); df_pvalue_onesided <- c(); df_pvalue_twosided <- c()
     df_sigs_onesided <- c(); df_sigs_twosided <- c(); vec_fitOK <- c()
     Nl <- max(N_start / 2, lb); Nu <- N_start/2+N_start
     if(verbose) cat(paste0("Initiating smart search to find simulation based N for power of ",
                            power_aim, " within ", steps, " steps\nand in total ",
                            R, " replications. Ns are drawn randomly...\n"))
     for(i in 1:nrow(Conditions))
     {

          # sample sample sizes
          Ns <- sample(seq(Nl, Nu, 1), size = Conditions$Reps[i], replace = TRUE,
                       prob = dnorm(x = seq(-2,2,4/(-1+length(seq(Nl, Nu, 1))))))
          Ns <- sort(Ns, decreasing = TRUE)
          if(verbose) cat(paste0("Step ", i, " of ", steps, ". Fitting ", length(Ns),
                                 " models with Ns in [", Nl, ", ", Nu,"].\n"))

          # simulate data,  fit models and evaluate significance of parameters of interest
          lavModel_attributes <- lavaan::lav_partable_attributes(lavModel)
          matrices <- get_matrices(lavModel, lavModel_attributes)

          if(CORES > 1L){
               cl <- parallel::makeCluster(CORES)
          parallel::clusterExport(cl = cl, varlist = ls(), envir = environment())
          parallel::clusterEvalQ(cl = cl, expr = {
               library(powerNLSEM)})
          }else{cl <- NULL}
          Fitted <- pbapply::pbsapply(cl = cl,
                                      X = seq_along(Ns), FUN = function(ni) sim_and_fit(n = Ns[ni], POI = POI, alpha = alpha,
                                                                                        lavModel = lavModel,
                                                                                        method = method,
                                                                                        lavModel_Analysis = lavModel_Analysis,
                                                                                        lavModel_attributes = lavModel_attributes,
                                                                                        matrices = matrices,
                                                                                        data_transformations = data_transformations,
                                                                                        prefix = ni,
                                                                                        FSmethod = FSmethod,
                                                                                        sim_seed = sim_seeds[ni]),
                                      simplify = FALSE)
          if(CORES > 1L) parallel::stopCluster(cl)

          sim_seeds <- sim_seeds[-c(1:length(Ns))] # delete used seeds



          ###### HERE a lot of redundant stuff is being done!!!
          ###### Rethink the use of Sigs and df in the original function and redo this part
          ###### ----------------------------------------------------------------------------------------------------------->


          df_est <- rbind(df_est, t(sapply(Fitted, "[[", "est")))
          df_se <- rbind(df_se, t(sapply(Fitted, "[[", "se")))
          df_pvalue_onesided <- rbind(df_pvalue_onesided, t(sapply(Fitted, "[[", "pvalue_onesided")))
          df_pvalue_twosided <- rbind(df_pvalue_twosided, t(sapply(Fitted, "[[", "pvalue_twosided")))
          df_sigs_onesided <- rbind(df_sigs_onesided, t(sapply(Fitted, "[[", "sigs_onesided")))
          df_sigs_twosided <- rbind(df_sigs_twosided, t(sapply(Fitted, "[[", "sigs_twosided")))
          temp_fitOK <- sapply(Fitted, "[[", "fitOK")
          vec_fitOK <- rbind(vec_fitOK, temp_fitOK)

          if(tolower(test) == "onesided")
          {
               Sigs <- data.frame(df_sigs_onesided, Ns); names(Sigs) <- c(colnames(df_sigs_onesided), "Ns")
          }else if(tolower(test) == "twosided")
          {
               Sigs <- data.frame(df_sigs_twosided, Ns); names(Sigs) <- c(colnames(df_sigs_twosided), "Ns")
          }
          Sigs$fitOK <- temp_fitOK
          df <- rbind(df, Sigs[Sigs$fitOK, ]); rm(Ns)

          ind_min <- which.min(colMeans(df, na.rm = TRUE))# find POI of lowest power

          ### run power model
          args <- names(formals(fit_power_model))
          args <- args[args!="..."]
          N_temp <- try(do.call("fit_power_model", mget(args)), silent = TRUE)
          if(inherits(N_temp, "try-error"))
          {
               N_temp <- list("Nnew" = Nnew, "Nl" = max(round(Nl/2), lb), "Nu" = Nu)
          }
          Nnew <- N_temp$Nnew; Nl <- N_temp$Nl; Nu <- N_temp$Nu


          Nfinal <- c(Nfinal, Nnew)
     }

     if(tolower(test) == "onesided")
     {
          SigDecisions <- data.frame(df_sigs_onesided, Ns)
          names(SigDecisions) <- c(colnames(df_sigs_onesided), "Ns")
     }else if(tolower(test) == "twosided")
     {
          SigDecisions <- data.frame(df_sigs_twosided, Ns)
          names(SigDecisions) <- c(colnames(df_sigs_twosided), "Ns")
     }
     out <- list("N" = Nnew, SigDecisions = SigDecisions,
                 "N_trials" = Nfinal,
                 df_est = df_est, df_se = df_se)

     return(out)
}

# generate sequence of Ns
get_Reps <- function(type = "u", R = 1000, steps = 10) {
     if(tolower(type) == "increasing")
     {
          temp <- R/(steps*(steps+1)/2)
          Ns <- seq(temp, temp*steps, temp) |> round()
          temp_diff <- R - sum(Ns)
          Ns[which.max(Ns)] <- Ns[which.max(Ns)] + temp_diff
     }else if(tolower(type) == "u"){
          temp <- (R/2)/((steps + steps%%2)/2*((steps + steps%%2)/2+1)/2)
          Ns <- seq(temp, temp*(steps + steps%%2)/2, temp) |> round()
          temp_diff <- R/2 - sum(Ns)
          Ns[which.max(Ns)] <- Ns[which.max(Ns)] + temp_diff
          Ns <- c(sort(Ns, decreasing = TRUE), Ns)
          if(length(Ns) != steps)
          {
               Ns <- c(Ns[1:((steps+1)/2-1)],
                       sum(Ns[((steps+1)/2):(((steps+1)/2+1))]),
                       Ns[1+((steps+1)/2+1):(steps)])
          }
     }else if(tolower(type) == "equal"){
          temp <- R/steps
          Ns <- rep(temp, steps) |> round()
          temp_diff <- R - sum(Ns)
          Ns[which.max(Ns)] <- Ns[which.max(Ns)] + temp_diff
     }
     return(Ns)
}

# find n from an glm-fit model
find_n_from_glm <- function(fit, pow = .8, alpha = .05,
                            uncertainty_method = "", Nmax = 10^4,
                            power_modeling_method)
{
     N_sequence <- 1:Nmax
     if(uncertainty_method == "exact")
     {
          alpha <- 1 # no influence on se
     }
     # predict linear model
     reg_fit <- predict(object = fit, newdata = data.frame("Ns" = N_sequence), se.fit = TRUE)
     # transform linear model to power (depending on the modeling method)
     if(tolower(power_modeling_method) == "logit")
     {
          power_trans_lb <- reg_fit$fit - qnorm(p = 1-alpha/2)*reg_fit$se.fit
          power <- exp(power_trans_lb)/(1+exp(power_trans_lb))
     }else if(tolower(power_modeling_method) == "probit")
     {
          power_trans_lb <- reg_fit$fit - qnorm(p = 1-alpha/2)*reg_fit$se.fit
          power <- pnorm(power_trans_lb)
     }
     minN <- suppressWarnings(min(N_sequence[power >= pow]))
     if(abs(minN) == Inf & Nmax == 10^6) return(Inf)
     if(abs(minN) == Inf) minN <- find_n_from_glm(fit = fit, pow = pow, alpha = alpha,
                                          uncertainty_method =  uncertainty_method, Nmax = 10^6)
     return(minN)
}

# fit power model
fit_power_model <- function(Nnew, Nl, Nu, Sigs, lb,
                            power_modeling_method, df, ind_min,
                            power_aim, alpha, i = 0, switchStep = Inf,
                            Conditions, uncertainty_method = "") {
     if(length(table(df$Ns)) > 1)
     {
          if(tolower(power_modeling_method) == "logit")
          {
               fit <- glm(df[,ind_min] ~ I(sqrt(Ns)), family = binomial(link = "logit"), data = df)
          }else if(tolower(power_modeling_method) == "probit")
          {
               fit <- glm(df[,ind_min] ~ I(sqrt(Ns)), family = binomial(link = "probit"), data = df)
          }else{
               stop("This power modeling method has not been implemented.")
          }
          Nnew_temp <- find_n_from_glm(fit = fit, pow = power_aim, alpha = alpha,
                                       uncertainty_method =  uncertainty_method,
                                       power_modeling_method = power_modeling_method)
          if(i <= switchStep)
          {
               Nl_temp <- find_n_from_glm(fit = fit, pow = .15, alpha = 1, uncertainty_method =  uncertainty_method,
                                          power_modeling_method = power_modeling_method)
               Nu_temp <- find_n_from_glm(fit = fit, pow = .85, alpha = alpha, uncertainty_method =  uncertainty_method,
                                          power_modeling_method = power_modeling_method)
          }else{
               Nl_temp <- find_n_from_glm(fit = fit, pow = max(power_aim - Conditions$Power_interval[i], .0001),
                                          alpha = 1, uncertainty_method =  uncertainty_method,
                                          power_modeling_method = power_modeling_method)
               Nu_temp <- find_n_from_glm(fit = fit, pow = min(power_aim + Conditions$Power_interval[i], .9999),
                                          alpha = alpha, uncertainty_method = uncertainty_method,
                                          power_modeling_method = power_modeling_method)
          }
     }else{
          fit <- NULL
          Nnew_temp <- Nl_temp <- Nu_temp <-  unique(df$Ns)
     }

     if(i == (switchStep+1)) # reset process after switch
     {
          Nl <- Nnew - 1; Nu <- Nnew + 1
     }

     # new iteration of Ns
     if(!is.null(fit))
     {
          Nnew <- evaluate_N(N_temp = Nnew_temp, N = Nnew, Sigs = Sigs,
                             ind_min = ind_min, fit = fit, lb = lb,
                             rel_tol = Conditions$Rel_tol[i])
          Nl <- evaluate_N(N_temp = Nl_temp, N = Nl, Sigs = Sigs,
                           ind_min = ind_min, fit = fit, lb = lb,
                           rel_tol = Conditions$Rel_tol[i])
          Nu <- evaluate_N(N_temp = Nu_temp, N = Nu, Sigs = Sigs,
                           ind_min = ind_min, fit = fit, lb = lb,
                           rel_tol = Conditions$Rel_tol[i])
     }else{
          Nnew <- Nnew_temp; Nl <- Nl_temp; Nu <- Nu_temp
     }

     if(Nu <= Nl){
          Nl <- round(Nu/2)
     }

  N_out <- list("Nnew" = Nnew, "Nl" = Nl, "Nu_temp" = Nu)
  return(N_out)
}


# evaluate resonablity of N
evaluate_N <- function(N_temp, N, Sigs, ind_min, fit, lb, rel_tol = .5)
{
     if(any(coef(fit)[-1] < 0))
     {
          N_temp <- ifelse(mean(unlist(Sigs[,ind_min]), na.rm = TRUE) < .5, Inf, -Inf)
     }

     if(abs(N_temp - N) > rel_tol*N){
          N <- round(N + sign(N_temp - N)*rel_tol*N)
     }else{
          N <- N_temp
     }
     # if variation is too low, change sample size drastically
     if(mean(unlist(Sigs[,ind_min]), na.rm = TRUE) > .99)
     {
          N <- round(N*rel_tol+1)
     }else if(mean(unlist(Sigs[,ind_min]), na.rm = TRUE) < .01)
     {
          N <- round(N/rel_tol)
     }
     return(max(lb, N))
}
