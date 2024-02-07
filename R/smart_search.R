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
                         alpha_power_modeling = .05,
                         lb = nrow(lavModel),
                         switchStep = round(steps/2),
                         constrainRelChange = TRUE,
                         CORES, verbose = TRUE,
                         FSmethod = "SL",
                         matchPI =TRUE,
                         PIcentering = "doubleMC",
                         liberalInspection = FALSE,
                         seeds,
                         test = "onesided",
                         ...)
{
     dotdotdot <- list(...)
     sim_seeds <- seeds
     Reps <- get_Reps(type = type, R = R, steps = steps)
     if(switchStep == Inf){
          Power_interval <- rep(0, steps)
     }else{
          Power_interval <- c(rep(0,switchStep), seq(.1, .01,
                                                     length.out = steps - switchStep))
     }
     if(constrainRelChange){
          Rel_tol <- c(seq(.5, 1, length.out = steps-1), Inf) # last step cannot be constrained!
     }else{
          Rel_tol <- rep(Inf, steps)
     }
     # let the precision increase in prediciting power relative to precision due increase replications
     if(alpha_power_modeling >= 1){
          Alpha_power_modeling <- rep(1, steps)
     }else{
          Alpha_power_modeling <- seq(1, alpha_power_modeling, length.out = steps)
     }
     Conditions <- data.frame(Reps, Power_interval, Rel_tol, Alpha_power_modeling)

     df_POI <- data.frame("matchLabel" = POI)
     fit_temp <- lavModel_Analysis[lavModel_Analysis$matchLabel %in% df_POI$matchLabel,,drop = FALSE]
     # sort by POI
     fit_temp <- merge(df_POI, fit_temp, by.x = "matchLabel", sort = FALSE)
     truth <- fit_temp$ustart; names(truth) <- POI

     Nfinal <- c(); Nnew <- N_start; df <- c()
     df_est <- c(); df_se <- c(); df_pvalue_onesided <- c(); df_pvalue_twosided <- c()
     df_sigs_onesided <- c(); df_sigs_twosided <- c(); vec_fitOK <- c()
     Nl <- ceiling(max(N_start / 2, lb)); Nu <- ceiling(N_start/2+N_start); Nall <- c()
     if(verbose) cat(paste0("Initiating smart search to find simulation based N for power of ",
                            power_aim, " within ", steps, " steps\nand in total ",
                            R, " replications. Ns are drawn randomly...\n"))
     for(i in 1:nrow(Conditions))
     {

          # sample sample sizes
          Ns <- sample(seq(Nl, Nu, 1), size = Conditions$Reps[i], replace = TRUE,
                       prob = dnorm(x = seq(-2,2,4/(-1+length(seq(Nl, Nu, 1))))))
          Ns <- sort(Ns, decreasing = TRUE); Nall <- c(Nall, Ns)
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
                                                                                        matchPI = matchPI,
                                                                                        PIcentering = PIcentering,
                                                                                        liberalInspection = liberalInspection,
                                                                                        sim_seed = sim_seeds[ni]),
                                      simplify = FALSE)
          if(CORES > 1L) parallel::stopCluster(cl)

          sim_seeds <- sim_seeds[-c(1:length(Ns))] # delete used seeds

          if(length(POI) == 1L)
          {
               temp_est <- sapply(Fitted, "[[", "est")
               temp_est <- as.matrix(temp_est)
               colnames(temp_est) <- POI; rownames(temp_est) <- NULL
               temp_se <- sapply(Fitted, "[[", "se")
               temp_se <- as.matrix(temp_se)
               colnames(temp_se) <- POI; rownames(temp_se) <- NULL
               }else{
               df_est <- rbind(df_est, t(sapply(Fitted, "[[", "est")))
               df_se <- rbind(df_se, t(sapply(Fitted, "[[", "se")))
          }
          vec_fitOK <- c(vec_fitOK, sapply(Fitted, "[[", "fitOK"))

          # compute p-values
          if(tolower(test) == "onesided")
          {
               trueMatrix <- matrix(rep(truth, each = nrow(df_est)),
                                    ncol = length(POI))
               pvalue <- pnorm(sign(trueMatrix)*as.matrix(df_est / df_se),
                               lower.tail = FALSE)
          }else if(tolower(test) == "twosided")
          {
               pvalue <- 2*pnorm(as.matrix(abs(df_est) / df_se),
                                 lower.tail = FALSE)
          }
          df <- data.frame(pvalue < alpha); df$Ns <- Nall
          df <- df[vec_fitOK, ]; names(df) <- c(POI, "Ns")
          ind_min <- which.min(colMeans(df, na.rm = TRUE))# find POI of lowest power
          relFreq_indMin <- mean(unlist(tail(df[, ind_min], n = length(Ns))),
                                 na.rm = TRUE)

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

     df_est <- data.frame(df_est); names(df_est) <- POI
     df_se <- data.frame(df_se); names(df_se) <- POI

     out <- list("N" = Nnew,
                 "N_trials" = Nfinal,
                 "est" = df_est, "se" = df_se,
                 "Ns" = Nall, "fitOK" = vec_fitOK,
                 "truth" = truth)

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
find_n_from_glm <- function(fit, pow = .8,
                            alpha_power_modeling, Nmax = 10^4,
                            power_modeling_method)
{
     N_sequence <- 1:Nmax
     if(class(fit)[1] == "WaldGLM")
     {
          if(all(is.na(fit$est))) return(NA)
          power <- Wald_pred_confint(out = fit, N_interest = N_sequence,
                                     alpha = alpha_power_modeling)$P_lb
     }else{
          # did the model converge
          if(!fit$converged) return(NA)

          # predict linear model
          reg_fit <- predict(object = fit, newdata = data.frame("Ns" = N_sequence),
                             se.fit = TRUE)
          # transform linear model to power (depending on the modeling method)
          if(tolower(power_modeling_method) == "logit")
          {
               power_trans_lb <- reg_fit$fit - qnorm(p = 1-alpha_power_modeling/2)*reg_fit$se.fit
               power <- exp(power_trans_lb)/(1+exp(power_trans_lb))
          }else if(tolower(power_modeling_method) == "probit")
          {
               power_trans_lb <- reg_fit$fit - qnorm(p = 1-alpha_power_modeling/2)*reg_fit$se.fit
               power <- pnorm(power_trans_lb)
          }
     }

     minN <- suppressWarnings(min(N_sequence[power >= pow]))
     if(abs(minN) == Inf & Nmax == 10^6) return(Inf)
     if(abs(minN) == Inf) minN <- find_n_from_glm(fit = fit, pow = pow,
                                                  alpha_power_modeling =  alpha_power_modeling,
                                                  power_modeling_method = power_modeling_method,
                                                  Nmax = 10^6)
     return(minN)
}

# fit power model
fit_power_model <- function(Nnew, Nl, Nu, lb, ind_min,
                            power_modeling_method,
                            df, relFreq_indMin,
                            power_aim, i = 0, switchStep = Inf,
                            Conditions, alpha_power_modeling) {
     if(length(table(df$Ns)) > 1)
     {
          if(tolower(power_modeling_method) == "wald")
          {
               fit <- suppressWarnings(fitWaldglm(sig = na.omit(df)[,ind_min],
                                                  Ns = na.omit(df)$Ns))
               # if Wald-GLM did not converge, retry with probit (one time)
               if(all(is.na(fit$est))) power_modeling_method <- "probit"
          }
          if(tolower(power_modeling_method) == "logit")
          {
               fit <- glm(df[,ind_min] ~ I(sqrt(Ns)), family = binomial(link = "logit"), data = df)
          }else if(tolower(power_modeling_method) == "probit")
          {
               fit <- glm(df[,ind_min] ~ I(sqrt(Ns)), family = binomial(link = "probit"), data = df)
          }else{
               stop("This power modeling method has not been implemented.")
          }
          Nnew_temp <- find_n_from_glm(fit = fit, pow = power_aim,
                                       alpha_power_modeling =  Conditions$Alpha_power_modeling[i],
                                       power_modeling_method = power_modeling_method)
          if(i <= switchStep)
          {
               Nl_temp <- find_n_from_glm(fit = fit, pow = .15, alpha_power_modeling = 1,
                                          power_modeling_method = power_modeling_method)
               Nu_temp <- find_n_from_glm(fit = fit, pow = .85, alpha_power_modeling = Conditions$Alpha_power_modeling[i],
                                          power_modeling_method = power_modeling_method)
          }else{
               Nl_temp <- find_n_from_glm(fit = fit, pow = max(power_aim - Conditions$Power_interval[i], .0001),
                                          alpha_power_modeling = 1,
                                          power_modeling_method = power_modeling_method)
               Nu_temp <- find_n_from_glm(fit = fit, pow = min(power_aim + Conditions$Power_interval[i], .9999),
                                          alpha_power_modeling = Conditions$Alpha_power_modeling[i],
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
          Nnew <- evaluate_N(N_temp = Nnew_temp, N = Nnew,
                             relFreq_indMin = relFreq_indMin,
                             fit = fit, lb = lb,
                             rel_tol = Conditions$Rel_tol[i])
          Nl <- evaluate_N(N_temp = Nl_temp, N = Nl,
                           relFreq_indMin = relFreq_indMin,
                           fit = fit, lb = lb,
                           rel_tol = Conditions$Rel_tol[i])
          Nu <- evaluate_N(N_temp = Nu_temp, N = Nu,
                           relFreq_indMin = relFreq_indMin,
                           fit = fit, lb = lb,
                           rel_tol = Conditions$Rel_tol[i])
     }else{
          Nnew <- Nnew_temp; Nl <- Nl_temp; Nu <- Nu_temp
     }

     if(Nu <= Nl){
          Nl <- round(Nu/2)
     }

     N_out <- list("Nnew" = Nnew, "Nl" = Nl, "Nu" = Nu)
     return(N_out)
}

# evaluate resonablity of N
evaluate_N <- function(N_temp, N, relFreq_indMin, fit, lb, rel_tol = .5)
{
     # handle possible non-convergence
     if(is.na(N_temp))  N_temp <- ifelse(relFreq_indMin < .5, Inf, -Inf)

     if(class(fit)[1] == "WaldGLM")
     {
          if(any(fit$est < 0)) N_temp <- ifelse(relFreq_indMin < .5, Inf, -Inf)
     }else{
          if(any(coef(fit)[-1] < 0)) N_temp <- ifelse(relFreq_indMin < .5, Inf, -Inf)
     }


     if(abs(N_temp - N) > rel_tol*N){
          N <- round(N + sign(N_temp - N)*rel_tol*N)
     }else{
          N <- N_temp
     }
     # if variation is too low, change sample size drastically
     if(relFreq_indMin > .99)
     {
          N <- round(N*rel_tol+1)
     }else if(relFreq_indMin < .01)
     {
          N <- round(N/rel_tol)
     }
     return(max(lb, N))
}




