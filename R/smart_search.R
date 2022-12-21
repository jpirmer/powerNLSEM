
smart_search <- function(POI,
                         method, lavModel,
                         lavModel_Analysis,
                         data_transformations,
                         search_method,
                         N_start = nrow(lavModel)*10, type = "u",
                         Ntotal = 1000, steps = 10,
                         power_aim = .8, alpha = .05,
                         lb = nrow(lavModel),
                         switchStep = round(steps/2),
                         CORES, verbose = T,
                         ...)
{
     dotdotdot <- list(...)

     Reps <- get_Reps(type = type, Ntotal = Ntotal, steps = steps)
     Power_interval <- c(rep(0,switchStep), seq(.1, .01, -(.1-.01)/(steps-switchStep-1)))
     Rel_tol <- seq(.5, 1, (1-.5)/(steps-1))
     Conditions <- data.frame(Reps, Power_interval, Rel_tol)
     Nfinal <- c(); Nnew <- N_start; df <- c()
     Nl <- max(N_start / 2, lb); Nu <- N_start/2+N_start
     if(verbose) cat(paste0("Initiating smart search to find simulation based N for power of ",
                            power_aim, " within ", steps, " steps\nand in total ",
                            Ntotal, " replications. Ns are drawn randomly...\n"))
     for(i in 1:nrow(Conditions))
     {

          # sample sample sizes
          Ns <- sample(seq(Nl, Nu, 1), size = Conditions$Reps[i], replace = T,
                       prob = dnorm(x = seq(-2,2,4/(-1+length(seq(Nl, Nu, 1))))))
          Ns <- sort(Ns, decreasing = T)
          if(verbose) cat(paste0("Step ", i, " of ", steps, ". Fitting ", length(Ns),
                                 " models with Ns in [", Nl, ", ", Nu,"].\n"))

          # simulate data,  fit models and evaluate significance of parameters of interest
          lavModel_attributes <- lavaan::lav_partable_attributes(lavModel)
          matrices <- get_matrices(lavModel, lavModel_attributes)

          cl <- parallel::makeCluster(CORES)
          parallel::clusterExport(cl = cl, varlist = ls(), envir = environment())
          parallel::clusterEvalQ(cl = cl, expr = {
               library(powerNLSEM)
          })
          Fitted <- pbapply::pbsapply(cl = cl,
                                      X = seq_along(Ns), FUN = function(ni) sim_and_fit(n = Ns[ni], POI = POI, alpha = alpha,
                                                                                        lavModel = lavModel,
                                                                                        method = method,
                                                                                        lavModel_Analysis = lavModel_Analysis,
                                                                                        lavModel_attributes = lavModel_attributes,
                                                                                        matrices = matrices,
                                                                                        data_transformations = data_transformations,
                                                                                        prefix = ni),
                                      simplify = T) |> t()
          parallel::stopCluster(cl)

          Sigs <- data.frame(Fitted, Ns); names(Sigs) <- c(colnames(Fitted), "Ns")

          df <- rbind(df, Sigs); rm(Ns)

          ind_min <- which.min(colMeans(Sigs, na.rm = T))# find POI of lowest power
          fit <- glm(df[,ind_min] ~ I(sqrt(Ns)), family = binomial(link = "logit"), data = df)

          Nnew_temp <- round(nlminb(start = 0, objective = find_n_from_glm, fit = fit, pow = power_aim, alpha = alpha)$par)
          if(i <= switchStep)
          {
               Nl_temp <- round(nlminb(start = 0, objective = find_n_from_glm, fit = fit, pow = .15, alpha = 1)$par)
               Nu_temp <- round(nlminb(start = 0, objective = find_n_from_glm, fit = fit, pow = .85, alpha = alpha)$par)
          }else{
               Nl_temp <- round(nlminb(start = 0, objective = find_n_from_glm, fit = fit,
                                       pow = max(power_aim - Conditions$Power_interval[i], .0001), alpha = 1)$par)
               Nu_temp <- round(nlminb(start = 0, objective = find_n_from_glm, fit = fit,
                                       pow = min(power_aim + Conditions$Power_interval[i], .9999), alpha = alpha)$par)
          }

          if(i == (switchStep+1)) # reset process after switch
          {
               Nl <- Nnew - 1; Nu <- Nnew + 1
          }

          # new iteration of Ns
          Nnew <- evaluate_N(N_temp = Nnew_temp, N = Nnew, Sigs = Sigs,
                             ind_min = ind_min, fit = fit, lb = lb,
                             rel_tol = Conditions$Rel_tol[i])
          Nl <- evaluate_N(N_temp = Nl_temp, N = Nl, Sigs = Sigs,
                           ind_min = ind_min, fit = fit, lb = lb,
                           rel_tol = Conditions$Rel_tol[i])
          Nu <- evaluate_N(N_temp = Nu_temp, N = Nu, Sigs = Sigs,
                           ind_min = ind_min, fit = fit, lb = lb,
                           rel_tol = Conditions$Rel_tol[i])
          if(Nu <= Nl){
               Nl <- round(Nu/2)
          }
          Nfinal <- c(Nfinal, Nnew)
     }

     out <- list("N" = Nnew, SigDecisions = df,
                 "N_trials" = Nfinal)

     return(out)
}

# generate sequence of Ns
get_Reps <- function(type = "u", Ntotal = 1000, steps = 10) {
     if(tolower(type) == "increasing")
     {
          temp <- Ntotal/(steps*(steps+1)/2)
          Ns <- seq(temp, temp*steps, temp) |> round()
          temp_diff <- Ntotal - sum(Ns)
          Ns[which.max(Ns)] <- Ns[which.max(Ns)] + temp_diff
     }else if(tolower(type) == "u"){
          temp <- (Ntotal/2)/((steps + steps%%2)/2*((steps + steps%%2)/2+1)/2)
          Ns <- seq(temp, temp*(steps + steps%%2)/2, temp) |> round()
          temp_diff <- Ntotal/2 - sum(Ns)
          Ns[which.max(Ns)] <- Ns[which.max(Ns)] + temp_diff
          Ns <- c(sort(Ns, decreasing = T), Ns)
          if(length(Ns) != steps)
          {
               Ns <- c(Ns[1:((steps+1)/2-1)],
                       sum(Ns[((steps+1)/2):(((steps+1)/2+1))]),
                       Ns[1+((steps+1)/2+1):(steps)])
          }
     }else if(tolower(type) == "equal"){
          temp <- Ntotal/steps
          Ns <- rep(temp, steps) |> round()
          temp_diff <- Ntotal - sum(Ns)
          Ns[which.max(Ns)] <- Ns[which.max(Ns)] + temp_diff
     }
     return(Ns)
}

# find n from an glm-fit model
find_n_from_glm <- function(n, fit, pow = .8, alpha = .05, uncertainty_method = "")
{
     if(uncertainty_method == "exact")
     {
          alpha <- 1 # no influence on se
     }
     logit_fit <- predict(object = fit, newdata = data.frame("Ns" = n), se.fit = T)
     logit_lb <- logit_fit$fit - qnorm(p = 1-alpha/2)*logit_fit$se.fit
     power <- exp(logit_lb)/(1+exp(logit_lb))
     log((power - pow)^2) # use log, otherwise might not minimize
}

# evaluate resonablity of N
evaluate_N <- function(N_temp, N, Sigs, ind_min, fit, lb, rel_tol = .5)
{
     if(any(coef(fit)[-1] < 0))
     {
          N_temp <- ifelse(mean(unlist(Sigs[,ind_min]), na.rm = T) < .5, Inf, -Inf)
     }

     if(abs(N_temp - N) > rel_tol*N){
          N <- round(N + sign(N_temp - N)*rel_tol*N)
     }else{
          N <- N_temp
     }
     # if variation is too low, change sample size drastically
     if(mean(unlist(Sigs[,ind_min]), na.rm = T) > .99)
     {
          N <- round(N*rel_tol+1)
     }else if(mean(unlist(Sigs[,ind_min]), na.rm = T) < .01)
     {
          N <- round(N/rel_tol)
     }
     return(max(lb, N))
}
