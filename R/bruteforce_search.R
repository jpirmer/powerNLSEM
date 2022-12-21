
bruteforce_search <- function(POI,
                              method, lavModel,
                              lavModel_Analysis,
                              data_transformations,
                              N_start = nrow(lavModel_Analysis)*10,
                              Ntotal = 1000,
                              power_aim = .8, alpha = .05,
                              lb = nrow(lavModel_Analysis),
                              CORES,
                              verbose = T,
                              Ns = NULL,
                              ...)
{
     dotdotdot <- list(...)
     if(is.null(Ns)){
          Nl <- max(N_start/4, lb); Nu <- N_start*3
          Ns <- sample(seq(Nl, Nu, 1), size = Ntotal, replace = T,
                       prob = dnorm(x = seq(-2,2,4/(-1+length(seq(Nl, Nu, 1))))))
     }else{
          if(length(table(Ns))==1)
          {
               Ns <- rep(Ns, Ntotal)
          }else{
               Ns <- sample(x = Ns, size = Ntotal, replace = T)
          }
     }; Nl <- min(Ns); Nu <- max(Ns)

     if(verbose) cat(paste0("Initiating brute force search to find simulation based N for power of ",
                            power_aim, " within ",
                            Ntotal, " replications...\n"))
     if(verbose) cat(paste0("Fitting ", length(Ns),
                            " models with Ns in [", Nl, ", ", Nu,"].\n"))


     # simulate data,  fit models and evaluate significance of parameters of interest ----
     df <- c()
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

     if(length(table(df$Ns)) > 1)
     {
          fit <- glm(df[,ind_min] ~ I(sqrt(Ns)), family = binomial(link = "logit"), data = df)
          Nnew_temp <- round(nlminb(start = 0, objective = find_n_from_glm, fit = fit, pow = power_aim, alpha = alpha)$par)
     }else{
          Nnew_temp <- unique(df$Ns)
     }

     # return ----
     out <- list("N" = Nnew_temp, SigDecisions = df,
                 "N_trials" = Nnew_temp)

     return(out)
}
