
bruteforce_search <- function(POI,
                              method, lavModel,
                              lavModel_Analysis,
                              data_transformations,
                              power_modeling_method,
                              N_start = nrow(lavModel_Analysis)*10,
                              Ntotal = 1000,
                              power_aim = .8, alpha = .05,
                              lb = nrow(lavModel_Analysis),
                              CORES,
                              verbose = T,
                              Ns = NULL,
                              uncertainty_method = "",
                              ...)
{
     i <- 1; switchStep <- 0; type <- "equal"; steps <- 1 # 1 Trial for now
     Reps <- get_Reps(type = type, Ntotal = Ntotal, steps = steps)
     # Power_interval <- c(rep(0,switchStep), seq(.1, .01, -(.1-.01)/(steps-switchStep-1)))
     # Rel_tol <- seq(.5, 1, (1-.5)/(steps-1))
     Rel_tol <- 2
     Power_interval <- .1
     Conditions <- data.frame(Reps, Power_interval, Rel_tol)
     dotdotdot <- list(...)
     if(is.null(Ns)){
          Nl <- max(N_start/4, lb); Nu <- N_start*3
          Ns <- sample(seq(Nl, Nu, 1), size = Ntotal, replace = T,
                       prob = dnorm(x = seq(-2,2,4/(-1+length(seq(Nl, Nu, 1))))))
     }else{
          if(length(table(Ns))==1)
          {
               N_start <- Ns; Nnew <- Ns
               Ns <- rep(Ns, Ntotal)
          }else if(length(Ns) == Ntotal){
               Ns <- Ns
          }else{
               Ns <- sample(x = Ns, size = Ntotal, replace = T)
          }
     }; Nl <- min(Ns); Nu <- max(Ns); Nnew <- round(mean(Ns));  Ns <- sort(Ns, decreasing = T)


     if(verbose) cat(paste0("Initiating brute force search to find simulation based N for power of ",
                            power_aim, " within ",
                            Ntotal, " replications...\n"))
     if(verbose) cat(paste0("Fitting ", length(Ns),
                            " models with Ns in [", Nl, ", ", Nu,"].\n"))


     # simulate data,  fit models and evaluate significance of parameters of interest ----
     df <- c()
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
                                                                                   prefix = ni),
                                 simplify = T) |> t()
     if(CORES > 1L) parallel::stopCluster(cl)

     Sigs <- data.frame(Fitted, Ns); names(Sigs) <- c(colnames(Fitted), "Ns")

     df <- rbind(df, Sigs); rm(Ns)

     ind_min <- which.min(colMeans(df, na.rm = T))# find POI of lowest power

     ### run power model
     args <- names(formals(fit_power_model))
     args <- args[args!="..."]
     N_temp <- do.call("fit_power_model", mget(args))
     Nnew <- N_temp$Nnew

     # return ----
     out <- list("N" = Nnew, SigDecisions = df,
                 "N_trials" = Nnew)

     return(out)
}
