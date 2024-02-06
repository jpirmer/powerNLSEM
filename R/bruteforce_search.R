#' @import stats
#' @import utils
#' @import parallel
#' @importFrom pbapply pbsapply
bruteforce_search <- function(POI, Ns = NULL, N_start = nrow(lavModel_Analysis)*10,
                              method, lavModel,
                              lavModel_Analysis,
                              data_transformations,
                              search_method,
                              power_modeling_method,
                              R = 1000,
                              power_aim = .8, alpha = .05,
                              lb = nrow(lavModel),
                              CORES, verbose = TRUE,
                              alpha_power_modeling = .05,
                              FSmethod = "SL",
                              matchPI =TRUE,
                              PIcentering = "doubleMC",
                              liberalInspection = FALSE,
                              seeds,
                              test = "onesided",
                              ...)
{
     i <- 1; switchStep <- 0; type <- "equal"; steps <- 1 # 1 Trial for now
     sim_seeds <- seeds
     Reps <- get_Reps(type = type, R = R, steps = steps)
     Rel_tol <- Inf
     Power_interval <- 0
     Conditions <- data.frame(Reps, Power_interval, Rel_tol)
     dotdotdot <- list(...)
     if(is.null(Ns)){
          Nl <- max(N_start/4, lb); Nu <- N_start*3
          Ns <- sample(seq(Nl, Nu, 1), size = R, replace = TRUE,
                       prob = dnorm(x = seq(-2,2,4/(-1+length(seq(Nl, Nu, 1))))))
     }else{
          if(length(table(Ns))==1)
          {
               N_start <- Ns; Nnew <- Ns
               Ns <- rep(Ns, R)
          }else if(length(Ns) == R){
               Ns <- Ns
          }else{
               Ns <- sample(x = Ns, size = R, replace = TRUE)
          }
     }; Nl <- min(Ns); Nu <- max(Ns); Nnew <- round(mean(Ns));  Ns <- sort(Ns, decreasing = TRUE)
     Nall <- Ns

     if(verbose) cat(paste0("Initiating brute force search to find simulation based N for power of ",
                            power_aim, " within ",
                            R, " replications...\n"))
     if(verbose) cat(paste0("Fitting ", length(Ns),
                            " models with Ns in [", Nl, ", ", Nu,"].\n"))



     df_POI <- data.frame("matchLabel" = POI)
     fit_temp <- lavModel_Analysis[lavModel_Analysis$matchLabel %in% df_POI$matchLabel,,drop = FALSE]
     # sort by POI
     fit_temp <- merge(df_POI, fit_temp, by.x = "matchLabel", sort = FALSE)
     truth <- fit_temp$ustart; names(truth) <- POI

     Nfinal <- c(); Nnew <- N_start; df <- c()
     df_est <- c(); df_se <- c(); df_pvalue_onesided <- c(); df_pvalue_twosided <- c()
     df_sigs_onesided <- c(); df_sigs_twosided <- c(); vec_fitOK <- c()

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

     df_est <- data.frame(df_est); names(df_est) <- POI
     df_se <- data.frame(df_se); names(df_se) <- POI

     out <- list("N" = Nnew,
            "N_trials" = Nfinal,
            "est" = df_est, "se" = df_se,
            "Ns" = Nall, "fitOK" = vec_fitOK,
            "truth" = truth)

     return(out)
}
