
est_power <- function(df, predictor_inds, criterion_inds, bootstrap, verbose, alpha, method = "logit")
{
     # unique predictors for predicted values (if bootstrapping is used)
     predictors_newdata <- data.frame(df[1:nrow(df), predictor_inds, drop = F])
     crit <- colnames(df[, criterion_inds, drop = F])
     pred <- colnames(df[, predictor_inds, drop = F])
     model_List <- sapply(crit, function(x) paste0(x, "~", paste(pred, collapse = "+", sep="")))


     if(method == "logit" | method == "probit")
     {
          power <-  glm_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
                                       indices =1:nrow(df), alpha = alpha,
                                       criterion_inds = criterion_inds, link = method, include_est = T)

          progressr::with_progress(Booted <- bootstrap_glm_loess(df = df, predictors_newdata = predictors_newdata, criterion_inds = criterion_inds,
                                                                 bootstrap = bootstrap, verbose = verbose, method = method, alpha = alpha))

     }else if(method == "loess"){
          power <-  loess_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
                                       indices = 1:nrow(df), alpha = alpha,
                                       criterion_inds = criterion_inds, include_est = T)

          progressr::with_progress(Booted <- bootstrap_glm_loess(df = df, predictors_newdata = predictors_newdata, criterion_inds = criterion_inds,
                                                                 bootstrap = bootstrap, verbose = verbose, method = method, alpha = alpha))
     }


     out <- list("power" = power[, seq(1, length(criterion_inds)*3, 3)], "CI" = power[, -seq(1, length(criterion_inds)*3, 3)],
                 "BootCI" = Booted)
     return(out)

}
