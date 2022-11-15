
est_power_glm <- function(df, predictor_inds, criterion_inds, bootstrap, verbose, alpha) {

  # unique predictors for predicted values (if bootstrapping is used)
  predictors_newdata <- data.frame(df[1:nrow(df), predictor_inds])
  crit <- colnames(df[, criterion_inds])
  pred <- colnames(df[, predictor_inds])
  model_List <- sapply(crit, function(x) paste0(x, "~", paste(pred, collapse = "+", sep="")))


  glm_power_function <- function(df, model_List, predictors_newdata, indices = NULL, alpha = .05)
  {
       if(is.null(indices)) indices <- 1:nrow(df)
       df_new <- data.frame(df[indices,])
       predicted_logit_se <- lapply(model_List, function(fo){predict(glm(fo, data = df_new,
                                                                         family = binomial(link = "logit")),
                                                                     newdata = predictors_newdata,
                                                                     se.fit = T)})
       UpperLowerList <- lapply(seq_along(predicted_logit_se), FUN = function(l){
            logit_se <- predicted_logit_se[[l]]
            logit_lower <- logit_se$fit - qnorm(p = 1-alpha/2)*logit_se$se.fit
            logit_upper <- logit_se$fit + qnorm(p = 1-alpha/2)*logit_se$se.fit
            lower <- exp(logit_lower)/(1 + exp(logit_lower))
            upper <- exp(logit_upper)/(1 + exp(logit_upper))
            cbind(lower, upper)
       })

       UpperLower <- data.frame(matrix(unlist(UpperLowerList), ncol = length(UpperLowerList)*2))
       names(UpperLower) <- sapply(colnames(criterions), function(x) paste0("p_",x,c("_lb", "_ub")))
       return(UpperLower)
  }




  if(!is.null(bootstrap))
  {
       if(verbose) cat("Bootstrapping confidence intervals for unique values predicting power per coefficient.")
       Booted <- glm_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
                                    indices = sample(x = 1:nrow(df), size = nrow(df), replace = T), alpha = alpha)
       Booted[,] <- 0 # empty object
       for(i in 1:bootstrap)
       {
            tempBoot <- glm_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
                                           indices = sample(x = 1:nrow(df), size = nrow(df), replace = T), alpha = alpha)
            Booted <- Booted + tempBoot
       }; Booted <- Booted/bootstrap
  }

}
