glm_power_function <- function(df, model_List, predictors_newdata, criterion_inds, indices = NULL, alpha = .05, link = "logit", include_est = F)
{
     criterion_names <- names(df)[criterion_inds]
     if(is.null(indices)) indices <- 1:nrow(df)
     df_new <- data.frame(df[indices,])
     predicted_glm_se <- lapply(model_List, function(fo){predict(glm(fo, data = df_new,
                                                                       family = binomial(link = link)),
                                                                   newdata = predictors_newdata,
                                                                   se.fit = T)})
     UpperLowerList <- lapply(seq_along(predicted_glm_se), FUN = function(l){
          glm_se <- predicted_glm_se[[l]]

          if(link == "logit")
          {
               est_logit <- glm_se$fit
               logit_lower <- est_logit - qnorm(p = 1-alpha/2)*glm_se$se.fit
               logit_upper <- est_logit + qnorm(p = 1-alpha/2)*glm_se$se.fit
               est <- exp(est_logit)/(1 + exp(est_logit))
               lower <- exp(logit_lower)/(1 + exp(logit_lower))
               upper <- exp(logit_upper)/(1 + exp(logit_upper))
          }else if(link == "probit"){
               est_probit <- glm_se$fit
               probit_lower <- glm_se$fit - qnorm(p = 1-alpha/2)*glm_se$se.fit
               probit_upper <- glm_se$fit + qnorm(p = 1-alpha/2)*glm_se$se.fit
               est <- pnorm(est_probit)
               lower <- pnorm(probit_lower)
               upper <- pnorm(probit_upper)
          }
          if(include_est){
               return(cbind( est, lower, upper))
          }else{
               return(cbind(lower, upper))
          }
     })
     if(include_est)
     {
          UpperLower <- data.frame(matrix(unlist(UpperLowerList), ncol = length(UpperLowerList)*3))
          names(UpperLower) <- sapply(criterion_names, function(x) paste0("p_",x,c("","_lb", "_ub"))) |> c()
     }else{
          UpperLower <- data.frame(matrix(unlist(UpperLowerList), ncol = length(UpperLowerList)*2))
          names(UpperLower) <- sapply(criterion_names, function(x) paste0("p_",x,c("_lb", "_ub"))) |> c()
     }
     return(UpperLower)
}
