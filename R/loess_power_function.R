loess_power_function <- function(df, model_List, criterion_inds, predictors_newdata, indices = NULL, alpha = .05, include_est = F)
{
     criterion_names <- names(df)[criterion_inds]
     if(is.null(indices)) indices <- 1:nrow(df)
     df_new <- data.frame(df[indices,])
     predicted_loess_se <- lapply(model_List, function(fo){predict(loess(fo, data = df_new),
                                                                   newdata = predictors_newdata,
                                                                   se = T)})
     UpperLowerList <- lapply(seq_along(predicted_loess_se), FUN = function(l){
          loess_se <- predicted_loess_se[[l]]
          est <- loess_se$fit
          lower <- est - qnorm(p = 1-alpha/2)*loess_se$se.fit
          upper <- est + qnorm(p = 1-alpha/2)*loess_se$se.fit
          if(include_est){
               return(cbind(est, lower, upper))
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
