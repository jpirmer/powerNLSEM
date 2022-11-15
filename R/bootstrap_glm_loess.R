bootstrap_glm_loess <- function(df, predictors_newdata, criterion_inds, bootstrap, verbose, method, alpha)
{
     Booted_out <- NULL
     if(!is.null(bootstrap))
     {
          pb <- progressr::progressor(along = bootstrap)

          if(verbose) cat("Bootstrapping confidence intervals for unique values predicting power per coefficient.")
          min_max_ind <- unique(c(apply(X = data.frame(predictors_newdata), 2, which.max),
                                  apply(X = data.frame(predictors_newdata), 2, which.min))) # to ensure loess has same support

          # if(method == "logit" | method == "probit")
          # {
          #      Booted <-  glm_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
          #                                    criterion_inds = criterion_inds,
          #                                    indices = c(min_max_ind, sample(x = 1:nrow(df),
          #                                                                    size = nrow(df)-length(min_max_ind),
          #                                                                    replace = T)),
          #                                    link = method, include_est = T)
          # }else if(method == "loess")
          # {
          #      Booted <- loess_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
          #                                     criterion_inds = criterion_inds,
          #                                     indices = c(min_max_ind, sample(x = 1:nrow(df),
          #                                                                     size = nrow(df)-length(min_max_ind),
          #                                                                     replace = T)), include_est = T)
          # }


          Booted <- c() # empty object
          for(i in 1:bootstrap)
          {
               if(method == "logit" | method == "probit")
               {
                    tempBoot <-  glm_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
                                                   criterion_inds = criterion_inds,
                                                   indices = c(min_max_ind, sample(x = 1:nrow(df),
                                                                                   size = nrow(df)-length(min_max_ind),
                                                                                   replace = T)),
                                                   link = method, include_est = T)
               }else if(method == "loess")
               {
                    tempBoot <- loess_power_function(df = df, model_List = model_List, predictors_newdata = predictors_newdata,
                                                   criterion_inds = criterion_inds,
                                                   indices = c(min_max_ind, sample(x = 1:nrow(df),
                                                                                   size = nrow(df)-length(min_max_ind),
                                                                                   replace = T)), include_est = T)
               }
               Booted <- cbind(Booted, unlist(tempBoot[, seq(1, length(criterion_inds)*3, 3)])); #pb()
          };
     }

     if(!is.null(Booted))
     {
          Booted <- apply(Booted, 1, function(x) quantile(x, probs = c(alpha/2, 1-alpha/2))) |> t()
          Booted_out <- c()
          for(i in 1:length(criterion_inds))
          {
               Booted_out <- cbind(Booted_out, Booted[(1:nrow(df))+(i-1)*nrow(df),1:2])
          }
          Booted_out <- data.frame(Booted_out); names(Booted_out) <- names(tempBoot[, -seq(1, length(criterion_inds)*3, 3)])
     }



     return(Booted_out)
}
