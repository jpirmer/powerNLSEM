getFS <- function(data, manifests, latent, FSmethod, endogene = FALSE)
{
     if(tolower(FSmethod) == "sl")
     {
          method <- ifelse(endogene, "Bartlett", "regression")
     }else{
          method <- FSmethod
     }

     model <- paste0(latent, " =~ ", paste(manifests, collapse = " + "))
     fit <- lavaan::sem(model, data)
     if(lav_object_post_check(fit) & fit@optim$converged)
     {
          FS <- lavaan::lavPredict(fit, method = method)
     }else{
          FS <- matrix(NA, ncol = length(latent), nrow = nrow(data))
          colnames(FS) <- latent
     }
     return(FS)
}
