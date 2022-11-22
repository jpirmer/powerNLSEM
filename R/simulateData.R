#' simulate data from lavModel object
#' @param n sample size
#' @param lavModel lavModel object
#' @param appendLVs logical whether latent variables should be observed. Default to FALSE. (For developmental purposes)
#' @export

simulateData <- function(n, lavModel, appendLVs = F) {

     # construct observed variables
     lavModel_attributes <- lavaan:::lav_partable_attributes(lavModel)
     matrices <- get_matrices(lavModel, lavModel_attributes)
     mat <- matrices[1:4]
     Order <- matrices$Order
     Zeta <- mvtnorm::rmvnorm(n = n, sigma = mat$Psi); colnames(Zeta) <- colnames(mat$Psi)
     Eps <- mvtnorm::rmvnorm(n = n, sigma = mat$Theta); colnames(Eps) <- colnames(mat$Theta)
     Eps <- cbind(Eps, matrix(0, nrow = nrow(Eps), ncol = length(matrices$vnames$ov.iv)))
     colnames(Eps)[colnames(Eps) == ""] <- matrices$vnames$ov.iv

     # start off with "latent" variables
     LV <- Zeta[, Order$lvov[Order$order == 1]]
     for(ord in 2:max(Order$order))
     {
          # construct nonlinear effects and append
          if(unique(Order$type[Order$order == ord]) == "nl")
          {
               tempOrder <- Order[Order==ord,]; NL <- c()
               for(i in 1:nrow(tempOrder))
               {
                    tempVarNames <- stringr::str_split(pattern = ":", tempOrder$lvov[i])[[1]]
                    NL <- cbind(NL, as.matrix(apply(LV[, tempVarNames], 1, prod)))
               }
               colnames(NL) <- tempOrder$lvov
               LV <- cbind(LV, NL)
          }
          # construct deoendent variable and append
          if(unique(Order$type[Order$order == ord]) == "dv")
          {
               tempOrder <- Order[Order==ord,]

               # catch only exisiting variables
               Zeta_temp <- Zeta[, tempOrder$lvov, drop = F]
               Beta_temp <- mat$Beta_dv[tempOrder$lvov, colnames(LV), drop = F]

               # construct DV from existing LV (and ovs) and "residual" Zeta_temp
               DV <- LV %*% t(Beta_temp) + Zeta_temp
               LV <- cbind(LV, DV)
          }
     }

     #  create manifest variables
     Lambda_temp <- mat$Lambda[matrices$vnames$ov, ]
     OV <- LV[, colnames(Lambda_temp)] %*% t(Lambda_temp) + Eps[, matrices$vnames$ov]
     OV <- as.data.frame(OV)
     if(appendLVs){
          LV <- as.data.frame(LV)
          OV <- cbind(OV, LV[, !(names(LV) %in% names(OV))])
     }
     return(OV)
}
