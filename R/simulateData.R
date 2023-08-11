#' simulate data from lavModel object
#' @importFrom mvtnorm rmvnorm
#' @importFrom stringr str_split
#' @param n sample size
#' @param lavModel lavModel object
#' @param appendLVs logical whether latent variables should be observed. Default to \code{FALSE}. (For developmental purposes)
#' @param lavModel_attributes attributes of the lavModel object. If \code{NULL}, this is computed from lavModel. Default to \code{NULL}.
#' @param matrices computed matrices for simulation. If \code{NULL}, this is computed from lavModel and lavModel_attributes. Default to \code{NULL}.
#' @param seed a seed for reproducability. Default to \code{NULL}.
#' @export

simulateNLSEM <- function(n, lavModel, appendLVs = FALSE, lavModel_attributes = NULL, matrices = NULL, seed = NULL) {

     if(!is.null(seed)) set.seed(seed)

     # construct observed variables
     if(is.null(lavModel_attributes)) lavModel_attributes <- lavaan::lav_partable_attributes(lavModel)
     if(is.null(matrices)) matrices <- get_matrices(lavModel, lavModel_attributes)
     mat <- matrices[1:4]
     Order <- matrices$Order
     Zeta <- mvtnorm::rmvnorm(n = n, sigma = mat$Psi); colnames(Zeta) <- colnames(mat$Psi)
     Eps <- mvtnorm::rmvnorm(n = n, sigma = mat$Theta); colnames(Eps) <- colnames(mat$Theta)
     Eps <- cbind(Eps, matrix(0, nrow = nrow(Eps), ncol = length(matrices$vnames$ov.iv)))
     colnames(Eps)[colnames(Eps) == ""] <- matrices$vnames$ov.iv

     # start off with "latent" variables
     LV <- Zeta[, Order$lvov[Order$order == 1], drop = FALSE]
     for(ord in 2:max(Order$order))
     {
          # construct nonlinear effects and append
          if(unique(Order$type[Order$order == ord]) == "nl")
          {
               tempOrder <- Order[Order==ord,, drop = FALSE]; NL <- c()
               for(i in 1:nrow(tempOrder))
               {
                    tempVarNames <- stringr::str_split(pattern = ":", tempOrder$lvov[i])[[1]]
                    NL <- cbind(NL, as.matrix(apply(LV[, tempVarNames, drop = FALSE], 1, prod)))
               }
               colnames(NL) <- tempOrder$lvov
               LV <- cbind(LV, NL)
          }
          # construct dependent variable and append
          if(unique(Order$type[Order$order == ord]) == "dv")
          {
               tempOrder <- Order[Order==ord,, drop = FALSE]

               # catch only exisiting variables
               Zeta_temp <- Zeta[, tempOrder$lvov, drop = FALSE]
               Beta_temp <- mat$Beta_dv[tempOrder$lvov, colnames(LV), drop = FALSE]

               # construct DV from existing LV (and ovs) and "residual" Zeta_temp
               DV <- LV %*% t(Beta_temp) + Zeta_temp
               LV <- cbind(LV, DV)
          }
     }

     #  create manifest variables
     Lambda_temp <- mat$Lambda[matrices$vnames$ov, , drop = FALSE]
     OV <- LV[, colnames(Lambda_temp), drop = FALSE] %*% t(Lambda_temp) + Eps[, matrices$vnames$ov, drop = FALSE]
     OV <- as.data.frame(OV)
     if(appendLVs){
          LV <- as.data.frame(LV)
          OV <- cbind(OV, LV[, !(names(LV) %in% names(OV)), drop = FALSE])
     }
     return(OV)
}
