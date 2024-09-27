#' powerNLSEM function
#' @import stats
#' @import parallel
#' @importFrom stringr str_replace_all
#' @param model Model in lavaan syntax. See documentation for help and examples.
#' @param POI Parameter Of Interest as a vector of strings. Must be in lavaan-syntax without any spaces. Nonlinear effects should have the same ordering as in model.
#' @param method Method used to fit to the data. Implemented methods are \code{"LMS"} (Klein & Moosbrugger, 2000) (requires an installation of \code{Mplus} and the \code{MplusAutomation} pacakge), \code{"UPI"} (Kelava & Brandt, 2009, Marsh et al., 2004) for the unconstrained product indicator approach, \code{"FSR"} (Ng and Chan, 2020) for the naïve factor score approach, and \code{"SR"}, for using scale means (i.e., scale regression/path modeling).
#' @param test Should the parameter be tested with a directed hypothesis (onesided) or with an undirected hypothesis (twosided, also equivalent to Wald-Test for single parameter). Default to \code{"onesided"}.
#' @param power_modeling_method Power modeling method used to model significant parameter estimates. Default to \code{"probit"} indicating glm with probit link function with sqrt(n) as predictor. Alternative is \code{"logit"}.
#' @param search_method String stating the search method. Default to \code{"adaptive"} (synonyme is \code{"smart"}). Alternative is \code{"bruteforce"}.
#' @param R Total number of models to be fitted. Higher number results in higher precision and longer runtime. Default to 2000.
#' @param power_aim Minimal power value to approximate. Default to \code{.8}.
#' @param alpha Type I-error rate for significance decision. Default to \code{.05}.
#' @param alpha_power_modeling Type I-error rate for confidence band around predicted power rate. Used to ensure that the computed \code{N} keeps the desired power value (with the given Type I-error rate \code{alpha_power_modeling} divided by 2). If set to 1, no confidence band is used. Default to \code{.05}.
#' @param CORES Number of cores used for parallelization. Default to number of available cores - 2.
#' @param verbose Logical whether progress should be printed in console. Default to \code{TRUE}.
#' @param seed Seed for replicability. Default to \code{NULL}, then a seed is drawn at random, which will also be saved in the output.
#' @param ... Additional arguments passed on to the search functions.
#' @seealso For further details for specific uses see corresponding functions: [power_search()] for all inputs possible, [UPI()] for specifics for the unconstrained product indicator approach, [LMS()] for the latent moderated structural equations approach, [FSR()] for factor score approaches, [SR()] for scale regression approaches.
#' @returns Returns an list object of class \code{powerNLSEM}.
#' @references Klein, A. G., & Moosbrugger, H. (2000). Maximum likelihood estimation of latent interaction effects with the LMS method. _Psychometrika, 65_(4), 457–474. \doi{10.1007/BF02296338}
#' @references Kelava, A., & Brandt, H. (2009). Estimation of nonlinear latent structural equation models using the extended unconstrained approach. _Review of Psychology, 16_(2), 123–132.
#' @references Lin, G. C., Wen, Z., Marsh, H. W., & Lin, H. S. (2010). Structural equation models of latent interactions: Clarification of orthogonalizing and double-mean-centering strategies. _Structural Equation Modeling, 17_(3), 374–391. \doi{10.1080/10705511.2010.488999}
#' @references Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of orthogonalizing powered and product terms: Implications for modeling interactions among latent variables. _Structural Equation Modeling, 13_(4), 497–519. \doi{10.1207/s15328007sem1304_1}
#' @references Marsh, H. W., Wen, Z. & Hau, K. T. (2004). Structural equation models of latent interactions: Evaluation of alternative estimation strategies and indicator construction. _Psychological Methods, 9_(3), 275–300. \doi{10.1037/1082-989X.9.3.275}
#' @references Ng, J. C. K., & Chan, W. (2020). Latent moderation analysis: A factor score approach. _Structural Equation Modeling: A Multidisciplinary Journal, 27_(4), 629–648. \doi{10.1080/10705511.2019.1664304}.
#' @references Irmer, J. P., Klein, A. G., & Schermelleh-Engel, K. (2024a). A General Model-Implied Simulation-Based Power Estimation Method for Correctly and Misspecfied Models: Applications to Nonlinear and Linear Structural Equation Models. _Behavior Research Methods._  \doi{10.31219/osf.io/pe5bj}
#' @references Irmer, J. P., Klein, A. G., & Schermelleh-Engel, K. (2024b). Estimating Power in Complex Nonlinear Structural Equation Modeling Including Moderation Effects: The `powerNLSEM R`-Package. _Behavior Research Methods._ \doi{10.3758/s13428-024-02476-3}
#' @examples
#' \donttest{
#' # write model in lavaan syntax
#' model <- "
#' # measurement models
#'           X =~ 1*x1 + 0.8*x2 + 0.7*x3
#'           Y =~ 1*y1 + 0.85*y2 + 0.78*y3
#'           Z =~ 1*z1 + 0.9*z2 + 0.6*z3
#'
#' # structural models
#'           Y ~ 0.3*X + .2*Z +  .2*X:Z
#'
#' # residual variances
#'          Y~~.7975*Y
#'          X~~1*X
#'          Z~~1*Z
#'
#' # covariances
#'          X~~0.5*Z
#'
#' # measurement error variances
#'          x1~~.1*x1
#'          x2~~.2*x2
#'          x3~~.3*x3
#'          z1~~.2*z1
#'          z2~~.3*z2
#'          z3~~.4*z3
#'          y1~~.5*y1
#'          y2~~.4*y2
#'          y3~~.3*y3
#' "
#' # run model-implied simulation-based power estimation
#' # for the effects: c("Y~X", "Y~Z", "Y~X:Z")
#' Result_Power <- powerNLSEM(model = model, POI = c("Y~X", "Y~Z", "Y~X:Z"),
#'                            method = "UPI", search_method = "adaptive",
#'                            steps = 10, power_modeling_method = "probit",
#'                            R = 1000, power_aim = .8, alpha = .05,
#'                            alpha_power_modeling = .05,
#'                            CORES = 1, seed = 2024)
#'
#' Result_Power
#' }
#'
#' @export

powerNLSEM <- function(model, POI,
                       method,
                       test = "onesided",
                       power_modeling_method = "probit",
                       search_method = "adaptive",
                       R = 2000,
                       power_aim = .8,
                       alpha = .05,
                       alpha_power_modeling = .05,
                       CORES = max(c(parallel::detectCores()-2, 1)),
                       verbose = TRUE, seed = NULL,
                       ...)
{
     call <- match.call()

     if(!is.null(seed)){
          set.seed(seed)
          seeds <- sample(1:10^9, size = R)
     }else{
          seeds <- NULL
     }
     t0 <- proc.time()

     # check whether needed packages are installed
     if(tolower(method) %in% c("sem", "upi")) rlang::check_installed("semTools")
     if(tolower(method) == "lms") rlang::check_installed("MplusAutomation")

     ### prepare model ----
     lavModel <- lavaan::lavMatrixRepresentation(lavaan::lavaanify(model = model,
                                                                    meanstructure = TRUE,auto.var = TRUE,
                                                                    auto.cov.lv.x = TRUE, auto.cov.y = TRUE,
                                                                    as.data.frame. = TRUE))
     lavModel <- add_varType(lavModel)
     temp <- check_model(lavModel); lavModel <- temp$lavModel; added_model_syntax <- temp$added_model_syntax
     Manifests <- handle_manifests(lavModel = lavModel, treat_manifest_as_latent = "ov")
     lavModel_Analysis <- Manifests$lavModel_Analysis
     data_transformations <- Manifests$data_transformations
     lavModel_Analysis$matchLabel <- stringr::str_replace_all(paste0(lavModel_Analysis$lhs,
                                             lavModel_Analysis$op,
                                             lavModel_Analysis$rhs), pattern = "_", replacement = ":")


     ### initialize ----
     dotdotdot <- list(...)
     if(!is.null(dotdotdot$Ns)){
          Ns <- dotdotdot$Ns
          search_method <- "bruteforce"
     }else{
          Ns <- NULL
     }
     if(!is.null(dotdotdot$N_start)){
          N_start <- dotdotdot$N_start
     }else{
          N_start <- nrow(lavModel[lavModel$op != "~1", ])*10
     }
     if(!is.null(dotdotdot$nlb)){
          nlb <- max(dotdotdot$nlb, nrow(lavModel))
     }else{
          nlb <-  nrow(lavModel[lavModel$op != "~1", ])*5
     }
     if(!is.null(dotdotdot$steps)){
          steps <- dotdotdot$steps
     }else{
          steps <- 10
     }
     if(!is.null(dotdotdot$switchStep)){
          switchStep <- dotdotdot$switchStep
     }else{
          switchStep <- round(steps/2)
     }
     if(!is.null(dotdotdot$distRj)){
          distRj <- dotdotdot$distRj
     }else{
          distRj <- "increasing"
     }
     if(!is.null(dotdotdot$FSmethod)){
          FSmethod <- dotdotdot$FSmethod
     }else{
          FSmethod <- "SL"
     }
     if(!is.null(dotdotdot$matchPI)){
          matchPI <- dotdotdot$matchPI
     }else{
          matchPI <- TRUE
     }
     if(!is.null(dotdotdot$PIcentering)){
          PIcentering <- dotdotdot$PIcentering
     }else{
          PIcentering <- "doubleMC"
     }
     if(!is.null(dotdotdot$liberalInspection)){
          liberalInspection <- dotdotdot$liberalInspection
     }else{
          liberalInspection <- FALSE
     }
     if(!is.null(dotdotdot$constrainRelChange)){
          constrainRelChange <- dotdotdot$constrainRelChange
     }else{
          constrainRelChange <- TRUE
     }
     if(!is.null(dotdotdot$pathLMS)){
          pathLMS <- dotdotdot$pathLMS
     }else{
          pathLMS <- tempdir()
     }

     # check input ------
     # check plausibility of input
     if(!(tolower(method) %in% c("lms", "upi", "sem", "fsr", "factorscores", "sr", "reg")))
          stop("Wrong input in 'method', should be 'LMS', 'UPI', 'SEM', 'FSR', or 'SR'.")
     if(!(tolower(test) %in% c("onesided", "twosided")))
          stop("Wrong input in 'test', should be 'onesided' or 'twosided'.")
     if(!(tolower(power_modeling_method) %in% c("probit", "wald", "logit")))
          stop("Wrong input in 'power_modeling_method', should be 'probit', 'wald', or 'logit'.")
     if(!(tolower(search_method) %in% c("smart", "bruteforce", "adaptive")))
          stop("Wrong input in 'test', should be 'adaptive' (or 'smart'), or 'bruteforce'.")
     if(R<0 | round(R) != R) stop("'R' needs to a natural number.")
     if(power_aim >= 1 | power_aim <= 0) stop("'power_aim' needs to be within (0, 1).")
     if(alpha >= 1 | alpha <= 0) stop("'alpha' needs to be within (0, 1).")
     if(alpha_power_modeling >= 1 | alpha_power_modeling <= 0) stop("'alpha_power_modeling' needs to be within (0, 1).")


     # check model and modeling approaches
     if(!("=~"  %in% lavModel$op)){
          if(!(tolower(method) %in% c("upi", "sem"))){
               warning("No latent variables present in model, advised to use method = 'sem' (or method = 'UPI', as they are identical).")
          }
     }
     # check cross-relations
     temp <- lavModel[lavModel$op == "=~",]
     if(nrow(temp) > 0L)
     {
          tempList <- lapply(unique(temp$lhs), FUN = function(x) temp[temp$lhs == x, ]$rhs)
          crossloadings <- FALSE
          for(i in 1:(length(tempList)-1))
          {
               mani1 <- tempList[[i]]
               if(any(sapply(tempList[-i], FUN = function(x) any(mani1 %in% x)))){
                    crossloadings <- TRUE
               }
          }
          if(crossloadings)
          {
               if(tolower(method) == "sr") warning("Cross-loadings influence the construction of scale scores.\nReconsider model!")
               if(tolower(method) == "upi") warning("Cross-loadings influence the construction of product indicators. There is not much research on this!\nReconsider model!")
               if(tolower(method) == "fsr") stop("Cross-loadings influence the construction of factor scores.\nCurrently implemented version of FSR() constructs factor scores latent-variable-wise, hence, cross-loadings are not taken into account, correctly.\nReconsider model!")
          }

          # residual covariances
          temp <- lavModel[lavModel$mat == "theta",]
          ResidCovMat <- matrix(FALSE, nrow = max(c(temp$row, temp$col)),
                                ncol = max(c(temp$row, temp$col)))
          for(i in 1:nrow(temp))
          {
               ResidCovMat[temp$row[i], temp$col[i]] <- ResidCovMat[temp$col[i], temp$row[i]] <- TRUE
          }
          diag(ResidCovMat) <- FALSE
          if(any(ResidCovMat))
          {
               if(tolower(method) == "sr") warning("Residual correlation influences the construction of scale scores.\nReconsider model!")
               if(tolower(method) == "upi") stop("Residual correlation influences the construction of product indicators.\nCurrently implemented version of UPI() does not handle residual correlations, yet!\nReconsider model!")
               if(tolower(method) == "fsr") stop("Residual correlation influences the construction of factor scores.\nCurrently implemented version of FSR() constructs factor scores latent-variable-wise, hence, cross-correlations are not taken into account, correctly.\nReconsider model!")
          }
     }


     # labels and parameter constraints
     if(any(lavModel$op == ":=") | any(lavModel$label != ""))
     {
          stop("Labels and parameter constrains are not implemented (yet).\nPlease reconsider your model!\nPower is expected to be smaller for more parameters, hence, the neglection of paramter constrains should be more conservative.")
     }

     # multiple groups
     if(any(lavModel$group != 1L)) stop("Multiple group analysis has not been implemented yet for any method!")

     # check input with regard to power modeling
     if(!(tolower(power_modeling_method) %in% c("probit", "wald"))) warning("Probit and/or Wald are the theoretical link between sqrt(n) and power.")
     if(test == "twosided"){
          if(power_modeling_method == "probit") warning("Probit can be used with test = 'twosided', but power_modeling_method = 'Wald' is adviced.")
          if(power_modeling_method == "logit") warning("Logit should NOT be used with test = 'twosided', but power_modeling_method = 'Wald' is adviced.")
     }
     if(test == "onesided"){
          if(power_modeling_method == "wald") warning("Wald can be used with test = 'onesided', but power_modeling_method = 'probit' is adviced.")
          if(power_modeling_method == "logit") warning("Logit should NOT be used with test = 'onesided', but power_modeling_method = 'probit' is adviced.")
     }

     ### begin ------
     POI <- stringr::str_replace_all(string = POI, pattern = " ", replacement = "")

     args <- names(formals(power_search))
     args <- args[args!="..."]

     if(!all(POI %in% lavModel_Analysis$matchLabel)){
          POI_missing <- POI[(!POI %in% lavModel_Analysis$matchLabel)]
          stop(paste0("Parameters Of Interest (POI) not given in model. Includes: ",
                      paste(POI_missing, collapse = ", ", sep = ""), "!"))
     }

     ### run power analysis ----
     out <- do.call("power_search", mget(args))
     calledArgs <- mget(args)

     # examine fitting performance via weighted average Bias, Relative Bias and RWMSE
     TruthMat <- matrix(rep(out$truth, each = R),
                            ncol = length(POI))
     Bias <- (out$est - TruthMat)[out$fitOK, ]
     RelBias <- Bias / TruthMat[out$fitOK, ]
     Ns <- out$Ns[out$fitOK]
     AvgWeightedBias <- t(as.matrix(Bias)) %*% as.matrix(Ns) /
          sum(Ns)
     AvgWeightedRelBias <- t(as.matrix(RelBias)) %*% as.matrix(Ns) /
          sum(Ns)
     AvgWeightedAbsBias <- t(as.matrix(abs(Bias))) %*% as.matrix(Ns) /
          sum(Ns)
     AvgWeightedRWMSE <- sqrt(t(as.matrix(Bias^2)) %*% as.matrix(Ns) /
                                  sum(Ns))

     Performance <- data.frame(rbind(t(AvgWeightedBias), t(AvgWeightedAbsBias),
                                     t(AvgWeightedRelBias), t(AvgWeightedRWMSE)))
     rownames(Performance) <- c("Bias", "absolute Bias", "relative Bias", "RWMSE")
     colnames(Performance) <- POI
     AveragePerformance <- rowMeans(Performance)

     if(mean(out$fitOK) < .5) warning("More than half of models did not converge.\nResults are not trustworthy.\nCheck your model formulation!")

     ### return results ----
     t <- proc.time()-t0
     out$power <- power_aim
     out$beta <- 1-power_aim
     out$alpha <- alpha
     out$alpha_power_modeling <- alpha_power_modeling
     out$method <- method
     out$search_method <- search_method
     out$power_modeling_method <- power_modeling_method
     out$test <- test
     out$convergenceRate <- mean(out$fitOK)
     out$Performance <- Performance
     out$AveragePerformance <- AveragePerformance
     out$seed <- seed
     out$model <- paste0(model, added_model_syntax, collapse = "\n")
     out$runtime <- t
     out$call <- call
     out$args <- calledArgs
     class(out) <- c("powerNLSEM", "list")
     return(out)
}
