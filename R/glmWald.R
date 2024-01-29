# These scripts fit a probit model with sqrt(n) as predictor for the two sided test
# which is equivalent to the one-parameter Wald test, i.e., this is a probit-regression for
# Wald-test distributed test results
LL_WaldProbit <- function(beta, Ns, sig, vec = FALSE)
{
     P <- pnorm(q = -beta[1] + beta[2]*sqrt(Ns))  +
          pnorm(q = beta[1] +  beta[2]*sqrt(Ns),
                lower.tail = FALSE)

     LL <- log(P)*sig + (1-sig)*log(1-P)
     if(vec) return(LL)
     return(-sum(LL))
}

fitWaldglm <- function(sig, Ns)
{
     betas <- abs(coef(glm(sig~sqrt(Ns), family = binomial("probit"))))
     trial <- 0
     repeat{
          start <- abs(betas + runif(n = 2, min = -.1, max = .1)*c(1/5, 1/10)*trial)
          est <- try(nlminb(start = start, lower = c(0, 0),
                            objective = LL_WaldProbit, Ns = Ns, sig = sig,
                            vec = FALSE,
                            control = list(rel.tol = 10^-12,
                                           abs.tol = 10^-30, eval.max = 10^3)),
                     silent = TRUE)
          if(!(inherits(est, "try-error") | est$convergence == 1)) break
          if(trial > 20) break
          trial <- trial + 1
     }

     if(inherits(est, "try-error") | est$convergence == 1){
          beta <- c(NA, NA)
          P <- rep(NA, length(Ns))
     }else{
          beta <- est$par
          P <- pnorm(q = -beta[1] + beta[2]*sqrt(Ns))  +
               pnorm(q = beta[1] +  beta[2]*sqrt(Ns),
                     lower.tail = FALSE)
     }
     out <- list("est" = beta, "P" = P, "sig" = sig, "Ns" = Ns)
     class(out) <- "WaldGLM"
     return(out)
}

predict_WaldGLM <- function(out, N_interest)
{
     if(class(out) != "WaldGLM") stop("Predict for WaldGLM object.")
     beta <- out$est
     P <- pnorm(q = -beta[1] + beta[2]*sqrt(N_interest))  +
          pnorm(q = beta[1] +  beta[2]*sqrt(N_interest),
                lower.tail = FALSE)
     return(P)
}

Wald_pred_confint <- function(out, N_interest, alpha)
{
     if(class(out) != "WaldGLM") stop("Confint for WaldGLM object.")
     beta <- out$est
     sig <- out$sig
     Ns <- out$Ns

     X <- cbind(1, sqrt(N_interest))
     Xm <- cbind(-1, sqrt(N_interest))
     Xmb <- Xm %*% beta
     Xb <- X %*% beta

     H <- 1/length(Ns)*numDeriv::hessian(func = LL_WaldProbit, x = beta,
                                         sig = sig, Ns = Ns, vec = FALSE)
     I <- cov(numDeriv::jacobian(func = LL_WaldProbit, x = beta,
                                 sig = sig, Ns = Ns, vec = TRUE))
     invI <- (solve(H) %*% I %*% solve(H)) / length(Ns)

     VarxHx <- apply(X, MARGIN = 1, function(x) t(x) %*% invI %*% as.matrix(x))
     VarxHx_minus <- apply(Xm, MARGIN = 1, function(x) t(x) %*% invI %*% as.matrix(x))

     Xmb_lb <- Xmb - qnorm(1-alpha/2)*sqrt(VarxHx_minus)
     Xmb_ub <- Xmb + qnorm(1-alpha/2)*sqrt(VarxHx_minus)

     Xb_lb <- Xb - qnorm(1-alpha/2)*sqrt(VarxHx)
     Xb_ub <- Xb + qnorm(1-alpha/2)*sqrt(VarxHx)

     P <-  pnorm(q = Xmb)  + pnorm(q = Xb, lower.tail = FALSE)
     P_lb <- pnorm(q = Xmb_lb)  + pnorm(q = Xb_ub, lower.tail = FALSE)
     P_ub <- pnorm(q = Xmb_ub)  + pnorm(q = Xb_lb, lower.tail = FALSE)

     return(data.frame("P" = P, "P_lb" = P_lb, "P_ub" = P_ub))
}
