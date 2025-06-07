#' Extract and Visualize GW Regression Model Residuals
#'
library(doParallel)
library(GWRM)

source("pgw.R")
source("dgw.R")

residuals.gw <- function (object, type = c("pearson", "response", "quantile"), rep = 19, envelope = FALSE, 
          title = "Simulated Envelope of Residuals", trace = FALSE, 
          parallel = TRUE, ncores = 2, ...) {
  
  type <- match.arg(type)
  
  resid.gw <- function(object, type) {
    
    #defino función necesaria
    pgw <- function (q, a, k, ro, lower.tail = TRUE) 
    {
      if (!(is.double(q) || is.integer(q)) || !(is.double(a) || 
                                                is.integer(a)) || !(is.double(k) || is.integer(k))|| 
          !(is.double(ro) || is.integer(ro))) 
        stop("Non-numeric argument to mathematical function")
      maximum <- max(length(q), length(a), length(k), length(ro))
      q <- rep(q, length.out = maximum)
      a <- rep(a, length.out = maximum)
      k <- rep(k, length.out = maximum)
      ro <- rep(ro, length.out = maximum)
      prob <- numeric(length = maximum)
      for (ind in seq_along(q)) {
        if (q[ind] < 0) {
          prob[ind] <- 0
          next
        }
        qlower <- trunc(q[ind])
        x <- 0:qlower
        probs <- dgw(x, a[ind], k[ind], ro[ind])
        if (lower.tail) {
          prob[ind] <- sum(probs)
        }
        else {
          prob[ind] <- 1 - sum(probs)
        }
      }
      return(prob)
    }
    
    #defino función necesaria
    dgw <- function(x, a, k, ro){
      
      if (mode(c(x, a, k, ro)) != "numeric") 
        stop("non-numeric argument to mathematical function")
      if (any(c(a, k, ro) <= 0)) 
        stop("a, k and rho must be greater than 0")
      
      lpmf <- lgamma(a + ro) + lgamma(k + ro) + lgamma(a + x) + lgamma(k + x) - lgamma(ro) - 
        lgamma(a) - lgamma(k) - lgamma(a + k + ro + x) - lgamma(x + 1)
      pmf <- exp(lpmf)
      return(pmf)
      
    }
    
    mu <- object$fitted.values
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
    a <- mu * (ro - 1)/k
    if (type == "pearson") {
      if (ro < 2) 
        stop("Variance is infinite")
      else {
        varianza <- ((k + ro - 1)/(ro - 2)) * (mu + (mu^2)/k)
        residuos <- sort(object$residuals/sqrt(varianza))
        return(residuos)
      }
    }
    if (type == "response") {
      residuos <- sort(object$residuals)
      return(residuos)
    }
    if (type == "quantile") {
      formula <- stats::as.formula(object$call[["formula"]])
      matriz.datos <- stats::model.frame(formula, data = object$data)
      y <- stats::model.response(matriz.datos)
      A <- ifelse(y == 0, 0, pgw(y - 1, a, k, ro))
      B <- pgw(y, a, k, ro)
      u <- stats::runif(length(y), A, B)
      qr <- stats::qnorm(u)
      return(sort(qr))
    }
  }
    
  if (envelope == FALSE) 
    return(resid.gw(object, type))
  else {
    st <- proc.time()
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
    mu <- object$fitted.values
    a <- mu * (ro - 1)/k
    n <- sum(object$W)
    object$data <- na.omit(object$data[all.vars(object$formula)])
    if (parallel) {
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      clusterExport(cl, list("object", "rgw", "type"), 
                    envir = environment())
      residuos.sim <- foreach(j = 1:rep, .combine = cbind, 
                              .multicombine = TRUE, .inorder = FALSE, .packages = c("GWRM"), 
                              .verbose = as.logical(trace)) %dopar% {
                                converged <- FALSE
                                varResponse <- GWRM:::getResponse(object$formula)
                                datos <- object$data[rep(1:nrow(object$data), 
                                                         object$W), ]
                                datos[varResponse] <- as.matrix(GWRM::rgw(n, a, k, 
                                                                    ro))
                                while (!converged) {
                                  fit <- try(GWRM::gw(object$formula, data = datos, 
                                                      k = object$k), silent = TRUE)
                                  if (fit$aic > 0 && fit$betaIIpars[2] > 2) 
                                    converged <- TRUE
                                  else {
                                    datos[varResponse] <- as.matrix(GWRM::rgw(n, a, 
                                                                        k, ro))
                                  }
                                }
                                as.matrix(resid.gw(fit, type))
                              }
      stopCluster(cl)
    }
    else {
      residuos.sim <- foreach(j = 1:rep, .combine = cbind, 
                              .multicombine = TRUE, .inorder = FALSE, .packages = c("GWRM"), 
                              .verbose = as.logical(trace)) %do% {
                                converged <- FALSE
                                varResponse <- GWRM:::getResponse(object$formula)
                                datos <- object$data[rep(1:nrow(object$data), 
                                                         object$W), ]
                                datos[varResponse] <- as.matrix(GWRM::rgw(n, a, k, 
                                                                    ro))
                                while (!converged) {
                                  fit <- try(GWRM::gw(object$formula, data = datos, 
                                                      k = object$k), silent = TRUE)
                                  if (fit$aic > 0 && fit$betaIIpars[2] > 2) 
                                    converged <- TRUE
                                  else {
                                    datos[varResponse] <- as.matrix(GWRM::rgw(n, a, 
                                                                        k, ro))
                                  }
                                }
                                as.matrix(resid.gw(fit, type))
                              }
    }
    residuos <- resid.gw(object, type)
    residuos.sim <- cbind(as.matrix(rep(residuos, object$W)), 
                          residuos.sim)
    minimos <- apply(residuos.sim[, 2:rep], 1, min)
    maximos <- apply(residuos.sim[, 2:rep], 1, max)
    t <- 1:n
    normal.score <- qnorm(t/(n + 1))
    xx <- c(normal.score, rev(normal.score))
    yy <- c(minimos, rev(maximos))
    plot(normal.score, residuos, type = "l", xlab = "Standard normal quantiles", 
         ylab = paste("Residuals ", "(", type[1], ")", sep = ""), 
         main = title)
    polygon(xx, yy, col = "gray", border = NA)
    lines(normal.score, residuos)
  }
  et <- proc.time()
  if (trace > 0) 
    cat(paste("\nOverall envelope simulation process took", 
              round((et - st)[3], 2), "seconds"))
  else cat("\n")
  ans <- list(type = type, residuals = residuos, sim.residuals = residuos.sim)
  class(ans) <- "residuals.gw"
  return(ans)
}