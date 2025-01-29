#' Extract and Visualize GW Regression Model Residuals
#'

residuals.gw <- function (object, type = c("pearson", "response", "quantile"), rep = 19, envelope = FALSE, 
          title = "Simulated Envelope of Residuals", trace = FALSE, 
          parallel = TRUE, ncores = 2, ...) {
  
  type <- match.arg(type)
  source("pgw.R")
  
  resid.gw <- function(object, type) {
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
      return(qr)
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
                                varResponse <- getResponse(object$formula)
                                datos <- object$data[rep(1:nrow(object$data), 
                                                         object$W), ]
                                datos[varResponse] <- as.matrix(rgw(n, a, k, 
                                                                    ro))
                                while (!converged) {
                                  fit <- try(GWRM::gw(object$formula, data = datos, 
                                                      k = object$k), silent = TRUE)
                                  if (fit$aic > 0 && fit$betaIIpars[2] > 2) 
                                    converged <- TRUE
                                  else {
                                    datos[varResponse] <- as.matrix(rgw(n, a, 
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
                                varResponse <- getResponse(object$formula)
                                datos <- object$data[rep(1:nrow(object$data), 
                                                         object$W), ]
                                datos[varResponse] <- as.matrix(rgw(n, a, k, 
                                                                    ro))
                                while (!converged) {
                                  fit <- try(GWRM::gw(object$formula, data = datos, 
                                                      k = object$k), silent = TRUE)
                                  if (fit$aic > 0 && fit$betaIIpars[2] > 2) 
                                    converged <- TRUE
                                  else {
                                    datos[varResponse] <- as.matrix(rgw(n, a, 
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