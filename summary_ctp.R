#'Summarizing CTP Fits
#'
#'These functions are all methods for class \code{"CTP"} or
#'\code{summary.ctp} objects.
#'
#'@param object an object of class \code{"CTP"}, usually, a result of a call
#'  to \code{ctp.fit}.
#'@param x an object of class \code{"summary.ctp"}, usually, a result of a
#'  call to \code{summary.ctp}.
#'@inheritParams stats::print.summary.ctp
#' @examples
#' ## Fit a CTP model
#'
#' fit <- ctp.fit(formula = Scored~Place+Type+Result+Rival, data = Germany)
#'
#' ## Obtain a summary of the fitted model
#'
#' summary(fit)
#'@name summary.ctp
NULL

#' @rdname summary.ctp
#' @export
summary.ctp <- function (object, ...) {
  # matrix of beta coefficients
  coef.p <- object$betas
  k <- length(coef.p)
  s.err_b  <- se.ctp(object)[1:k]
  if (is.null(s.err_b)) {
    coef.table <- matrix(coef.p, ncol = 1)
    rownames(coef.table) <- names(object$coefficients)
    colnames(coef.table) <- c("Estimate")
  } else {
    zvalue <- coef.p / s.err_b
    pvalue <- 2 * stats::pnorm(abs(zvalue), lower.tail = FALSE)
    coef.table <- cbind(coef.p, s.err_b, zvalue, pvalue)
    rownames(coef.table) <- names(object$coefficients)
    colnames(coef.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }
  # matrix of a, b parameters
  coef <- object$parameters
  s.err <- unname(se.ctp(object)[(k+1):(k+2)])
  if (is.null(s.err)) {
    coef.table.p <- matrix(coef, ncol = 1)
    rownames(coef.table.p) <- c("a","b")
    colnames(coef.table.p) <- c("Estimate")
  } else {
    zvalue <- coef / s.err
    pvalue <- 2 * stats::pnorm(abs(zvalue), lower.tail = FALSE)
    coef.table.p <- cbind(coef, s.err, zvalue, pvalue)
    dimnames(coef.table.p) <- list(c("a","b"),c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  }
  keep <- match(c("call", "aic"),
                names(object),
                0L
  )
  result <- object[keep]
  others <- list(
    coeff_betas = coef.table,
    coeff_par = coef.table.p
  )
  result <- c(result, others)
  class(result) <- "summary.ctp"
  result
}

vcov.ctp <- function(object) {
  library(hypergeo)
  library(pracma)

  formula <- object$model
  ro.mu <- stats::model.frame(formula, data = object$data)
  matrix.mu <- stats::model.matrix(attr(ro.mu, "terms"), ro.mu)

  # n <- nrow(matrix.mu)
  # betas <- object$betas
  y     <- object$response
  a     <- object$parameters[1]
  b     <- object$parameters[2]
  c     <- complex(real = a,imaginary = b)
  ro    <- as.vector(object$ro)
  mu    <- object$fitted.values
  w     <- object$weights

  H11 <- crossprod(w*(2*Re(mapply(psi,1,ro+c))-mapply(psi,1,2*a+ro+y)-mapply(psi,1,ro))*(ro-1)*matrix.mu,(ro-1)*matrix.mu)+
    crossprod(w*(2*Re(mapply(psi,ro+c))-mapply(psi,2*a+ro+y)-mapply(psi,ro))*matrix.mu,(ro-1)*matrix.mu)
  H12 <- crossprod(w*(-2*(2*a/mu + 1)*Re(mapply(psi,1,ro+c))+(2*a/mu + 2)*mapply(psi,1,2*a+ro+y)+(2*a/mu)*mapply(psi,1,ro)),(ro-1)*matrix.mu)+
    crossprod(w*((2*a/mu)*(-2*Re(mapply(psi,ro+c))+mapply(psi,2*a+ro+y)+mapply(psi,ro))),matrix.mu)
  H13 <- crossprod(w*(2*Im(mapply(psi,1,ro+c))-(2*b/mu)*(2*Re(mapply(psi,1,ro+c))-mapply(psi,1,2*a+ro+y)-mapply(psi,1,ro))),(ro-1)*matrix.mu)+
    crossprod(w*((2*b/mu)*(-2*Re(mapply(psi,ro+c))+mapply(psi,2*a+ro+y)+mapply(psi,ro))),matrix.mu)
  H22 <- sum(w*(2*Re(mapply(psi,1,c+y)-mapply(psi,1,c)+(2*a/mu + 1)^2*mapply(psi,1,ro+c)+(2/mu)*mapply(psi,ro+c)) - 
    (2*a/mu + 2)^2*mapply(psi,1,ro+2*a+y) - (2*a/mu)^2*mapply(psi,1,ro) - (2/mu)*(mapply(psi,ro+2*a+y) + mapply(psi,ro))))
  H23 <- sum(w*(-2*Im(mapply(psi,1,c+y)-mapply(psi,1,c)+(2*a/mu + 1)*mapply(psi,1,ro+c)) + 
    (2*b/mu)*(2*(2*a/mu + 1)*Re(mapply(psi,1,ro+c)) - (2*a/mu+2)*mapply(psi,1,ro+2*a+y) - (2*a/mu)*mapply(psi,1,ro))))
  H33 <- sum(w*(-2*Re(mapply(psi,1,c+y)-mapply(psi,1,c))-(2*b/mu)*(4*Im(mapply(psi,1,ro+c))+(2*b/mu)*(mapply(psi,1,ro+2*a+y)+mapply(psi,1,ro))) + 
    (2/mu)*(2*Re(mapply(psi,ro+c))-mapply(psi,2*a+ro+y)-mapply(psi,ro)) + ((2*b/mu)^2-1)*2*Re(mapply(psi,1,ro+c))))
  
  H <- as.matrix(cbind(rbind(H11,H12,H13),c(H12,H22,H23),c(H13,H23,H33)))
  rownames(H) <- c(names(object$coefficients),"a","b")
  colnames(H) <- c(names(object$coefficients),"a","b")
  tryCatch(solve(-H),
           error = function(e) NULL
  )
}

se.ctp <- function(object) {

  V <- vcov.ctp(object)
  tryCatch(sqrt(diag(V)),
           error = function(e) NULL
  )
}

#' @rdname summary.ctp
#' @export
print.summary.ctp <- function (x, digits = max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
  )
  cat("Mean model coefficients (with log link):\n")
  if (ncol(x$coeff_betas) > 1) {
    stats::printCoefmat(x$coeff_betas,
                        digits = digits,
                        signif.stars = signif.stars,
                        na.print = "NA",
                        ...
    )
  } else {
    print(x$coeff_betas)
  }
  cat("\nParameters:\n")
  if (ncol(x$coeff_par) > 1) {
    stats::printCoefmat(x$coeff_par,
                        digits = digits,
                        signif.stars = signif.stars,
                        na.print = "NA",
                        ...
    )
  } else {
    print(x$coeff_par)
  }
  cat("\n")
  cat("AIC:", format(signif(x$aic, digits)), "\n")
  invisible(x)
}


#' #' @export
#' vcov.ctp <- function(object, ...) {
#'   if (! is.null(object$matrix.mu)) {
#'     matrix.mu <- object$matrix.mu
#'   } else if (! is.null(object$model.mu)) {
#'     matrix.mu <- stats::model.matrix(attr(object$model.mu, "terms"),
#'                                      object$model.mu)
#'   } else {
#'     formula.mu <- stats::as.formula(object$call[["formula.mu"]])
#'     a.mu <- stats::model.frame(formula.mu, data = object$data)
#'     matrix.mu <- stats::model.matrix(attr(a.mu, "terms"), a.mu)
#'   }
#' 
#'   n <- nrow(matrix.mu)
#'   betas   <- object$betas
#'   lambdas <- object$lambdas
#'   gammas  <- object$gammas
#'   offset  <- ifelse(is.null(object$offset), rep.int(0, n), object$offset)
#'   mu      <- exp(offset + matrix.mu %*% betas)
#' 
#'   In_beta  <- t(matrix.mu) %*% diag(as.vector((mu^2) /
#'               variances_hp(lambdas, gammas, object$maxiter_series, object$tol))) %*% matrix.mu
#'   solve(In_beta)
#' }
