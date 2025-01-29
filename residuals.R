#' Extract and Visualize cTP Model Residuals
#'
#' residuals is a method which extracts model residuals from a \code{"CTP"}
#' object, commonly returned by \code{\link{ctp.fit}}. Optionally, it produces a
#' half normal plot with a simulated envelope of the residuals.
#'
#' The response residuals (\eqn{r_i=y_i - \mu_i}{r[i]=y[i] - \mu[i]}), Pearson
#' residuals (\eqn{r^P_i = r_i/\sigma_i}{r[i]^P = r[i]/\sigma[i]}) or randomized
#' quantile residuals are computed. The randomized quantile residuals are
#' obtained computing the cumulative probabilities that the fitted model being
#' less than \emph{y} and less or equal than \emph{y}. A random value from a
#' uniform distribution between both probabilities is generated and the value of
#' the residual is the standard normal variate with the same cumulative
#' probability. Four replications of the quantile residuals are recommended
#' because of the random component (see Dunn and Smyth, 1996 for more details).
#'
#' The functions \code{\link{plot.ctp}} generate a residuals against fitted values
#' plot and a Normal Q-Q plot.
#'
#' The Normal Q-Q plot may show an unsatisfactory pattern of the Pearson
#' residuals of a fitted model: then we are led to think that the model is
#' incorrectly specified.

#'The half normal plot with simulated envelope indicates that, under the
#'distribution of the response variable, the model is fine when only a few
#'points fall off the envelope.
#'
#'@param object an object of class \code{"CTP"}, typically the result of a call
#' to \code{\link{ctp.fit}}.
#'@param type the type of residuals which should be returned. The alternatives
#'  are: "pearson" (default), "response" and "quantile". Can be abbreviated.
#'@param envelope a logical value indicating whether the envelope should be
#'  computed.
#'@param rep	number of replications for envelope construction. Default is 19,
#'  that is the smallest 95 percent band that can be built.
#'@param title	a string indicating the main title of the envelope.
#'@param ... further arguments passed to or from other methods.
#'
#'@return Residual values.
#'
#'@seealso \code{\link{plots}}
#'
#'@references Peter K. Dunn and Gordon K. Smyth (1996). "Randomized quantile
#'  residuals". Journal of Computational and Graphical Statistics, 5(3), pp.
#'  236-244.
#'
#'  A. C. Atkinson (1981). "Two graphical displays for outlying and influential
#'  observations in regression". Biometrika, 68(1), pp. 13–20.
#'
#'@name residuals
NULL

#' @rdname residuals
#' @examples
#' ## Estimate a CTP model
#' fit <- ctp.fit(formula = Scored~Place+Type+Result+Rival, data = Germany)
#'
#' ## Compute residuals
#'
#' r <- residuals.ctp(fit)
#' @export
residuals.ctp <- function(object, type = c("pearson", "response", "quantile"),
                            envelope = FALSE, rep = 19,
                            title = "Simulated Envelope of Residuals", ...) {
  stopifnot(is.logical(envelope))
  stopifnot(is.numeric(rep), length(rep) == 1, rep >= 1)
  stopifnot(is.character(title), length(title) == 1)

  type <- match.arg(type)

  if (! envelope)
    return(res_ctp(object, type))

  residuals_sim <- matrix(0, nrow = length(object$residuals), ncol = rep)
  pb <- utils::txtProgressBar(min = 1, max = rep, initial = 1, style = 3)
  for (x in 1:rep) {
    utils::setTxtProgressBar(pb, x)
    tmp <- simulation_ctp(object, type)
    residuals_sim[, x] <- tmp
  }

  residuals_sim <- apply(residuals_sim, 2, sort)
  minima <- apply(residuals_sim, 1, min)
  maxima <- apply(residuals_sim, 1, max)

  resi <- sort(res_ctp(object, type))
  n <- length(resi)
  t <- 1:n

  normal_score <- stats::qnorm(t / (n + 1))

  xx <- c(normal_score, rev(normal_score))
  yy <- c(minima, rev(maxima))

  type_r <- paste0(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)))
  graphics::plot(normal_score,
                 resi,
                 type = "l",
                 xlab = "Standard normal quantiles",
                 ylab = paste(type_r, "residuals"),
                 main = title)
  graphics::polygon(xx, yy, col = "gray", border = NA)
  graphics::lines(normal_score, resi)

  structure(
    list(
      type = type,
      residuals = resi,
      sim_residuals = residuals_sim
    ),
    class = "residuals.ctp"
  )
}

res_ctp <- function(object, type) {
  a <- object$parameters[1]
  b <- object$parameters[2]
  ro <- object$ro
  gama <- ro + 2*a

  if (type == "pearson") {
    mu       <- object$fitted.values
    variance <- mu * (mu+gama-1)/(ro-2)
    # variance <- variances_ctp(lambda, gamma, object$maxiter_series, object$tol)
    r        <- as.vector(object$residuals / sqrt(variance))
    names(r) <- seq(r)
    return(r)
  }
  if (type == "response")
    return(object$residuals)

  # type == "quantile"
  y <- object$response

  # if (! is.null(object$y)) {
  #   y <- object$y
  # } else if (! is.null(object$model.mu)) {
  #   y <- stats::model.response(object$model.mu)
  # } else {
  #   formula.mu <- stats::as.formula(object$call[["formula.mu"]])
  #   a.mu <- stats::model.frame(formula.mu, data = object$data)
  #   y <- stats::model.response(a.mu)
  # }

  minp <- matrix(0, length(y))
  maxp <- matrix(0, length(y))
  for(ind in (1:length(y))){
  minp[ind] <- ifelse(y[ind] == 0, 0, cpd::pctp(y[ind] - 1, a, b, gama[ind]))
  maxp[ind] <- cpd::pctp(y[ind], a, b, gama[ind])}
  u <- stats::runif(length(y), minp, maxp)
  qr <- stats::qnorm(u)
  return(qr)
}

simulation_ctp <- function(object, type) {
  repeat {
    a <- object$parameters[1]
    b <- object$parameters[2]
    ro <- object$ro
    gama <- ro + 2*a

    response_n <- get_response(object$model)
    data <- object$dataset
    data[response_n] <- sapply(seq_len(nrow(data)),
                               function(i) cpd::rctp(1, a, b, gama[i]))
    the_call <- object$call
    the_call[["data"]] = data
    fit <- tryCatch(eval(the_call),
                    error = function(e) TRUE,
                    warning = function(e) TRUE)
    if (is.logical(fit))
      next
    if(any(is.nan(res_ctp(fit, type))))
      next
    break
  }
  return(res_ctp(fit, type))
}

get_response <- function(formula) {
  tt <- stats::terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1]
  response <- attr(tt, "response")
  vars[response]
}
