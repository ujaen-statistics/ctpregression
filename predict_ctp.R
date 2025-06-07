#' Predict Method for CTP Fits
#'
#' Obtains predictions from a fitted \code{"CTP"} object.
#'
#' @param object	a fitted object of class inheriting from \code{"CTP"}.
#' @param newdata	optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the fitted linear predictors are used.
#' @param type the type of prediction required. The default is on the scale of
#'   the linear predictors; the alternative \code{"response"} is on the scale of
#'   the response variable.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A vector with the prediction means.
#'
#' @examples
#' data(Germany)
#' ## Estimate a CTP model
#' fit <- ctp.fit(formula = Scored~Place+Type+Result+Rival, data = Germany)
#'
#' ## As the newdata parameter is not used the fitted values are obtained
#' predict(fit, type = "response")
#' @export
#' 
predict.CTP <- function (object = NULL, newdata = NULL, ...) 
{
  tt <- terms(object)
  if (!inherits(object, "glm_ctp")) 
    warning("calling predict.ctp(<fake-glm_ctp-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    nobs <- nrow(as.matrix(m))
    X <- model.matrix(Terms, m)
    offset <- rep(0, nobs)
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }
  ncovars <- ncol(X)
  beta <- object$coefficients[1:(ncovars)]
  #a <- object$parameters[1]
  #b <- object$parameters[2]
  
  if (is.null(object$offset)) 
    fits <- exp(X %*% beta)
  else fits <- exp(offset + X %*% beta)
  predictor <- cbind(fits)
  colnames(predictor) <- c("fit")
  ans <- data.frame(predictor)
  return(ans)
}