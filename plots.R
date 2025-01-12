#' Plot Diagnostics for CTP Objects
#'
#' Two plots are currently available: a plot of residuals against fitted values
#' and a Normal Q-Q plot.
#'
#' @param x \code{CTP} object, typically the result of \code{\link{ctp.fit}}.
#' @param type the type of residuals. The alternatives are: "quantile"
#'   (default), "pearson" and "response". Can be abbreviated.
#' @param ask logical; if TRUE, the user is asked before each plot, see
#'   \code{\link[graphics]{par}}(ask=.).
#' @param ... other parameters to be passed through to plotting functions.
#' @name plots
NULL

#' @rdname plots
#' @examples
#' ## Fit the CTP model
#' fit <- ctp.fit(formula = Favor~Lugar+Oficial+Resultado+CM, data = Alemania1)
#' oldpar <- par(mfrow = c(1, 2))
#'
#' ## Plot diagnostics
#'
#' plot.ctp(fit)
#' par(oldpar)
#' @export
plot.ctp <- function(x, type = c("quantile", "pearson", "response"),
                       ask = prod(graphics::par("mfcol")) < 2 &&
                             grDevices::dev.interactive(),
                       ...) {
  type <- match.arg(type)
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  grDevices::dev.hold()
  graphics::plot(x$fitted.values,
                 residuals.ctp(x, type = type),
                 xlab = "Fitted values",
                 ylab = "Residuals",
                 main = "Residuals vs Fitted"
  )
  grDevices::dev.flush()
  stats::qqnorm(residuals.ctp(x, type = type))
  stats::qqline(residuals.ctp(x, type = type))
  invisible()
}
