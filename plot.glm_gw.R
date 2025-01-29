plot.gw <- function(x, type = c("quantile", "pearson", "response"),
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
                   residuals.gw(x, type = type),
                   xlab = "Fitted values",
                   ylab = "Residuals",
                   main = "Residuals vs Fitted"
    )
    grDevices::dev.flush()
    stats::qqnorm(residuals.gw(x, type = type))
    stats::qqline(residuals.gw(x, type = type))
    invisible()
  }
