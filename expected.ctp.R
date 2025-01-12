#'Expected Probabilities and Frequencies for the CTP Model
#'
#'The \code{ctp_expected} functions calculate the
#'probability distribution of the count response variable Y for each observation
#'and obtain the corresponding expected frequencies. It is an informal way of
#'assessing the fit of the CTP model by comparing the predicted
#'distribution of counts with the observed distribution.
#'
#'The average expected probabilities are computed as \deqn{\bar(Pr)(y=k) =
#'\frac{1}{n} \sum_{i=1}^n \widehat{Pr}(y_i = k | x_i)}
#'
#'The expected frequencies are obtained by multiplying by n.
#'
#'Two measures are offered for summarizing the comparison between expected and
#'observed frequencies: the sum of the absolute value of differences and the sum
#'of the square of differences (similar to the Pearson statistic of goodness of
#'fit).
#'
#'@param object a fitted object of class inheriting from \code{"ctp.fit"}.
#'
#'@return A list containing the following components:
#'
#'  \item{\code{frequencies}}{the expected counts for the CTP fit.}
#'  \item{\code{observed_freq}}{the observed distribution.}
#'  \item{\code{probabilities}}{the expected distribution for the CTP fit.}
#'  \item{\code{dif}}{sum of the absolute value of differences between
#'  \code{frequencies} and \code{observed_freq}.} \item{\code{chi2}}{sum of the
#'  square of differences between \code{frequencies} and \code{observed_freq}.}
#'
#'@references
#'
#'J. M. Hilbe (2011). Negative Binomial Regression. (2nd ed.). Cambridge
#'University Press.
#'
#'M. Scott Long and Jeremy Freese (2014). Regression Models for Categorical
#'Dependent Variables using STATA. (3rd ed.). Stata Press.
#'
#'@name expected
NULL

#' @rdname expected
#' @examples
#' ## Fit a CTP model
#'
#' fit <- ctp.fit(formula = Favor~Lugar+Oficial+Resultado+CM, data = Alemania1)
#'
#' ## Compute the expected probabilities and the frequencies
#'
#' ctp_expected(fit)
#' @export
expected.ctp <- function(object) {

#    stopifnot(inherits(object, "glm_hP"))

  library(cpd)

  a <- object$parameters[1]
  b <- object$parameters[2]
  ro <- object$ro
  gama <- ro + 2*a

  y <- object$response

  # if (! is.null(object$response)) {
  #   y <- object$response
  # } else if (! is.null(object$model)) {
  #   y <- stats::model.response(object$model)
  # } else {
  #   formula <- stats::as.formula(object$call[["formula"]])
  #   ro.mu <- stats::model.frame(formula, data = object$data)
  #   y <- stats::model.response(ro.mu)
  # }

  row <- max(y) + 1
  col <- length(ro)
  i <- 0:(row - 1)
  f0 <- matrix(0, col)
  prob <- matrix(0, nrow = row + 1, ncol = col)
  dist <- matrix(0, nrow = row + 1, ncol = col)
  for (j in seq_len(col)) {
    prob[1:row, j] <- dctp(i, a, b, gama[j])
    #prob[row + 1, j] <- pctp(row - 1, a, b, gama[j], lower.tail = FALSE)
    prob[row + 1, j] <- 1-sum(dctp(0:(row-1), a, b, gama[j]))
    dist[, j] <- object$weights[j] * prob[, j]
  }

  estimated_freq <- rowSums(dist)
  estimated_prob <- estimated_freq / sum(estimated_freq)

  observed_freq <- table(factor(y, levels = 0:row))
  observed_prob <- observed_freq / sum(observed_freq)

#  dif <- sum(abs(observed_prob - estimated_prob))
  dif <- sum(abs(observed_freq - estimated_freq))
  
  chi2 <- sum((observed_freq[1:(row-1)] - estimated_freq[1:(row-1)]) ^ 2 / estimated_freq[1:(row-1)]) +
    (sum(observed_freq[row]) - sum(estimated_freq[row:(row+1)])) ^ 2 / sum(estimated_freq[row:(row+1)])
    
  barplot(observed_freq, xlab = "Y", ylab = "Frequencies")
  lines(estimated_freq, x = c(1:(row+1)), col = "blue")
  legend("topright", legend=c(paste("Dif = ", round(dif, 3)), paste(expression(chi ^ 2), "=", round(chi2, 3))))


    list(
    frequencies = estimated_freq,
    observed_freq = observed_freq,
    probabilities = estimated_prob,
    dif = dif,
    chi2 = chi2
  )
}
