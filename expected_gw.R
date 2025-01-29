#'Expected Probabilities and Frequencies for the GW Regression Model
#'
#'The \code{gw_expected} functions calculate the
#'probability distribution of the count response variable Y for each observation
#'and obtain the corresponding expected frequencies. It is an informal way of
#'assessing the fit of the GW regression model by comparing the predicted
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
#'@param object a fitted object of class inheriting from \code{"gw.fit"}.
#'
#'@return A list containing the following components:
#'
#'  \item{\code{frequencies}}{the expected counts for the GW fit.}
#'  \item{\code{observed_freq}}{the observed distribution.}
#'  \item{\code{probabilities}}{the expected distribution for the GW fit.}
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
#' ## Fit a GW regression model
#'
#' fit <- gw(formula = Favor~Lugar+Oficial+Resultado+CM, data = Germany)
#'
#' ## Compute the expected probabilities and the frequencies
#'
#' gw_expected(fit)
#' @export
expected.gw <- function(object) {

  source("dgw.R")

  mu <- object$fitted.values
  k <- object$betaIIpars[1]
  ro <- object$betaIIpars[2]
  a <- mu * (ro - 1)/k
  
  y <- object$Y

  row <- max(y) + 1
  col <- length(a)
  i <- 0:(row - 1)
  f0 <- matrix(0, col)
  prob <- matrix(0, nrow = row + 1, ncol = col)
  dist <- matrix(0, nrow = row + 1, ncol = col)
  for (j in seq_len(col)) {
    prob[1:row, j] <- dgw(i, a[j], k, ro)
    prob[row + 1, j] <- 1-sum(dgw(0:(row-1), a[j], k, ro))
    dist[, j] <- object$W[j] * prob[, j]
  }

  estimated_freq <- rowSums(dist)
  estimated_prob <- estimated_freq / sum(estimated_freq)

  observed_freq <- table(factor(y, levels = 0:row))
  observed_prob <- observed_freq / sum(observed_freq)

  dif <- sum(abs(observed_prob - estimated_prob))

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
