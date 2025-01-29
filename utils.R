#Goodness of fit measures

MPB <- function(fit){
  return(sum(fit$weights * (fit$response - fit$fitted.values)) / sum(fit$weights))
}

MAD <- function(fit){
  return(sum(fit$weights * abs(fit$response - fit$fitted.values)) / sum(fit$weights))
}

MSPE <- function(fit){
  return(sum(fit$weights * (fit$response - fit$fitted.values)^2) / sum(fit$weights))
}

pearson<- function(fit){
  return(sum(residuals.ctp(fit)^2))
}