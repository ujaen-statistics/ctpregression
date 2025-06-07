simulation_ctp <- function(object, type) {
  repeat {
    a <- object$parameters[1]
    b <- object$parameters[2]
    ro <- object$ro
    gama <- ro + 2*a
    
    response_n <- get_response(object$model)
    #data <- object$dataset
    data <- object$data
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
