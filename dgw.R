#pmf of a UGW(a, k, ro)

dgw <- function(x, a, k, ro){
  
  if (mode(c(x, a, k, ro)) != "numeric") 
    stop("non-numeric argument to mathematical function")
  if (any(c(a, k, ro) <= 0)) 
    stop("a, k and rho must be greater than 0")
  
  lpmf <- lgamma(a + ro) + lgamma(k + ro) + lgamma(a + x) + lgamma(k + x) - lgamma(ro) - 
    lgamma(a) - lgamma(k) - lgamma(a + k + ro + x) - lgamma(x + 1)
    pmf <- exp(lpmf)
    return(pmf)

  }