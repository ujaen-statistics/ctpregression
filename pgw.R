#Permite que a, k y ro sean vectores
pgw <- function (q, a, k, ro, lower.tail = TRUE) 
{
  if (!(is.double(q) || is.integer(q)) || !(is.double(a) || 
                                            is.integer(a)) || !(is.double(k) || is.integer(k))|| 
      !(is.double(ro) || is.integer(ro))) 
    stop("Non-numeric argument to mathematical function")
  maximum <- max(length(q), length(a), length(k), length(ro))
  q <- rep(q, length.out = maximum)
  a <- rep(a, length.out = maximum)
  k <- rep(k, length.out = maximum)
  ro <- rep(ro, length.out = maximum)
  prob <- numeric(length = maximum)
  for (ind in seq_along(q)) {
    if (q[ind] < 0) {
      prob[ind] <- 0
      next
    }
    qlower <- trunc(q[ind])
    x <- 0:qlower
    probs <- dgw(x, a[ind], k[ind], ro[ind])
    if (lower.tail) {
      prob[ind] <- sum(probs)
    }
    else {
      prob[ind] <- 1 - sum(probs)
    }
  }
  return(prob)
}