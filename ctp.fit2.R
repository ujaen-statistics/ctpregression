#Funcion de ajuste del modelo de regresion basado en la CTP con a y b reales fijos y
#ro = gamma - 2 * a > 1 dependiendo de las covariables a traves de la media
#
#Se utiliza la funcion nloptr de la libreria nloptr
#--------------------------------------------------

ctp.fit2 <- function(formula,w=NULL,astart=NULL,bstart=NULL,betastart=NULL,iters=10000,data){

  library(hypergeo)
  library(gsl)
  library(nloptr)
  library(pracma)
  library(cpd)

  #Design matrix

  ro.mu <- model.frame(formula, data=data)
  y <- model.extract(ro.mu, "response")
  offset <- model.extract(ro.mu, "offset")
  
  if (is.null(w)) w <- rep(1,length(y))
  
  if (is.null(offset)){
    offset<-rep(1,length(y))
	  covoffset<-FALSE}
	else covoffset<-TRUE

  if (terms(formula)[[3]]==1){
	  matrizmu<-matrix(1,c(length(y),1))
	  namescovars.mu<-c("(Intercept)")
		}
  else {
    matrizmu<-model.matrix(terms(formula),model.frame(terms(formula),data=data,na.action=NULL))
    namescovars.mu<-dimnames(matrizmu[0,])[[2]]
  }

  
  ncovars.mu<-ncol(matrizmu)
  
  
  if (is.null(betastart)) {
    parameters <- list(formula = formula, data = data,
                       family = "poisson")
    betastart <- do.call(stats::glm, parameters)$coefficients
  }

   if ((is.null(astart)) && (is.null(bstart))){
     ctp <- cpd::fitctp(y)
     astart <- ctp$coefficients[1]
     bstart <- ctp$coefficients[2]
   }

  #if (is.null(betastart)) betastart <- rep(0,ncovars.mu)
  
  # if (is.null(astart)) astart <- 0.1
  # 
  # if (is.null(bstart)) bstart <- 0.1
  

  #Unidad imaginaria

  i<-sqrt(as.complex(-1))
  
  #Log-likelihood

  loglik <- function(p){
    beta <- p[1:(ncovars.mu)]
    a <- p[ncovars.mu + 1]
    b <- p[ncovars.mu + 2]
    mu <- offset*exp(matrizmu%*%beta)
    ro <- (a ^ 2 + b ^ 2) / mu + 1
   
    -sum(w*(2*Re(complex_gamma(ro+a+b*i,log=TRUE)+complex_gamma(a+b*i+y,log=TRUE)-complex_gamma(a+b*i,log=TRUE))-lgamma(ro)-lgamma(ro+2*a+y))) 
    #-sum(w*(2*Re(mapply(lngamma_complex,ro+a,b)+mapply(lngamma_complex,a+y,b)-mapply(lngamma_complex,a,b))-lgamma(ro)-lgamma(ro+2*a+y))) 
    }
  
  # Gradient of the log-likelihood
  
  loglik_grad <- function(p){
    beta <- p[1:(ncovars.mu)]
    a <- p[(ncovars.mu+1)]
    b <- p[(ncovars.mu+2)]
    mu <- offset*exp(matrizmu%*%beta)
    ro <- as.vector((a ^ 2 + b ^ 2) / mu + 1)
    grad_a <- - sum(w*(2*Re(mapply(pracma::psi,a+b*i+y)-pracma::psi(a+b*i)+(2*a/mu+1)*mapply(pracma::psi,ro+a+b*i))-(2*a/mu+2)*mapply(pracma::psi,2*a+ro+y)-2*a/mu*mapply(pracma::psi,ro)))
    grad_b <- - sum(w*(-2*Im(mapply(pracma::psi,a+b*i+y)-pracma::psi(a+b*i)+mapply(pracma::psi,ro+a+b*i))+2*b/mu*(2*Re(mapply(pracma::psi,ro+a+b*i))-mapply(pracma::psi,2*a+ro+y)-mapply(pracma::psi,ro))))
    grad_betas <- - colSums(w*((ro-1)*as.vector(-2*Re(mapply(pracma::psi,ro+a+b*i))+mapply(pracma::psi,2*a+ro+y)+mapply(pracma::psi,ro))*matrizmu))
    
    c(grad_betas, grad_a, grad_b)
  }
  
  # Inequality constraint
  
  ineq <- function(p) {
    beta <- p[1:(ncovars.mu)]
    a <- p[(ncovars.mu+1)]
    b <- p[(ncovars.mu+2)]
    mu <- offset*exp(matrizmu%*%beta)
     - a ^ 2 - b ^ 2 - (2 * a + 1) * mu
  }
  
  # Gradient of the inequality constraint
  
  grad_ineq <- function(p){
    beta <- p[1:(ncovars.mu)]
    a <- p[(ncovars.mu+1)]
    b <- p[(ncovars.mu+2)]
    mu <- offset*exp(matrizmu%*%beta)
    grad_ineq_betas <- matrix(0, ncol = ncol(matrizmu), nrow = nrow(matrizmu))
    for(j in 1:ncol(matrizmu)) grad_ineq_betas[, j] <- - (2 * a + 1) * mu * matrizmu[, j]
    
    cbind(grad_ineq_betas, - 2 * (a + mu), matrix(- 2 * rep(b, nrow(matrizmu))))
    
  }
  
  
  # Initial values
  
  p0 <- c( betastart, astart, bstart )
  
  # Lower and upper bounds of control (en realidad aqui no son necesarias porque no hay)

  lb <- rep(-Inf, ncovars.mu + 2)
  ub <- rep(Inf, ncovars.mu + 2)
  
  # Options for the optimization process

  #local_opts <- list( "algorithm" = "NLOPT_LD_LBFGS",
  #                    "xtol_rel" = 1.0e-7 )

  local_opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                      "xtol_rel" = 1.0e-7 )
  
  opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 10000,
                "local_opts" = local_opts )

  fit <- nloptr( x0 = p0,
                 eval_f = loglik,
                 eval_grad_f = loglik_grad,
                 eval_g_ineq = ineq,
                 eval_jac_g_ineq = grad_ineq,
                 lb = lb,
                 ub = ub,
                 opts = opts )

  fit$pars <- fit$solution

  results <- list(nloptr = fit, 
                  offset = unname(stats::model.extract(ro.mu, "offset")),
                  response = y,
                  aic = 2 * (fit$objective + sum(w*lfactorial(y))) + (ncovars.mu + 2) * 2, 
                  bic = 2 * (fit$objective + sum(w*lfactorial(y))) + (ncovars.mu + 2) * log(sum(w)), 
                  loglik = fit$objective+sum(w*lfactorial(y)),
                  call = match.call(),
                  model = formula, 
                  df.residual = sum(w) - (ncovars.mu + 2),
                  df.null = sum(w) - 2, 
                  betas = fit$pars[1:(ncovars.mu)],
                  parameters = fit$pars[(ncovars.mu + 1):(ncovars.mu + 2)],
                  data = data, 
                  weights = stats::setNames(w, seq(w)),
                  code = fit$status)

  mu <- as.vector(offset*exp(matrizmu %*% fit$pars[1:ncovars.mu]))
  results$fitted.values <- mu
  names(results$fitted.values) <- seq(results$fitted.values)
  results$linear.predictors <- as.vector(log(offset) + matrizmu %*% 
                                           fit$pars[1:ncovars.mu])
  names(results$linear.predictors) <- seq(results$linear.predictors)
  results$coefficients <- fit$pars[1:ncovars.mu]
  names(results$coefficients) <- colnames(matrizmu)
  results$residuals <- stats::setNames(y - results$fitted.values, 
                                       seq(y))
  results$ro <- 1+(fit$pars[(ncovars.mu + 1)]^2+fit$pars[(ncovars.mu + 2)]^2)/mu
  results$gama <- results$ro + 2*fit$pars[(ncovars.mu + 1)]
#  results$y <- y
  results$matrix.mu <- matrizmu
  results$model.mu <- ro.mu
  class(results) <- "glm_ctp"
  results
   
  }