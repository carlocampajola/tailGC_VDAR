#### estimation function for univariate DAR(p) with max likelihood ####
#' @export
estimateDARpLL <- function(x, p, nu0=0.5, gamma0 = NULL, chi0 = 0.5, lr=0.01, maxiter=100, verbose=T, grad_tol=1, n_small_grad=10, deltaf_tol = 0.001){
  if(is.null(gamma0) & p > 1) gamma <- rep(1/p, p-1) else if(is.null(gamma0)) gamma <- 0.999 else gamma <- gamma0
  nu <- nu0
  chi <- chi0
  
  xt <- xtmp <- list()
  for(t in (p+1):length(x)){
    xt[[t-p]] <- x[t]
    xtmp[[t-p]] <- x[t - p:1]
  }
  parvec <- c(nu, chi, gamma)
  
  # opt <- optim(par=parvec, fn=function(par, xt, xtmp, p) -loglik_DARp(xt, xtmp, p, par),
  #              gr = function(par, xt, xtmp, p) -grad_DARp(xt, xtmp, p, par),
  #              xt=xt, xtmp=xtmp, p=p,
  #              method="L-BFGS-B", lower=0.0001, upper=0.9999,
  #                 control = list(trace = 3, maxit = 200))
  opt <- gdesc_Nesterov(par=parvec, fn=function(par, xt, xtmp, p) -loglik_DARp(xt, xtmp, p, par),
                        gr = function(par, xt, xtmp, p) grad_DARp(xt, xtmp, p, par),
                        xt = xt, xtmp = xtmp, p = p,
                        lower = 0.0001, upper = 0.9999, lr = lr, maxiter = maxiter, verbose=verbose,
                        grad_tol = grad_tol, n_small_grad=n_small_grad, deltaf_tol=deltaf_tol)
  optpars <- opt$par
  nuopt <- optpars[1]
  chiopt <- optpars[2]
  gammaopt <- optpars[3:length(optpars)]
  
  logLx <- -opt$value
  
  output <- list(nu = nuopt, chi = chiopt, gamma = gammaopt,
                 logLx=logLx)
  return(output)
}

#### loglikelihood term t of univariate DAR(p) ####
loglikt_DARp <- function(xt, xtmp, p, nu, gamma, chi){
  if(length(xtmp) != p) stop(paste0("ERROR: incorrect dimension of xtmp argument \n Expected ", p, ", got ", dim(xtmp), " instead"))
  if(length(nu) != 1 | (length(gamma) != p-1 & p>1) | length(chi) != 1) stop("ERROR: incorrect dimension of model parameters")
  
  if(p > 1){
    copyterm <- sum(gamma*1*(xt == xtmp[p:2])) + (1-sum(gamma))*1*(xt==xtmp[1])
  } else {
    copyterm <- 1*(xt == xtmp)
  }
  
  nocopyterm <- (chi^xt) * ((1 - chi)^(1 - xt))
  out <- log(nu*copyterm + (1-nu)*nocopyterm)
  
  return(out) 
}

#### total loglikelihood (wrapper for the sum) ####
## xt and xtmp are lists of length t-p containing time t and times t-p:t-1 values of x, respectively
loglik_DARp <- function(xt, xtmp, p, pars){
  nu <- pars[1]
  chi <- pars[2]
  gamma <- pars[3:length(pars)]
  out <- sum(mapply(loglikt_DARp, xt, xtmp, MoreArgs = list(p=p, nu=nu, gamma=gamma, chi=chi)))
  return(out)
}

#### gradient term t of univariate DAR(p) ####
gradt_DARp <- function(xt, xtmp, p, nu, gamma, chi){
  if(length(xtmp) != p) stop(paste0("ERROR: incorrect dimension of xtmp argument \n Expected ", p, ", got ", dim(xtmp), " instead"))
  if(length(nu) != 1 | (length(gamma) != p-1 & p>1) | length(chi) != 1) stop("ERROR: incorrect dimension of model parameters")
  
  if(p > 1){
    copyterm <- sum(gamma*1*(xt == xtmp[p:2])) + (1-sum(gamma))*1*(xt==xtmp[1])
  } else {
    copyterm <- 1*(xt == xtmp)
  }
  
  nocopyterm <- (chi^xt) * ((1 - chi)^(1 - xt))
  den <- (nu*copyterm + (1-nu)*nocopyterm)
  dnu <- (copyterm - nocopyterm)/den
  dchi <- (1-nu)*(2*xt - 1)/den
  if(p > 1){
    dgamma <- nu*(1*(xt == xtmp[p:2]) - 1*(xt==xtmp[1]))/den
  } else {
    dgamma <- 0
  }
  
  out <- c(dnu, dchi, dgamma)
  
  return(out)
}

#### total gradient (wrapper for the sum) of univariate DAR(p) ####
grad_DARp <- function(xt, xtmp, p, pars){
  nu <- pars[1]
  chi <- pars[2]
  gamma <- pars[3:length(pars)]
  preout <- mapply(gradt_DARp, xt, xtmp, MoreArgs = list(p=p, nu=nu, gamma=gamma, chi=chi))
  out <- rowSums(preout)
  return(out)
}
