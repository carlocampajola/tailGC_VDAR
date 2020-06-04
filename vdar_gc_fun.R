#### VDAR GRANGER CAUSALITY TEST FUNCTIONS

#### test function for bivariate time series ####
## Z = Tx2 numeric array, p = order of BiDAR(p) model if specified, conf_level = confidence level to reject null, method = "LL" or "YW" to fit BiDAR (YW not implemented yet)
GCtailLR <- function(Z, p = 1, conf_level = 0.05, method = "LL"){
  if(dim(Z)[2] != 2 | is.null(dim(Z))) stop("ERROR in GCtailLR: Not a bivariate time series")
  x <- Z[,1]
  y <- Z[,2]
  if(method == "LL"){
    BiDAR <- estimateBiDARpLL(Z,p,...)
    outXy <- BiDAR[[1]]
    outYx <- BiDAR[[2]]
    outDARx <- estimateDARpLL(x,p,...)
    outDARy <- estimateDARpLL(y,p,...)
  } else if(method == "YW"){ ### STILL NOT IMPLEMENTED
    BiDAR <- estimateBiDARpYW(Z,p)
    outXy <- BiDAR[[1]]
    outYx <- BiDAR[[2]]
    outDARx <- estimateDARpYW(x,p)
    outDARy <- estimateDARpYW(y,p)
  } else {
    stop("ERROR in GCtailLR: method argument not valid")
  }
  logLxcondy <- outXy$LogLx
  logLycondx <- outYx$LogLx
  logLx0 <- outDARx$logL
  logLy0 <- outDARy$logL
  
  statLRtestX <- - 2*(logLx0 - logLxcondy)
  statLRtestY <- - 2*(logLy0 - logLycondx)
  pyTGCx <- pchisq(statLRtestX, p, lower.tail = F)
  pxTGCy <- pchisq(statLRtestY, p, lower.tail = F)
  
  yTGCx <- ifelse(pyTGCx < conf_level, T, F)
  xTGCy <- ifelse(pxTGCy < conf_level, T, F)
  
  output <- list(YtoX = list(stat = statLRtestX, p.value = pyTGCx, result = yTGCx), 
                 XtoY = list(stat = statLRtestY, p.value = pxTGCy, result = xTGCy),
                 call = list(Z=Z, p=p, conf_level=conf_level, method=method))
  
  return(output)
}

#### estimation function for BiDAR(p) with max likelihood estimation ####
estimateBiDARpLL <- function(Z, p, nu0 = rep(0.5,2), lambda0 = rep(0.5,2), gamma0 = NULL, chi0 = rep(0.5,2), 
                             lr=0.01, maxiter=100, verbose=T, grad_tol=1, n_small_grad=10, deltaf_tol=0.0001){
  if(is.null(gamma0)) gamma <- array(1/p, dim=c(2,2,p-1)) else gamma <- gamma0
  if(p==1) gamma <- c(0.5,0.5)
  nu <- nu0
  lambda <- lambda0
  chi <- chi0
  Zt <- list()
  Ztmp <- list()
  for(t in (p+1):nrow(Z)){
    Zt[[t-p]] <- Z[t,]
    Ztmp[[t-p]] <- Z[t - p:1,]
  }
  parvec <- c(nu, lambda, chi, as.vector(gamma))
  opt <- gdesc_Nesterov(par=parvec, fn=function(par, Zt, Ztmp, p) -loglik_BiDARp(Zt, Ztmp, p, par), 
                        gr=function(par, Zt, Ztmp, p) grad_BiDARp(Zt, Ztmp, p, par),
                        Zt=Zt, Ztmp=Ztmp, p=p,
                        lower=0.0001, upper=0.9999, lr = lr, maxiter = maxiter, verbose=verbose,
                        grad_tol = grad_tol, n_small_grad=n_small_grad, deltaf_tol=deltaf_tol)
  optpars <- opt$par
  nuopt <- optpars[1:2]
  lambdaopt <- optpars[3:4]
  chiopt <- optpars[5:6]
  gammaopt <- array(optpars[7:length(optpars)], dim=c(2,2,p-1))
  
  logLx <- loglik_BiDARp(Zt, Ztmp, p, optpars, submod="X")
  logLy <- loglik_BiDARp(Zt, Ztmp, p, optpars, submod="Y")
  
  output <- list(nu = nuopt, lambda = lambdaopt, chi = chiopt, gamma = gammaopt,
                 logLx=logLx, logLy=logLy)
}

#### loglikelihood term t of BiDAR(p) ####
## Zt = vector of dimension 2 with current observation, Ztmp = p x 2 matrix with past observations, p = order of model, 
loglikt_BiDARp <- function(Zt, Ztmp, p, nu, lambda, gamma, chi, submod = "full"){
  if(p > 1){
    if(any(dim(Ztmp) != c(p,2))) stop(paste0("ERROR: incorrect dimension of Ztmp argument \n Expected (", p, ",2), got ", dim(Ztmp), " instead"))
    if(length(nu) != 2 | length(lambda) != 2 | (any(dim(gamma) != c(2,2,p-1)) & p>1) | length(chi) != 2) stop("ERROR: incorrect dimension of model parameters")
    
    out <- 0
    
    if(submod != "Y"){
      selfcopy1 <- (1 - lambda[1])*drop(sum(gamma[1,1,]*1*(Zt[1] == Ztmp[p:2,1])) + (1 - sum(gamma[1,1,]))*1*(Zt[1] == Ztmp[1,1]))
      crosscopy1 <- lambda[1]*drop(sum(gamma[1,2,]*1*(Zt[1] == Ztmp[p:2,2])) + (1 - sum(gamma[1,2,]))*1*(Zt[1] == Ztmp[1,2]))
      copyterm1 <- selfcopy1 + crosscopy1
      out <- out + log(drop(nu[1]*copyterm1 + (1 - nu[1])*chi[1]^Zt[1] * (1 - chi[1])^(1 - Zt[1])))
    }
    if(submod != "X"){
      selfcopy2 <- (1 - lambda[2])*drop(sum(gamma[2,2,]*1*(Zt[2] == Ztmp[p:2,2])) + (1 - sum(gamma[2,2,]))*1*(Zt[2] == Ztmp[1,2]))
      crosscopy2 <- lambda[2]*drop(sum(gamma[2,1,]*1*(Zt[2] == Ztmp[p:2,1])) + (1 - sum(gamma[2,1,]))*1*(Zt[2] == Ztmp[1,1]))
      copyterm2 <- selfcopy2 + crosscopy2
      out <- out + log(drop(nu[2]*copyterm2 + (1 - nu[2])*chi[2]^Zt[2] * (1 - chi[2])^(1 - Zt[2])))
    }
  } else {
    out <- 0
    
    if(submod != "Y"){
      selfcopy1 <- (1 - lambda[1])*1*(Zt[1] == Ztmp[1])
      crosscopy1 <- lambda[1]*1*(Zt[1] == Ztmp[2])
      copyterm1 <- selfcopy1 + crosscopy1
      out <- out + log(drop(nu[1]*copyterm1 + (1 - nu[1])*chi[1]^Zt[1] * (1 - chi[1])^(1 - Zt[1])))
    }
    if(submod != "X"){
      selfcopy2 <- (1 - lambda[2])*1*(Zt[2] == Ztmp[2])
      crosscopy2 <- lambda[2]*1*(Zt[2] == Ztmp[1])
      copyterm2 <- selfcopy2 + crosscopy2
      out <- out + log(drop(nu[2]*copyterm2 + (1 - nu[2])*chi[2]^Zt[2] * (1 - chi[2])^(1 - Zt[2])))
    }
  }
  return(out)
}

#### total loglikelihood (just a wrapper for the sum) ####
## Zt and Ztmp are lists of length t-p containing time t and times t-p:t-1 values of Z, respectively
loglik_BiDARp <- function(Zt, Ztmp, p, pars, submod = "full"){
  nu <- pars[1:2]
  lambda <- pars[3:4]
  chi <- pars[5:6]
  gamma <- array(pars[7:length(pars)], dim=c(2,2,p-1))
  out <- sum(mapply(loglikt_BiDARp, Zt, Ztmp, MoreArgs = list(p=p, nu=nu, lambda=lambda, gamma=gamma, chi=chi, submod = submod)))
  return(out)
}

#### gradient term t of BiDAR(p) ####
gradt_BiDARp <- function(Zt, Ztmp, p, nu, lambda, gamma, chi){
  if(p > 1){
    if(any(dim(Ztmp) != c(p,2))) stop(paste0("ERROR: incorrect dimension of Ztmp argument \n Expected (", p, ",2), got ", dim(Ztmp), " instead"))
    if(length(nu) != 2 | length(lambda) != 2 | (any(dim(gamma) != c(2,2,p-1)) & p>1) | length(chi) != 2) stop("ERROR: incorrect dimension of model parameters")
    
    selfcopy1 <- drop(sum(gamma[1,1,]*1*(Zt[1] == Ztmp[p:2,1])) + (1 - sum(gamma[1,1,]))*1*(Zt[1] == Ztmp[1,1]))
    crosscopy1 <- drop(sum(gamma[1,2,]*1*(Zt[1] == Ztmp[p:2,2])) + (1 - sum(gamma[1,2,]))*1*(Zt[1] == Ztmp[1,2]))
    selfcopy2 <- drop(sum(gamma[2,2,]*1*(Zt[2] == Ztmp[p:2,2])) + (1 - sum(gamma[2,2,]))*1*(Zt[2] == Ztmp[1,2]))
    crosscopy2 <- drop(sum(gamma[2,1,]*1*(Zt[2] == Ztmp[p:2,1])) + (1 - sum(gamma[2,1,]))*1*(Zt[2] == Ztmp[1,1]))
  } else {
    selfcopy1 <- 1*(Zt[1] == Ztmp[1])
    crosscopy1 <- 1*(Zt[1] == Ztmp[2])
    selfcopy2 <- 1*(Zt[2] == Ztmp[2])
    crosscopy2 <- 1*(Zt[2] == Ztmp[1])
  }
  
  A1 <- (1 - lambda[1])*selfcopy1
  A2 <- lambda[1]*crosscopy1
  A <- A1 + A2
  C1 <- (1 - lambda[2])*selfcopy2
  C2 <- lambda[2]*crosscopy2
  C <- C1 + C2
  
  B <- chi[1]^Zt[1] * (1 - chi[1])^(1 - Zt[1])
  D <- chi[2]^Zt[2] * (1 - chi[2])^(1 - Zt[2])
  den1 <- drop(nu[1]*A + (1 - nu[1])*B)
  den2 <- drop(nu[2]*C + (1 - nu[2])*D)
  
  dnu <- dlam <- dchi <- c(0,0)
  dnu[1] <- (A - B)/den1
  dnu[2] <- (C - D)/den2
  dlam[1] <- nu[1]*(crosscopy1 - selfcopy1)/den1
  dlam[2] <- nu[2]*(crosscopy2 - selfcopy2)/den2
  dchi[1] <- (1 - nu[1])*(2*Zt[1]-1)/den1
  dchi[2] <- (1 - nu[2])*(2*Zt[2]-1)/den2
  
  if(p > 1){
    dgamma <- array(0, dim=c(2,2,p-1))
    dgamma[1,1,] <- nu[1]*(1-lambda[1])*(1*(Zt[1] == Ztmp[p:2,1]) - 1*(Zt[1] == Ztmp[1,1]))/den1
    dgamma[1,2,] <- nu[1]*lambda[1]*(1*(Zt[1] == Ztmp[p:2,2]) - 1*(Zt[1] == Ztmp[1,2]))/den1
    dgamma[2,1,] <- nu[2]*lambda[2]*(1*(Zt[2] == Ztmp[p:2,1]) - 1*(Zt[2] == Ztmp[1,1]))/den2
    dgamma[2,2,] <- nu[2]*(1-lambda[2])*(1*(Zt[2] == Ztmp[p:2,2]) - 1*(Zt[2] == Ztmp[1,2]))/den2
  } else dgamma <- c(0,0)
  
  out <- c(dnu, dlam, dchi, as.vector(dgamma))
  
  return(out)
}

#### total gradient of BiDAR(p) ####
grad_BiDARp <- function(Zt, Ztmp, p, pars){
  nu <- pars[1:2]
  lambda <- pars[3:4]
  chi <- pars[5:6]
  gamma <- array(pars[7:length(pars)], dim=c(2,2,p-1))
  preout <- mapply(gradt_BiDARp, Zt, Ztmp, MoreArgs = list(p=p, nu=nu, lambda=lambda, gamma=gamma, chi=chi))
  out <- rowSums(preout)
  return(out)
}

#### estimation function for univariate DAR(p) with max likelihood ####
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
  
  logLx <- opt$value
  
  output <- list(nu = nuopt, chi = chiopt, gamma = gammaopt,
                 logLx=logLx)
  return(output)
}

#### loglikelihood term t of univariate DAR(p) ####
loglikt_DARp <- function(xt, xtmp, p, nu, gamma, chi){
  if(length(xtmp) != p) stop(paste0("ERROR: incorrect dimension of xtmp argument \n Expected ", p, ", got ", dim(xtmp), " instead"))
  if(length(nu) != 1 | (length(gamma) != p-1 & p>1) | length(chi) != 1) stop("ERROR: incorrect dimension of model parameters")
  
  copyterm <- sum(gamma*1*(xt == xtmp[p:2])) + (1-sum(gamma))*1*(xt==xtmp[1])
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
  
  copyterm <- sum(gamma*1*(xt == xtmp[p:2])) + (1-sum(gamma))*1*(xt==xtmp[1])
  nocopyterm <- (chi^xt) * ((1 - chi)^(1 - xt))
  den <- (nu*copyterm + (1-nu)*nocopyterm)
  dnu <- (copyterm - nocopyterm)/den
  dchi <- (1-nu)*(2*xt - 1)/den
  dgamma <- nu*(1*(xt == xtmp[p:2]) - 1*(xt==xtmp[1]))/den
  
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

#### sampling function for DAR(p) model ####
sample_darp <- function(x, p, nu, gamma, chi){
  nur <- rbinom(1, 1, nu)
  if(p == 1) gammar <- 1 else gammar <- drop(rmultinom(1, 1, c(gamma, 1-sum(gamma))))
  chir <- rbinom(1, 1, chi)
  
  nur * sum(gammar * x[p:1]) + (1 - nur) * chir
}

#### sampling function for BiDAR(p) model ####
sample_BiDARp <- function(Z, p, nu, lambda, gamma, chi){
  nur <- c(rbinom(1, 1, nu[1]), rbinom(1, 1, nu[2]))
  lambdar <- c(rbinom(1, 1, lambda[1]), rbinom(1, 1, lambda[2]))
  if(p==1) gammar <- c(1,1) else gammar <- array(0, dim=c(2,2,p))
  if(lambdar[1] == 1){
    if(p > 1) gammar[1,2,] <- drop(rmultinom(1,1,c(gamma[1,2,], 1-sum(gamma[1,2,]))))
  } else {
    if(p > 1) gammar[1,1,] <- drop(rmultinom(1,1,c(gamma[1,1,], 1-sum(gamma[1,1,]))))
  }
  if(lambdar[2] == 1){
    if(p > 1) gammar[2,1,] <- drop(rmultinom(1,1,c(gamma[2,1,], 1-sum(gamma[2,1,]))))
  } else {
    if(p > 1) gammar[2,2,] <- drop(rmultinom(1,1,c(gamma[2,2,], 1-sum(gamma[2,2,]))))
  }
  chir <- c(rbinom(1,1,chi[1]), rbinom(1,1,chi[2]))
  
  out <- c(0,0)
  if(p == 1){
    out[1] <- nur[1]*(lambdar[1]*Z[2] + (1 - lambdar[1])*Z[1]) + (1 - nur[1])*chir[1]
    out[2] <- nur[2]*(lambdar[2]*Z[1] + (1 - lambdar[2])*Z[2]) + (1 - nur[2])*chir[2]
  } else {
    out[1] <- nur[1]*(lambdar[1]*sum(gammar[1,2,]*Z[p:1,2]) + (1-lambdar[1])*sum(gammar[1,1,]*Z[p:1,1]) ) + (1 - nur[1]) * chir[1]
    out[2] <- nur[2]*(lambdar[2]*sum(gammar[2,1,]*Z[p:1,1]) + (1-lambdar[2])*sum(gammar[2,2,]*Z[p:1,2]) ) + (1 - nur[2]) * chir[2]
  }
  
  return(out)
}


#### Nesterov Grad Descent ####
gdesc_Nesterov <- function(par, fn, gr, ..., lower=NA, upper=NA, lr=0.01, maxiter=100, verbose=T, grad_tol=0.0005, n_small_grad=10, deltaf_tol =0.0005){
  alpha <- lr
  prevv <- prevx <- prevy <- x <- y <- v <- par
  value <- initvalue <- fn(par, ...)
  small_grad <- 0
  for(iter in 1:maxiter){
    theta <- 2/(1+iter)
    thetap <- 2/(2+iter)
    if(any(par < lower | par > upper)){
      par <- prevpar 
      alpha <- alpha/2
    } else {
      prevpar <- prevx <- par
      prevy <- y
      prevv <- v
    }
    grad <- gr(prevy, ...)
    
    prop <- drop(grad/sqrt(sum(grad^2)))
    v <- prevv + alpha*prop/theta
    x <- (1 - theta)*prevx + theta*v
    y <- (1 - thetap)*x + thetap*v
    
    par <- x
    prevalue <- value
    value <- fn(par, ...)
    
    deltaf <- 1 - value/prevalue
    
    if(verbose) print(paste(c("iteration", iter, "- par =", round(par,3), "- fn =", round(value,2), "- fn rel. var =", round(deltaf,5), "- gr_norm =", round(sqrt(sum(grad^2)),2)), collapse=" "))
    if(sqrt(sum(grad^2)) < grad_tol | deltaf < deltaf_tol) small_grad <- small_grad + 1 else small_grad <- 0
    if(small_grad > n_small_grad) break
  }
  percimprove <- 100*value/initvalue
  cat(paste("Converged after", iter, "iterations - fn value", round(value, 3), 
            "\n final loglikelihood is", round(percimprove,2), "% of starting value",
            "\n final gradient norm", round(sqrt(sum(grad^2)),2))[1])
  return(list(par=par, value=value))
}