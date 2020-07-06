#### estimation function for BiDAR(p) with max likelihood estimation ####
#' @export
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
