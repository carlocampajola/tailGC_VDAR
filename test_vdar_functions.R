##### test vdar_gc_fun functions


x <- c(1,0,0)
nu0 <- 0.7
chi0 <- 0.8
Tmax <- 3000

p <- 3
gamma0 <- c(0.6,0.2)

for(t in 1:Tmax){
  samp <- sample_darp(x[length(x) + 1 - p:1], p, nu0, gamma0, chi0)
  x <- c(x, samp)
}

out <- estimateDARpLL(x, p, lr=0.0001, maxiter=10000)

#####
## sample BiDAR(p)

source("vdar_gc_fun.R")
p <- 3
Z <- matrix(rbinom(p*2, 1, 0.5), ncol=2)
nu0 <- c(0.7, 0.4)
chi0 <- c(0.5, 0.5)
Tmax <- 3000

lambda0 <- c(0.6,0.4)
selfgamma <- 1/p + seq(1/(2*p), -1/(2*p), length.out =  p-1)
crossgamma <- 1/(2*p) + seq(1/(3*p), -1/(3*p), length.out =  p-1)
gamma0 <- array(0, dim=c(2,2,p-1))

for(ord in 1:(p-1)){
  gamma0[1,1,ord] <- gamma0[2,2,ord] <- selfgamma[ord]
  gamma0[1,2,ord] <- gamma0[2,1,ord] <- crossgamma[ord]
}

for(t in 1:Tmax){
  samp <- sample_BiDARp(Z[nrow(Z) + 1 - p:1,], p, nu0, lambda0, gamma0, chi0)
  Z <- rbind(Z, unname(samp))
}

## estimate from random starting points

out <- estimateBiDARpLL(Z, p, lr=0.001, maxiter=1000)

GCtailLR(Z)
