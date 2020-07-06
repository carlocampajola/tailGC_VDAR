#' sampling function for DAR(p) model
#'
#' This function samples an observation from a DAR(p) model starting from a set of initial configurations and given parameters.
#'
#' @param x A vector of length \code{p} containing the initial configurations.
#' @param p The order of the DAR(p) model
#' @param nu The value of the copy probability parameter \nu. Must be in (0,1).
#' @param gamma A vector of length \code{p-1} containing the probabilities of copying from the \code{p-1} past values of \code{x}, starting from the last element. Each must be in (0,1) and the sum has to be smaller than 1.
#' @param chi The value of the independent sampling probability \chi. Must be in (0,1).
#' @return A sample from the DAR(p) conditional probability distribution.
#' @export
sample_darp <- function(x, p, nu, gamma, chi){
  nur <- rbinom(1, 1, nu)
  if(p == 1) gammar <- 1 else gammar <- drop(rmultinom(1, 1, c(gamma, 1-sum(gamma))))
  chir <- rbinom(1, 1, chi)
  
  nur * sum(gammar * x[p:1]) + (1 - nur) * chir
}

#' sampling function for BiDAR(p) model
#'
#' This function samples an observation from a BiDAR(p) model starting from a set of initial configurations and given parameters.
#'
#' @param x A matrix with \code{p} rows and 2 columns containing the initial configurations.
#' @param p The order of the BiDAR(p) model
#' @param nu A vector of length 2 containing the values of the copy probability parameters \nu. Each must be in (0,1).
#' @param gamma An array of dimension \code{c(2,2,p-1)} containing the probabilities of copying from the \code{p-1} past values of \code{x}, starting from the last row. Each element \code{gamma[i,j,l]} is the probability that variable \code{i} copies from variable \code{j} at lag \code{l}. Each must be in (0,1) and the sum over the last dimension has to be smaller than 1.
#' @param chi A vector of length 2 containing the values of the independent sampling probability \chi. Each must be in (0,1).
#' @return A sample from the DAR(p) conditional probability distribution.
#' @export
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
