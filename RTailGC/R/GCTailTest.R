#' Granger Causality in Tail Test function for a multivariate time series of tail hits.
#'
#' Function to test the presence of Granger Causality in tail relations by means of the likelihood ratio test proposed in Mazzarisi et al. (2020).
#' @param Z A TxN multivariate time series of 0 and 1 indicating whether a tail event occurred for variable \code{i=1,...,N} at time \code{t=1,...,T}.
#' @param maxp Maximum order of the BiDAR(p) model to test Tail Granger Causality for. Defaults to 1.
#' @param conf_level Desired confidence level for the likelihood ratio test. Defaults to 0.05.
#' @param method Estimation method for the BiDAR(p) models. Defaults to "LL" standing for log-likelihood maximization (no other option in this version).
#' @return A list with elements
#' \describe{
#'	\item{GC_adjacency}{An adjacency matrix for the resulting Tail Granger Causality graph. A 1 in position \code{[i,j]} indicates a rejected TGC test from variable \code{i} to variable \code{j}.}
#'	\item{p.values}{A matrix of p-values for the tests in the same shape as the \code{GC_adjacency} argument.}
#' 	\item{order_mat}{A matrix reporting the order at which the Tail Granger Causality was found between the variables. This uses the AIC criterion for optimal BiDAR(p) selection.}
#' }
#' @export
GCTailTest <- function(Z, maxp=1, conf_level = 0.05, method = "LL"){
  N <- ncol(Z)
  GCtail_mat <- order_mat <- pval_mat <- matrix(0, N, N)
  
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      Zsel <- Z[,c(i,j)]
      order_mat[i,j] <- order_mat[j,i] <- AIC_GCtail(Zsel, maxp=maxp)$selection
      res <- GCtailLR(Zsel, order_mat[i,j], conf_level=conf_level, method=method)
      GCtail_mat[i,j] <- 1*res$XtoY$result
      GCtail_mat[j,i] <- 1*res$YtoX$result
      pval_mat[i,j] <- res$XtoY$p.value
      pval_mat[j,i] <- res$YtoX$p.value
    }
  }
  
  out <- list(GC_adjacency = GCtail_mat, p.values = pval_mat, order_mat = order_mat)
}

#### AIC criterion for BiDAR model selection ####
AIC_GCtail <- function(Z, maxp=5){
  aicvec <- rep(0, maxp)
  for(p in 1:maxp){
    BiDAR <- estimateBiDARpLL(Z,p)
    npar <- 6 + 4*(p-1)
    aicvec[p] <- 2*npar - 2*(BiDAR$logLx + BiDAR$logLy)
  }
  out <- list(selection=which.min(aicvec), criteria=aicvec)
  return(out)
}

#### test function for bivariate time series for given order ####
## Z = Tx2 numeric array, p = order of BiDAR(p) model if specified, conf_level = confidence level to reject null, method = "LL" or "YW" to fit BiDAR (YW not implemented yet)
GCtailLR <- function(Z, p = 1, conf_level = 0.05, method = "LL"){
  if(dim(Z)[2] != 2 | is.null(dim(Z))) stop("ERROR in GCtailLR: Not a bivariate time series")
  x <- Z[,1]
  y <- Z[,2]
  if(method == "LL"){
    BiDAR <- estimateBiDARpLL(Z,p)
    outDARx <- estimateDARpLL(x,p)
    outDARy <- estimateDARpLL(y,p)
  } else if(method == "YW"){ ### STILL NOT IMPLEMENTED
    BiDAR <- estimateBiDARpYW(Z,p)
    outXy <- extractModelRes(BiDAR, 1)
    outYx <- extractModelRes(BiDAR, 2)
    outDARx <- estimateDARpYW(x,p)
    outDARy <- estimateDARpYW(y,p)
  } else {
    stop("ERROR in GCtailLR: method argument not valid")
  }
  logLxcondy <- BiDAR$logLx
  logLycondx <- BiDAR$logLy
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
                 call = list(p=p, conf_level=conf_level, method=method),
                 Z=Z)
  
  return(output)
}
