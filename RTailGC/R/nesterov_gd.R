#' Nesterov Gradient Descent
#'
#' This is a generic implementation of the Nesterov accelerated gradient descent algorithm. It takes a vector of parameters par as input, together with a function fn to be minimized and its gradient gr. Optional arguments are the lower and upper bounds for the par values, the learning rate lr, maximum number of iterations, stopping criteria and whether verbose reporting should be enabled.
#'
#' @param par A numeric vector with the starting values of the parameters.
#' @param fn A function returning the value of the cost function to be minimized taking as first argument the vector of parameters in par.
#' @param gr A function returning the gradient of the cost function taking as first argument the vector of parameters in par.
#' @param ... Additional arguments for fn and gr.
#' @param lower A scalar indicating a lower bound for the parameters in par. NA implies no lower bound. Defaults to NA.
#' @param upper Same as lower, but for the upper bound. Defaults to NA.
#' @param lr The learning rate of the algorithm. Defaults to 0.01.
#' @param maxiter Maximum number of gradient descent steps. Defaults to 100.
#' @param verbose Whether progress should be reported with text messages. Defaults to TRUE.
#' @param grad_tol threshold for the modulus of the gradient below which the algorithm stops after n_small_grad iterations. Defaults to 0.0005.
#' @param n_small_grad Number of iterations after which to stop if gradient modulus falls below the grad_tol level. Defaults to 10.
#' @param deltaf_tol Threshold for the relative decrease of the cost function below which the algorithm stops. Defaults to 0.0005.
#' @return A list containing an argument \code{par} with the optimized values of the parameters and an argument \code{value} with the value of the objective function at the optimized coordinates.


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
            "\n final gradient norm", round(sqrt(sum(grad^2)),2),"\n")[1])
  return(list(par=par, value=value))
}
