#' Compute the approximate Hessian of the majorization.
#' 
#' \code{ddg} computes the Hessian of the majorization of the proximity function.
#' 
#' @param x non-anchor point
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param hgrad Handle for output mapping Jacobian
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10
#' p <- 2
#' x <- matrix(rnorm(p),p,1)
#' v <- 1
#' w <- 1
#' A <- matrix(rnorm(n*p),n,p)
#' hgrad <- function(x) {return(t(A))}
#' sol <- ddg(x,v,w,hgrad)
ddg = cmpfun(function(x,v,w,hgrad) {
  p <- length(x)
  dh <- as.matrix(hgrad(x))
  return(diag(sum(v),p,p) + sum(w)*dh%*%t(dh))
})
