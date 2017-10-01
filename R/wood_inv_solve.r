#' Compute the inverse approximate Hessian of the majorization using the Woodbury inversion formula.
#
#' \code{wood_inv_solve} computes the inverse of the Hessian term of the majorization of the proximity function
#' using the Woodbury formula. The function \code{mmqn_step} invokes \code{wood_inv_solve} instead of {ddg} if the argument
#' \code{woodbury=TRUE}. This should be used when p << n.
#'
#' 
#' @param x non-anchor point
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param hgrad Handle for output mapping Jacobian
#' @param df Right hand side
#' @export
wood_inv_solve = cmpfun(function(x,v,w,hgrad,df) {
  n <- length(x)
  vv <- sum(v); ww <- sum(w)
  dh <- as.matrix(hgrad(x))
  p <- dim(dh)[2]  
  return( (1/vv)*df - (ww/vv^2)*dh %*% solve(diag(p) + (ww/vv)*t(dh) %*% dh , t(dh)%*%df ) )
})
