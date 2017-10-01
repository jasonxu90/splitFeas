#' Backtracking Line Search
#' 
#' Given a descent direction \code{backtrack} computes a step size that ensures sufficient decrease in an objective.
#' 
#' @param x Current iterate
#' @param dx Descent direction
#' @param f objective function
#' @param df gradient of objective function
#' @param alpha sufficient decrease parameter
#' @param beta sufficient decrease parameter
#' @export
backtrack = cmpfun(function(x,dx,f,df,alpha=0.01,beta=0.8) {
  t <- 1
  g <- df(x)
  u <- alpha*sum(g*dx)
  k <- 1
  repeat {
    if (f(x + t*dx) <= f(x) + t*u) break
    t <- beta*t
    print(paste0("backtrack ",k))
    k <- k + 1
  }
  return(t)
})