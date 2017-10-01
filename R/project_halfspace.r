#' Projection onto a halfspace
#' 
#' \code{project_halfspace} computes the Euclidean projection of a point onto a closed half-space.
#' The function returns the projection onto the set
#       C = {y : t(a)%*%y <= b}
#' 
#' @param x Point to project
#' @param a is the normal vector
#' @param b is the threshold
#' @export
#' @examples
#' set.seed(12345)
#' p <- 3
#' a <- matrix(rnorm(p),p,1)
#' a <- a/norm(a,'f')
#' b <- runif(1)
#' x <- matrix(rnorm(p),p,1)
#' y <- project_halfspace(x,a,b)
project_halfspace = cmpfun(function(x, a, b) {
  x <- as.matrix(x)
  y <- x
  a <- as.matrix(a)
  if (t(a)%*%x > b) {
    y = x - a%*%((t(a)%*%x - b)/norm(a,'f')**2)
  }
  return(y)
})
