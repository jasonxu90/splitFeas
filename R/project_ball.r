#' Projection onto a ball
#' 
#' \code{project_ball} computes the Euclidean projection of a point onto a ball.
#' 
#' @param x Point to project
#' @param center Center of the sphere
#' @param r Radius of the sphere
#' @export
#' @examples
#' set.seed(12345)
#' p <- 3
#' cnt <- rnorm(p)
#' r <- runif(1)
#' x <- rnorm(p)
#' y <- project_ball(x,cnt,r)
project_ball = cmpfun(function(x, cnt, r) {
  d = x - cnt
  dl = norm(as.matrix(d),'f')
  if (dl <= r) {
    return(x)
  } else {
    d = d / dl
    return(cnt + r*d)
  }
})
