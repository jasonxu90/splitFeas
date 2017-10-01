#' Project onto a cube
#' 
#' \code{project_cube} computes the Euclidean projection of a point onto a cube.
#' 
#' @param x Point to project
#' @param center Center of the square
#' @param r Half the length of a side
#' @export
#' @examples
#' set.seed(12345)
#' p <- 3
#' center <- matrix(rnorm(p),p,1)
#' r <- runif(1)
#' x <- matrix(rnorm(p),p,1)
#' y <- project_cube(x,center,r)
project_cube = cmpfun(function(x, cnt, r) {
  y1 = min(max(cnt[1]-r, x[1]), cnt[1]+r)
  y23 = project_square(x[-1], cnt[-1], r)
  y = matrix(c(y1, y23), ncol=1)
  return(y)
})
