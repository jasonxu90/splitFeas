#' Project onto a square
#' 
#' \code{project_square} computes the Euclidean projection of a point onto a square.
#' 
#' @param x Point to project
#' @param center Center of the square
#' @param r Half the length of a side
#' @export
#' @examples
#' set.seed(12345)
#' p <- 2
#' center <- matrix(rnorm(p),p,1)
#' r <- runif(1)
#' x <- matrix(rnorm(p),p,1)
#' y <- project_square(x,center,r)
project_square = cmpfun(function(x, cnt, r) {
  d = length(cnt)
  if (norm(as.matrix(x-cnt),'I') <= r) y = x
  if (x[1] - cnt[1] > r & x[2] - cnt[2] > r) y = matrix(c(cnt[1] + r, cnt[2] + r), d, 1)
  if (x[1] - cnt[1] < -r & x[2] - cnt[2] > r) y = matrix(c(cnt[1] - r, cnt[2] + r), d, 1)
  if (x[1] - cnt[1] < -r & x[2] - cnt[2] < -r) y = matrix(c(cnt[1] -r ,cnt[2] - r), d, 1)
  if (x[1] - cnt[1] > r & x[2] -cnt[2] < -r) y = matrix(c(cnt[1] + r, cnt[2] - r), d, 1)
  if (abs(x[1] - cnt[1]) <= r & x[2] - cnt[2] > r) y = matrix(c(x[1], cnt[2]+r), d, 1)
  if (abs(x[1] - cnt[1]) <= r & x[2] - cnt[2] < -r) y = matrix(c(x[1], cnt[2]-r), d, 1)
  if (x[1] - cnt[1] > r & abs(x[2] - cnt[2]) <= r) y = matrix(c(cnt[1]+r, x[2]), d, 1)
  if (x[1] - cnt[1] < -r & abs(x[2] - cnt[2]) <= r) y = matrix(c(cnt[1]-r, x[2]), d, 1)
  return(y)
})
