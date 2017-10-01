#' Compute the gradient of the majorization.
#' 
#' \code{dg} computes the gradient of the majorization of the proximity function.
#' 
#' @param x non-anchor point
#' @param xa Anchor point
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param h Function handle for output mapping
#' @param hgrad Handle for output mapping Jacobian
#' @export

dg = cmpfun(function(x,xa,v,w,plist1,plist2,h,hgrad) {
  p <- length(x)
  SPC <- double(p)
  Sum2 <- sum(w)*h(x)
  
  SPQ <- double(length(Sum2))
  n1 <- length(plist1)
  for (i in 1:n1) {
    PCi <- as.matrix(plist1[[i]](xa))
    SPC <- SPC + v[i]*PCi
  }
  Sum1 <- sum(v)*x - SPC
  n2 <- length(plist2)
  for (j in 1:n2) {
    PQj <- as.matrix(plist2[[j]](h(xa)))
    SPQ <- SPQ + w[j]*PQj
  }
  Sum2 <- Sum2 - SPQ
  return(Sum1 + hgrad(x)%*%Sum2)
#  for (i in 1:n) {
#    z <- as.matrix(plist1[[i]](xa))
#    u <- u + v[i]*(x-z)
#  }
#  n <- length(plist2)
#  for (i in 1:n) {
#    z <- hgrad(x)%*%(h(x) - plist2[[i]](h(xa)))
#    u <- u + w[i]*z
#  }
#  return(u)
})
