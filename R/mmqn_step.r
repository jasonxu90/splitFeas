#' MM-quasi-Newton step
#' 
#' \code{mmqn_step} computes a single step.
#' 
#' @param x Current iterate
#' @param f objective function
#' @param df gradient of objective function
#' @param v weights for first set of constraints
#' @param w weights for second set of constraints
#' @param plist1 list of projection functions for first set of constraints; each takes a single point and returns its projection
#' @param plist2 list of projection functions for second set of constraints; each takes a single point and returns its projection
#' @param h Function handle for output mapping
#' @param hgrad Handle for output mapping Jacobian
#' @param woodbury Boolean: TRUE to use the Woodbury inversion formula
#' @export
mmqn_step = cmpfun(function(x,v,w,plist1,plist2,f,df,h,hgrad,woodbury=TRUE) {
  df <- dg(x,x,v,w,plist1,plist2,h,hgrad)
  if(woodbury==TRUE){
#    H <- wood_inv(x,v,w,hgrad)
#    dx <- -H %*% df
    dx <- -wood_inv_solve(x,v,w,hgrad,df)
  }else{
    H <- ddg(x,v,w,hgrad)
    dx <- -solve(H,df)
  }
  t <- backtrack(x,dx,f,df)
  return(x + t*dx)
})