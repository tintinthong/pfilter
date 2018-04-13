#' @title Show method
#'
#' @description Show object of pframe_1d class
#'
#' @param x is object, i.e. pframe_1d
#'
#' @return Description of the object you will be returning.
#' @details In genetic algorithms, these functions are isomorphic to mutation/ prediction
#' @seealso \code{\link{particle}} This part tells you to refer to other code

#' @note You can also leave an important note here
#' @references Links can be useful \url{http://en.wikipedia.org/wiki/Fermat's_little_theorem}
#' @author Justin Thong \email{justinthong93@gmail.com}
#' @examples
#' a+b

#'@rdname show
#' @export
show<-function(obj,...){
  UseMethod("show",obj)
}

#'@rdname show
#' @export
show.pframe_1d <- function(obj,dist_f=rnorm,...) {

  if((obj$t)>obj$tau+1) stop("Cannot show over number of timesteps")
  #can extend X with this

  f<-match.fun(dist_f)
  #obj$f$f<-f
  obj$X[obj$t+2,]<- obj$X_hist[obj$t+2,]<-obj$X[obj$t+1,]+f(obj$N,...) #t is 0 at x0
  obj$t<-obj$t+1


  #it's best to restrict to supported dist_f's such as norm,beta,...

  #assign parameters into object
  obj$f$f<-f
  obj$f$dist_f<-dist_f

  return(obj)
}

#'@rdname show
#' @export
show.default <- function(obj) {

  stop("Object not of class pframe_1d")

}

# t=0   x0      X[1, ]
# t=1   x1  y1    X[2, ]   W[1,]
# t=2   x2  y2    X[3, ]   W[2, ]
# t=3   x3  y3    X[4, ]   W[3, ]



