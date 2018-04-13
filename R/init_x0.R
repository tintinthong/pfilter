#' @title Initialise prior distributions for pframe_1d
#'
#' @description  \code{init_x0} is a method
#' to initialise the prior state of particles in \code{x0} in pframe_1d object.
#'
#' @param obj Object of pframe_1d class
#'
#' @return Updates or initialise pframe_1d with prior distribution.
#'
#' @seealso \code{\link{project}}, \code{\link{update}} and \code{\link{resample}}
#'
#' @note \code{init_x0} is not needed if all arguments are provided by the constructor pframe_1d.
#' This function is a convenience for user. A proficient user will be able to update or initialise variables by hand.
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#'
#' @examples
#' init_x0(pframe_1d())
#'
#'@rdname init_x0
#' @export
init_x0<-function(obj,...){
  UseMethod("init_x0",obj)
}
#'@rdname init_x0
#' @param obj Object of pframe_1d
#' @param dist_x0 Prior distribution (gaussian,t, gamma, etc)
#' @param ... optional arguments
#' @export
init_x0.pframe_1d <- function(obj,dist_x0="gaussian",
                           params_x0=list(mean=0,sd=1),...) {

  if(is.null(obj$x0)){

    #choose distributions
    obj$prior$f_x0<-switch(dist_x0,
           gaussian = "rnorm",
           t = "rt",
           gamma = "rgamma"
    )

    #assign object values
    obj$x0<-do.call(obj$prior$f_x0,c(list(n=obj$N),params_x0))
    obj$X[1,]<-obj$X_hist[1,]<-obj$x0
    obj$prior$dist_x0<-dist_x0
    obj$prior$params_x0<-params_x0

    obj$m[1]<-mean(obj$x0)
    obj$conf[1,]<-quantile(obj$x0,probs=c(0.025,0.975))
    obj$sd_x[1]<-sd(obj$x0)


  }

  return(obj)
}

#'@rdname init_x0
#' @export
init_x0.default <- function(obj) {
  stop("Object not of class pframe_1d")
}



