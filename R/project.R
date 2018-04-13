#' @title Project particles
#'
#' @description Method to project particles forward using state model f(xn+1|xn). Given a set of particles and
#' state model, predicts where particles may be.
#'
#' @param obj Object of pframe_1d class
#'
#' @return Object of pframe_1d class
#'
#' @details In genetic algorithms, these functions are isomorphic to mutation/ prediction
#'
#' @seealso \code{\link{particle}}
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#'
#' @examples
#' project(init_x0(pframe_1d()))

#'@rdname project
#' @export
project<-function(obj,...){
  UseMethod("project",obj)
}

#'@rdname project
#' @export
project.pframe_1d <- function(obj,dist_f,params_f,...) {

  #check error
  if((obj$t)>obj$tau+1) stop("Cannot project over number of timesteps")
  #can extend to more timesteps for real time

  #add new project model
  if(!missing(dist_f) && !missing(params_f)){
    message("dist_f, params_f, (f) has changed for object")
    obj$project$dist_f<-dist_f
    obj$project$params_f<-params_f
    obj$project$f<-switch(dist_f,
                           gaussian = "rnorm",
                           t = "rt",
                           gamma = "rgamma"
    )
  }

  #project
  obj$X[obj$t+2,]<- obj$X_hist[obj$t+2,]<-obj$X[obj$t+1,]+do.call(obj$project$f,c(list(n=obj$N),obj$project$params_f))
  obj$t<-obj$t+1 #increase timestep

  return(obj)
}

#'@rdname project
#' @export
project.default <- function(obj) {

  stop("Object not of class pframe_1d")

}





