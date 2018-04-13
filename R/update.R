#' @title Update weights
#'
#' @description Each particle is given a different weight based on how close each particle is to
#' the data point by using the likelihood.
#'
#' @param obj Object of pframe_1d class
#'
#' @return Object of pframe_1d class
#'
#' @details In genetic algorithms, these functions are part of selection/updating
#'
#' @note One must becareful of particle degeneracy. Occasionally, all weight is given to one particle
#' only. This usually occurs when the state model does not conform with the data.
#'
#' @seealso \code{\link{particle}}
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#'
#' @examples
#' update(project(init_x0(pframe_1d())))

#'@rdname update
#' @export
update<-function(obj,...){
  UseMethod("update",obj)
}

#'@rdname update
#' @export
update.pframe_1d <- function(obj,dist_g, params_g,...) {

  #normalize weights
  norm<-function(w){
    return(w/sum(w))
  }

  #add new update model
  if(!missing(dist_g) && !missing(params_g)){
    message("dist_g, params_g, (g) has changed for object")
    obj$update$dist_g<-dist_g
    obj$update$params_g<-params_g
    obj$update$g<-switch(dist_g,
                          gaussian = "dnorm",
                          t = "dt",
                          gamma = "dgamma"
    )
  }

  #update weights
  obj$W[obj$t,]<-norm(do.call(obj$update$g,c(list(x=obj$y[obj$t],mean=obj$X[obj$t+1,]),obj$project$params_g)))

  return(obj)
}

#'@rdname update
#' @export
update.default <- function(obj) {

  stop("Object not of class pframe_1d")

}



