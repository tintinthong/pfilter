#' @title Resample particles
#'
#' @description Resample particles based on computed weights.
#'
#' @param obj Object of pframe_1d class
#'
#' @return Object of pframe_1d class
#'
#' @details In genetic algorithms, these functions are part of selection/updating
#'
#' @seealso \code{\link{particle}},\code{\link{update}}

#' @note Since resampling is computationally costly, it be chosen to be done only "sometimes".
#' This is why the effective sample size is used, \code{N_eff}
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#'
#' @examples
#' resamp(update(project(init_x0(pframe_1d()))))

#'@rdname resamp
#' @export
resamp<-function(obj,...){
  UseMethod("resamp",obj)
}

#'@rdname resamp
#' @export
resamp.pframe_1d <- function(obj,...) {
  #if effective sample size below threshold then don't resamp
  #n_resamp=obj$N,type="multinomial",smooth=FALSE
  #check if resample is TRUE

  #resample  indexes
  ind<-sample(1:obj$N,obj$N, replace=TRUE, prob=obj$W[obj$t,])

  if(!attr(obj, "smooth")){
    obj$X[obj$t,]<-obj$X[obj$t,ind]
  }else{
    obj$X<-obj$X[,ind]
  }

  #assign values
  obj$resamp_mat[obj$t,]<-ind

  return(obj)
}

#'@rdname resamp
#' @export
resamp.default <- function(obj) {
  stop("Object not of class pframe_1d")
}



