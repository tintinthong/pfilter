#' @title Running Particle Filter
#'
#' @description Wrapper for project, update, resample to run through whole entire particle filter,
#' given the full \code{y} is know.
#'
#' @param obj Object of pframe_1d class
#'
#' @return Object of pframe_1d class
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#'
#' @examples
#' project(init_x0(pframe_1d()))
#'
#'@rdname particle
#' @export
particle<-function(obj,...){
  UseMethod("particle",obj)
}

#'@rdname particle
#' @import dplyr
#' @export
#pframe_1d x0=rnorm(20,0,1), y=c(1,1,1)
#project obj
#update obj dist_g
#obj,n_resamp=obj$N,type="multinomial",smooth=FALSE,
particle.pframe_1d <- function(obj,...) {

  #initialise x0
  obj<-obj %>% init_x0()

  i<-1
  while(i<=obj$tau){
    #must always project and update in that order. Resample is optional
    obj %>% project() %>%update() %>% resamp()%>% metric()->obj
    i<-i+1
  }

  #add conditional for resampling

  return(obj)

}

#'@rdname particle
#' @export
particle.default <- function(obj) {

  stop("Object not of class pframe_1d")

}




