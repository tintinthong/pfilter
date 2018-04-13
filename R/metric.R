#' @title Compute metrics
#'
#' @description Compute metrics from values of object
#'
#' @param obj Object of pframe_1d class
#'
#' @return Object of pframe_1d class
#'
#' @seealso \code{\link{particle}}
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#' @examples
#'metric(resamp(update(project(init_x0(pframe_1d())))))

#'@rdname metric
#' @export
metric<-function(obj,...){
  UseMethod("metric",obj)
}

#'@rdname metric
#' @export
metric.pframe_1d <- function(obj,...) {


  #assign values
  x<-obj$X[obj$t+1,]
  w<-obj$W[obj$t,]

  #compute metrics
  obj$m[obj$t+1]<-sum(w*x)
  obj$N_eff[obj$t]<-1/sum(w^2)
  obj$conf[obj$t+1,]<-quantile(x,probs=c(0.025,0.975))
  obj$sd_x[obj$t+1]<-sd(x)
  #obj$mse[obj$t]<-(obj$x-obj$m[obj$t+1])^2 #has to be known states

  return(obj)
}

#'@rdname metric
#' @export
metric.default <- function(obj) {
  stop("Object not of class pframe_1d")
}



