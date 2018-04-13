#' Constructor for pframe_1d class
#'
#' \code{pframe_1d} is a constructor function that generates an instance of the S3 pframe_1d class.
#' By passing in a set of parameters (or using defaults), a pframe_1d object is generated
#' where filtering methods can be applied on. pframe_1d is a frameholder that creates empty variables, vectors and matrices;
#' some of which are expected, whereas others are optional.
#'
#' @param N Number of particles
#' @param tau Number of time points
#' @param y Full data vector
#' @param resample Boolean for resampling
#' @param smooth Boolean for smooth
#' @param dist_f Distribution for project method (gaussian,t, gamma, etc)
#' @param params_f list of parameters for project method
#' @param dist_g Likelihood for update method (gaussian,t, gamma, etc)
#' @param params_g list of parameters for update method
#'
#' @return An object instantiated from pframe_1d class. Methods include init ,project, update, resample and
#' metric.
#'
#' @import ArgumentCheck
#'
#' @seealso \code{\link{project}}, \code{\link{update}} and \code{\link{resample}}
#'
#' @note \code{y} is the full data vector, therefore, we assume that filtering is not essential. Therefore,
#' it is preferred to run smoothing if \code{y} is available.
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#'
#' @examples
#' obj1<-pframe_1d()

#' @export
pframe_1d <- function(N=100,y=c(1,1,1),resample=TRUE,smooth=TRUE,sequential=TRUE,
                       dist_f="gaussian", params_f=list(),dist_g="gaussian",
                       params_g=list(sd=1),resamp_type="multinomial",...) {

  #assign variables
  tau<-length(y) #this may not be the case
  x0<-NULL
  t<-0

  #placeholder variables
  W<-matrix(NA,ncol=N,nrow=tau) #float
  X<-matrix(NA,ncol=N,nrow=tau+1) #float
  X_hist<-matrix(NA,ncol=N,nrow=tau+1) #float
  resamp_mat<-matrix(NA,ncol=N,nrow=tau) #float

  # t=0   x0      X[1, ]
  # t=1   x1  y1    X[2, ]   W[1,]
  # t=2   x2  y2    X[3, ]   W[2, ]
  # t=3   x3  y3    X[4, ]   W[3, ]

  #metrics
  m<-rep(NA,tau+1) #float
  conf<-matrix(rep(NA,2*(tau+1)),ncol=2) #float
  N_eff<-rep(NA,tau) #float
  sd_x<-rep(NA,tau+1)
  mse<-rep(NA,tau)


  #decide what is f and what is g
  #choose distributions
  f<-switch(dist_f,
            gaussian = "rnorm",
            t = "rt",
            gamma = "rgamma"
  )
  g<-switch(dist_g,
            gaussian = "dnorm",
            t = "dt",
            gamma = "dgamma"
  )

  #assigning values to object
  obj<- list(

    #variables
    x0=x0,
    y=y,
    tau=tau,
    t=t,
    N=N,
    X=X,
    X_hist=X_hist,
    W=W,
    resamp_mat=resamp_mat,

    #models
    project=list(dist_f=dist_f,f=f,params_f=params_f),
    update=list(dist_g=dist_g,g=g,params_g=params_g),
    prior=list(dist_x0=NULL,f_x0=NULL,params_x0=NULL), #init_x0 initialises
    resamp=if(resample) list(resamp_type=resamp_type,smooth=smooth),

    #metrics
    m=m,
    conf=conf,
    N_eff=N_eff,
    sd_x=sd_x,
    mse=mse

  )

  #assign attributes (or meta classes)
  attr(obj, "Dim")<-1
  attr(obj, "class") <- "pframe_1d"
  attr(obj, "resample")<- ifelse(resample,TRUE, FALSE)
  attr(obj, "smooth")<- ifelse(smooth,TRUE, FALSE)

  #output
  return(obj)
}




