

#load libraries
source("resample.R")
source("random-walk.R")

set.seed(123)
#setup
tau<-100 #no of states, excluding known state x0 
N<-500 #no of particles 
sigma<-1 # fixed sigma for process model
sigma.meas<-1 # fixed sigma.meas for measurement model
N.thr<-as.integer(0.5*N)
d<-2 #no of dimensions

#generate true state and data
x0<-0
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.y=sigma.meas)

#matrix of true states
x<-matrix(NA,d*tau,nrow=d,ncol=tau) # matrix of true states(rows is the dimension, cols is time)



 
#-----------------------OLD-----------------

  
  x[1,]<-cumsum(rnorm(tau)) #propagate true state using independent process model
  x[2,]<-cumsum(rnorm(tau)) # propagate true state using independent process model
  
  y<-x+rnorm(d*tau) #generate measurements from true states

#-----------------------NEW-----------------
source("resample.R")
source("random-walk.R")
#generate state

X<-rand.walk() #100 time states including x0
y<-rand.y(X) #generate measurement(measurement not meant to be designed for state x0)


#make particle filter 

x.pf<-array(NA, c(d,tau+1,N)) #(dimension, time state, particles)
w.pf<-array(NA,c(d,tau+1,N))  #(dimension, time state, particles)
 
N.eff<-matrix(NA,nrow=d,ncol=tau) #(dimension, effective sample size )
w.var<-matrix(NA,nrow=d,ncol=tau) #(dimension, weights )
  
#x.pf<-matrix(rep(NA,N*(tau+1)),nrow=tau+1)
#w.pf<-matrix(rep(NA,N*(tau+1)),nrow=tau+1)
#w.var<-array(NA,c(d,tau)) #w.var<-rep(NA,tau+1)
  
#initialising first time step
x.pf[1,1,]<-rnorm(N)
x.pf[2,1,]<-rnorm(N)
w.pf[1,1,]<-1/N
w.pf[2,1,]<-1/N
  
#proportional weight
w.tilde<-matrix(rep(NA,d*N),nrow=d)
  
  for(t in 2:(tau+1)){
    x.pf[1,t,]<-rnorm(N,x.pf[1,t-1,]) #project particles forward by matrix multiplication
    x.pf[2,t,]<-rnorm(N,x.pf[2,t-1,])
    #m<-mean(x.pf[y,])
    
    #w.tilde is col of particles at time step t
    w.tilde[1,]<-w.pf[1,t-1,]*dnorm(y[1,t-1],x.pf[1,t,]) #y is an index of x
    w.tilde[2,]<-w.pf[2,t-1,]*dnorm(y[2,t-1],x.pf[2,t,])
    
    #normalise weights
    w.pf[1,t,]<-w.tilde[1,]/sum(w.tilde[1,])
    w.pf[2,t,]<-w.tilde[2,]/sum(w.tilde[2,])
    
    #compute variance of importance weights
    w.var[1,t-1]<-var(w.pf[1,t,])
    w.var[2,t-1]<-var(w.pf[2,t,])
    
    N.eff[1,t-1]<-1/ sum(w.pf[1,t,]^2)
    N.eff[2,t-1]<-1/ sum(w.pf[2,t,]^2)
    
    if(sum(N.eff[,t-1]<N.thr)==0){
      s1<-sample(1:N,size=N,replace=TRUE,prob=w.pf[1,t,])
      s2<-sample(1:N,size=N,replace=TRUE,prob=w.pf[2,t,])
      #x.pf<-x.pf[,s]
      x.pf[1,t,]<-x.pf[1,t,s1] #resample paths use t if not resampling paths
      x.pf[2,t,]<-x.pf[2,t,s2]
      w.pf[1,t,]<-1/N
      w.pf[2,t,]<-1/N
    
    }
  
  
  
}
  




