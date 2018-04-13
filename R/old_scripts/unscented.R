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

#generate true state and data
x0<-0
x<-process.1D(tau=tau,x0=x0,sigma=sigma)
y<-extended.meas.1D(x0=x0,x=x,sigma.y=sigma.meas)
#filter should contain x0 and x for all t
#filter should not contain y for all t but not for x0 

#set up filter
x0.pf<-rnorm(N,mean=x0,sd=sigma)
x.pf <- matrix(NA,ncol=N,nrow=tau)

#set up diagnostics
diagnostic.init(tau)

# 1. Initialize
x.pf[1, ] <- x0.pf+rnorm(N,sd=sigma) #project particles from x0 to the first time step
w.tilde<- dnorm(y[1],mean=x.pf[1,],sd=sigma.meas) #the previous weight would be equal weights 1/N depending on how x0 was chosen
w<-w.tilde/sum(w.tilde)
x.pf[1,]<-mult.resample(N=N,x=x.pf[1,],w=w) 
m[1]<-sum(w*x.pf[1,])
MSE[1]<-mean((x[1]-x.pf[1,])^2)
N.eff[1]<-1/sum(w^2)
var.w[1]<-var(w)
conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))

for (t in 2:(tau)) {
  
  #projec  state forward
  x.pf[t, ] <-x.pf[t-1,]+rnorm(N,sd=sigma)
  
  #compute weights
  w.tilde <- dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas) #when u add w resampling does not help
  w <- w.tilde/sum(w.tilde) #we run out of numbers to represent
  
  #effective sample size
  N.eff[t]<-1/sum(w^2)
  
  
  #print(N.eff[t])
  if(N.eff[t]<N.thr){
    x.pf[t,]<-mult.resample(N=N,x=x.pf[t,],w=w) 
    #x.pf[t,]<-strat.resample(x=x.pf[t,],w=w)
    #x.pf[t,]<-sys.resample(x=x.pf[t,],w=w)
    #x.pf[t,]<-res.resample(x=x.pf[t,],w=w) #there is a bug here
    w<-rep(1/N,N)
  }
  
  #calculate diagnostics
  m[t]<-sum(w*x.pf[t,])
  MSE[t]<-mean((x[t]-x.pf[t,])^2)
  var.w[t]<-var(w) #this does not mean anything
  conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
  
}

#diagnostic plots
plot(x)
points(y,col="red")

plot(x)
lines(m,col="red")

plot(x)
points(x.pf[,1],col="red")
