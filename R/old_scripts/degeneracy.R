#bootstrap 
#- always set working directory to file location
#- modules meant to be sourced ie functions are dependent on source("module.R"). Do not attempt to use just one function

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
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.y=sigma.meas)


#filter should contain x0 and x for all t
#filter should not contain y for all t but not for x0 

#set up filter
x0.pf<-rnorm(N,mean=x0,sd=sigma) #this is not correct sigma
x.pf <- matrix(NA,ncol=N,nrow=tau)
w.pf<-matrix(NA,ncol=N,nrow=tau) #primary weights
w0.pf<-rep(1/N,N)

#set up diagnostics
diagnostic.init(tau)

# 1. Initialize
x.pf[1, ] <- x0.pf+rnorm(N,sd=sigma) #project particles from x0 to the first time step
w.tilde<- dnorm(y[1],mean=x.pf[1,],sd=sigma.meas)*w0.pf #the previous weight would be equal weights 1/N depending on how x0 was chosen
w.pf[1,]<-w.tilde/sum(w.tilde)
#x.pf[1,]<-mult.resample(N=N,x=x.pf[1,],w=w) 
#w<-rep(1/N,N)
m[1]<-sum(w.pf[1,]*x.pf[1,])
MSE[1]<-mean((x[1]-x.pf[1,])^2)
N.eff[1]<-1/sum(w^2)
var.w[1]<-var(w)
conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))

for (t in 2:(tau)) {
  
  #projec  state forward
  x.pf[t, ] <-x.pf[t-1,]+rnorm(N,sd=sigma)
  
  #compute weights
  w.tilde <- dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas)*w.pf[t-1,] #when u add w resampling does not help
  w.pf[t,] <- w.tilde/sum(w.tilde) #we run out of numbers to represent
  
  #effective sample size
  N.eff[t]<-1/sum(w^2)
  
  
  #print(N.eff[t])
  #if(N.eff[t]<N.thr){
  #  x.pf[t,]<-mult.resample(N,x=x.pf[t,],w=w)
    
    #x.pf[t,]<-strat.resample(x=x.pf[t,],w=w)
    #x.pf[t,]<-sys.resample(x=x.pf[t,],w=w)
    #x.pf[t,]<-res.resample(x=x.pf[t,],w=w) #there is a bug here
  #  w<-rep(1/N,N)
  #}
  
  #calculate diagnostics
  m[t]<-sum(w.pf[t,]*x.pf[t,])
  MSE[t]<-mean((x[t]-x.pf[t,])^2)
  var.w[t]<-var(w) #this does not mean anything
  conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
  
}

#diagnostic plots
plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
points(y,col="red")
legend(20, 10, legend=c("States, x", "Measurements, y"),
       col=c("black", "red"),pch=c(1,1), cex=0.8)

plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
lines(m,col="blue")

hist(w.pf[1,],ylim=c(0,360))
hist(w.pf[2,],ylim=c(0,360))
hist(w.pf[3,],ylim=c(0,360))
hist(w.pf[4,],ylim=c(0,360))

plot(MSE)





#for t=28, all elements of w become very very small, such that there is only 1 element with non-zero value
# this is even with re sampling or no resampling
# we end up dividing a very small number, by another very small number to get normaility (this is not a big problem)
# but this turns out bad, 
# after we project forward, the likelihood of the non-zero weighted particle possibly becomes close to 0 
#this happens when we used peaked likelihood, and multiply the weight by previous weight value

#when likelihood is peaked then variance of weights get really high really quickly
#for better estimation, for marginal with multi-modularity, we need likelihood to be a little broader

#if(t==28){
#  print(dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas)*w)
#  print(paste("hi",sum(w^2)))
#  #divide a very small number, with a very small number
#}
