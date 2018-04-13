#--PLOT VARIANCE OF WEIGHTS THROUGH TIME & PLOT WEIGHTS THROUGH TIME--
# Using SIS(no resampling)

#load library
source("filter_mod.R")

#set model variables
N<-30
x0<-0
tau<-100
sigma<-1
sigma.meas<-3

#generate states and simulated data y
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.meas=sigma.meas)

#run filter
obj<-particle.filter(N=N,x=x,y=y,x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=0)

#extract info
x.pf<-obj$x.pf.out
w.mat<-obj$w.mat.out

#plot variance of weights
plot(apply(w.mat,1,var),type="l",ylab="Variance",xlab="Time,t",main="Variance of weights")

#plot weights
plot(w.mat[,1],type="l",ylim=c(0,1),ylab="Weight Value",xlab="Time,t",main="Weights")
for(i in 2:N){
  lines(w.mat[,i])
}



