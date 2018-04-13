#--PLOT MSE OF FILTERED VS SMOOTH--
# M monte carlo runs
# Test 2 cases sigma.meas=3 and sigma.meas=1

#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-100
sigma<-3
sigma.meas<-3 #sigma.meas=1
N<-500
M<-100

#generate states and simulated data y
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y.mat<-matrix(NA,nrow=tau,ncol=M)
for(k in 1:M){
  y.mat[,k]<-rand.y.1D(x,sigma.meas=sigma.meas)
}

#compute MSE for smoothing and filtered
MSE<-matrix(NA,tau,M)
MSE.filter<-matrix(NA,tau,M)
means<-matrix(NA,tau,M)
means.filter<-matrix(NA,tau,M)
neff<-matrix(NA,tau,M)
for(k in 1:M){
  obj<-particle.filter.path(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=1)
  x.pf<-obj$x.pf.out
  means.filter[,k]<-obj$m.out
  means[,k]<-rowMeans(x.pf)
  MSE[,k]<-(rowMeans(x.pf)-x)^2
  MSE.filter[,k]<-(obj$m.out-x)^2
  neff[,k]<-obj$N.eff.out
}

#plot MSE for smoothing and filtered
plot(rowMeans(MSE),main="MSE",ylab="MSE",xlab="Time,t",type="l")
lines(rowMeans(MSE.filter),col="red",type="l")
legend(10, 14, legend=c("Filtered", "Joint Smoothing"),
       col=c( "red","black"),pch=c(1,1), cex=0.8)
#legend(75, 1.8, legend=c("Filtered", "Joint Smoothing"),
#       col=c( "red","black"),pch=c(1,1), cex=0.8)



