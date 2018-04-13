#--COMPUTE VARIANCE FOR DIFF RESAMPLING THRESHOLDS--
# more resampling, should increase variance
# two ways: compute var of M simulations of mean estimate or compute central limit variance
# non-convincing results
# can obtain results by comparing resampling vs no resampling 

#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-100
sigma<-3
sigma.meas<-1
N<-500

#generate states and measurements
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.meas=sigma.meas)

#set up variables 
M<-10
mean_vec<-rep(NA,M)
N.thr.vec<-seq(0.1,1,0.01)

#calculate variance of mean estimate at t=40 for diff resampling times
mat<-matrix(NA,nrow=M,ncol=length(N.thr.vec))
for(k in 1:length(N.thr.vec)){
  for(j in 1:M){
    obj.prior<-particle.filter(N=N,x=x,y=y,x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=N.thr.vec[k])
    mean_vec[j]<-obj.prior$m.out[40]
  }
  mat[,k]<-mean_vec
}

#calculate central limit variance at t=40 at diff resampling times
var_vec<-rep(NA,length(N.thr.vec))
for(j in 1:length(N.thr.vec)){
  obj.prior<-particle.filter(N=N,x=x,y=y,x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=N.thr.vec[k])
  var_vec[j]<-obj.prior$var.s.out[40]
}

#plot both variance
plot(apply(mat,2,var),type="l",main="Monte Carlo Variance",xlab="Resampling Threshold",ylab="Variance of mean")
lines(var_vec/N,col="red",lty=c(1,1), cex=0.8)
legend(60, 0.020, legend=c("Simulated variance","Central Limit variance"),
       col=c("black", "red"),lty=c(1,1), cex=0.8)




