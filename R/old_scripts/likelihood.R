#--LIKELIHOOD WORKING FILE--
#- always set working directory to file location
#- modules meant to be sourced ie functions are dependent on source("module.R"). Do not attempt to use just one function
#\\1) bootstrap for random walk
#\\2)bootstrap for population model

#\\always set working directory to file location
#\\ modules meant to be sourced ie functions are dependent on source("module.R")- Do not attempt to use just one function 


#load library
source("filter_mod.R")

#---------------------------------------RANDOM WALK EXAMPLE------------------------------------------

#set model variables
x0<-0
tau<-100
sigma<-1
sigma.meas<-1
N<-500

#generate states and measurements
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.meas=sigma.meas)


#run particle filter
obj.like<-auxiliary.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,
                                resample.type="standard",N.thr.per=0.5)


#extract info
m<-obj.like$m.out
x.pf<-obj.like$x.pf.out
conf<-obj.like$conf.out
N.eff<-obj.like$N.eff.out
res.point<-obj.like$res.point.out

#main plot
plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
points(y,col="red")
legend(20, 10, legend=c("States, x", "Measurements, y"),
       col=c("black", "red"), pch=c(1,1), cex=0.8)


#plot mean over states
plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
lines(m,col="blue")


#plot credible intervals

#ggplot
#library(ggplot2)
#df <- data.frame(t=1:tau,mu=x,L=conf[,1],U=conf[,2])
#ggplot(df, aes(x = t, y = mu)) +geom_point(size = 1) +geom_errorbar(aes(ymax = U, ymin = L))

#base package
plot(x,main=paste("Random Walk 1D (","N", "particles)"),xlab="t")
for(i in 1:tau){
  segments(i,conf[i,1],i,conf[i,2])
}

#plot Effective
plot(N.eff,type="b")

#plot resampled points
plot(res.point,type="b")

#compute RMSE
M<-50
MSE<-matrix(NA,tau,M)
for(k in 1:M){
  obj.randw<-likelihood.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)
  MSE[,k]<-obj.randw$MSE.k.out
}
RMSE<-sqrt(colSums(MSE)/tau)
ARMSE<-mean(RMSE)


