#--SMOOTHING DEGENERACY PLOT--
# as t increases, paths coalesce
# seeds must be set in order

#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-200
sigma<-1
sigma.meas<-1
N<-500

#generate states and simulated data y
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.meas=sigma.meas)

#smoothed estimate
obj<-particle.filter.path(N=N,x=x[1:100],y=y[1:100],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=1)
obj$x.pf.out

#plot ten
plot(x,xlab="Time,t",ylab="State,x",main=paste("Bootstrap-Random Walk(",N,"particles)"))
for(i in 1:N){
  lines(obj$x.pf.out[,i],col="blue")
}
legend(150, 10, legend=c("States, x", "Mean Estimate"),
       col=c("black", "blue"),pch=c(1,NA),lty=c(NA,1), cex=0.8)


#generate states and simulated data y
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.meas=sigma.meas)

#smoothed estimate
obj2<-particle.filter.path(N=N,x=x[1:200],y=y[1:200],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=1)
obj2$x.pf.out

#plot ten
plot(x,xlab="Time,t",ylab="State,x",main=paste("Bootstrap-Random Walk(",N,"particles)"))
for(i in 1:N){
  lines(obj2$x.pf.out[,i],col="blue")
}
legend(150, 10, legend=c("States, x", "Mean Estimate"),
       col=c("black", "blue"),pch=c(1,NA),lty=c(NA,1), cex=0.8)












