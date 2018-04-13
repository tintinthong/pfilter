
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
#x<-process.pop(tau=tau,x0=x0,sigma=sigma)
#y<-meas.pop(x,sigma.meas=sigma.meas)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.meas=sigma.meas)

#run particle filter
#obj.randw<-particle.filter.pop.path(N=N,x=x,y=y,x0=5,sigma=sigma,sigma.meas=sigma.meas,
#                           resample.type="standard",N.thr.per=0.5)
obj.randw<-particle.filter.path(N=N,x=x,y=y,x0=5,sigma=sigma,sigma.meas=sigma.meas,
                                    resample.type="standard",N.thr.per=1)
#obj.randw<-particle.filter(N=N,x=x,y=y,x0=5,sigma=sigma,sigma.meas=sigma.meas,
#                                    resample.type="standard",N.thr.per=0)


#extract info
m<-obj.randw$m.out
x.pf<-obj.randw$x.pf.out
conf<-obj.randw$conf.out
N.eff<-obj.randw$N.eff.out
res.point<-obj.randw$res.point.out

plot(x)#,ylim=c(-10,15) 
points(y,col="red")

for(i in 1:N){
  points(x.pf[,i])
}

un<-rep(NA,tau)
for(t in 1:tau){
  un[t]<-length(unique(x.pf[t,]))
}
plot(un,type="b")

hist(x.pf[1,])
hist(x.pf[25,])
hist(x.pf[50,])
hist(x.pf[tau,])




plot(x)
lines(m)
plot(N.eff,type="b")




#------------------------------------------------------------------------
#main plot
plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
points(y,col="red")
legend(20, 10, legend=c("States, x", "Measurements, y"),
       col=c("black", "red"), pch=c(1,1), cex=0.8)

#plot mean over states
plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
lines(m,col="blue")
for(i in 1:N){
  lines(x.pf[,i])
}

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
  y<-rand.y.1D(x,sigma.meas=sigma.meas)
  obj.randw<-particle.filter(N=N,x=x,y=y,x0=-5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
  MSE[,k]<-obj.randw$MSE.k.out
}
RMSE<-sqrt(colMeans(MSE))
ARMSE<-mean(RMSE)
ARMSE
MSE<-rowMeans(MSE)


#plot MSE

plot(MSE)