#Rejection Filter
#//can use rejection as oppose to importance sampling 
#//the ideas are similar because weights are a way to compare posterior from the proposal 

#load library
source("filter_mod.R")


#---------------------------------------RANDOM WALK EXAMPLE------------------------------------------

#set model variables
x0<-0
tau<-100
sigma<-1
sigma.meas<-1
N<-100
Mnorm<-0.4 #about the max of a normal distribution
#\\weird:don't know why Mnorm close to 0.4 gives me an instance a being very small
#\\it is just random 

#generate states and measurements
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.meas=sigma.meas)

#run particle filter
obj.randw<-rejection.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,Mnorm)

#extract info
m<-obj.randw$m.out
x.pf<-obj.randw$x.pf.out
conf<-obj.randw$conf.out
N.eff<-obj.randw$N.eff.out
res.point<-obj.randw$res.point.out

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


#compute RMSE
M<-50
MSE<-matrix(NA,tau,M)
for(k in 1:M){
  print(k)
  obj.randw<-rejection.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,Mnorm)
  MSE[,k]<-obj.randw$MSE.k.out
}
RMSE<-sqrt(colSums(MSE)/tau)
ARMSE<-mean(RMSE)
ARMSE

#plot RMSE
plot(RMSE)



#f<-function(zp){
#  dnorm(y[1],x.old[i],sd=sigma.meas)*dnorm(zp,x.old[i], sd=sigma)/dnorm(zp,x.old[i],sd=sigma) 
#}#unsure whether denominator is x.old or x0.pf[]
#M<-0.4
#optim(x.old[i],f,method="Brent",lower=max(x.old)-3*sigma,upper=max(x.old)+3*sigma)$par



#---------------------------------------POPULATION EXAMPLE------------------------------------------

#set model variables
x0<-0 
tau<-100
sigma<-1
sigma.meas<-3
N<-100
Mnorm<-0.4

#generate states and measurements
set.seed(123)
x<-process.pop(tau=tau,x0=x0,sigma=sigma)
y<-meas.pop(x=x,sigma.meas=sigma.meas)


#run particle filter
obj.pop<-rejection.filter.pop(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,Mnorm)

#extract info
m<-obj.pop$m.out
x.pf<-obj.pop$x.pf.out
conf<-obj.pop$conf.out
N.eff<-obj.pop$N.eff.out
res.point<-obj.pop$res.point.out

#diagnostic plots

#main plot
plot(x,main=paste("Population (",N, "particles)"),xlab="t")
points(y,col="red")
legend(40, -10, legend=c("States, x", "Measurements, y"),
       col=c("black", "red"),pch=c(1,1), cex=0.8)

#plot mean over states
plot(x,main=paste("Population (",N, "particles)"),xlab="t")
lines(m,col="blue")

#plot confidence intervals

#ggplot
#library(ggplot2)
#df <- data.frame(t=1:tau,mu=x,L=conf[,1],U=conf[,2])
#ggplot(df, aes(x = t, y = mu)) +geom_point(size = 1) +geom_errorbar(aes(ymax = U, ymin = L))

#base package
plot(x,main=paste("Population (","N", "particles)"),xlab="t")
for(i in 1:tau){
  segments(i,conf[i,1],i,conf[i,2])
}
plot(density(x.pf[4,]))
hdi(x.pf[4,])



#compute RMSE
M<-50
MSE<-matrix(NA,tau,M)
for(k in 1:M){
  print(k)
  obj.pop<-rejection.filter.pop(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,Mnorm)
  MSE[,k]<-obj.pop$MSE.k.out
}
RMSE<-sqrt(colSums(MSE)/tau)
ARMSE<-mean(RMSE)




