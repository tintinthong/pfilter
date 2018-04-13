#--SIS PLOT--
#1) Random Walk mu=1 sigma=1

#load library
source("filter_mod.R")

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
obj.randw<-particle.filter(N=N,x=x,y=y,x0=5,sigma=sigma,sigma.meas=sigma.meas,
                           resample.type="standard",N.thr.per=0)

#extract info
m<-obj.randw$m.out

#plot states, measurements  and mean estimates
plot(x,main=paste("Random Walk 1D (",N, "particles) " ),xlab="Time,t")
points(y,col="red")
lines(m,col="blue")
legend(20, 10, legend=c("States, x", "Measurements,y", "Mean Estimate,m"),
       col=c("black", "red","blue"), pch=c(1,NA,NA),lty=c(NA,1,1), cex=0.8)
