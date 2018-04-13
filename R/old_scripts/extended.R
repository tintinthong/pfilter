#load libraries
source("resample.R")
source("random-walk.R")
source("bootstrap_mse.R")

#set model variables
x0<-0
sigma<-1
sigma.meas<-1
tau<-100

#generate states and measurements
set.seed(123)
x<-process.1D(tau=tau,x0=x0,sigma=sigma)
y<-extended.meas.1D(x0=x0,x=x,sigma.y=sigma.meas)

#particle filter
obj.mult<-particle.filter.pop(N=500,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="multinomial",N.thr.per=1)

m<-obj.mult$m.out

#plot state and measurement 
plot(x,main="Plot States and Measurements",xlab="Time,t")
points(y,col="red")
legend(20, -10, legend=c("States, x", "Measurements, y"),
       col=c("black", "red"),pch=c(1,1), cex=0.8)


#plot line of mean estimate 
plot(x,main="Plot States and Estimate",xlab="Time,t")
lines(m,col="blue")
legend(20, -10, legend=c("States, x","Mean Estimate,m"),
       col=c("black","blue"),pch=c(1,NA),lty=c(NA,1), cex=0.8)



