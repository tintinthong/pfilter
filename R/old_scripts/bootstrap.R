#--BOOTSTRAP WORKING FILE--
# bootstrap for random walk +credible intervals

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
                          resample.type="standard",N.thr.per=1)

#extract info
m<-obj.randw$m.out
conf<-obj.randw$conf.out


#plot
plot(x,main=paste("Random Walk 1D (",N, "particles) " ),xlab="Time,t",ylim=c(-5,10))
points(y,col="red")
lines(m,col="blue")
lines(conf[,1],lty=3)
lines(conf[,2],lty=3)
legend(10, 10, legend=c("States, x", "Measurements, y","Mean, m"),
       col=c("black", "red","blue"), pch=c(1,1,NA),lty=c(NA,NA,1), cex=0.8)











