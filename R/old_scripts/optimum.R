#load library
source("resample.R")
source("random-walk.R")
source("filter.R")

#set model variables
x0<-0
sigma<-1
sigma.meas<-1
tau<-100

#generate states and measurements
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.y=sigma.meas)

#run particle filter
obj.mult<-optimum.filter(N=500,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="multinomial",N.thr.per=0)

#plot
m<-obj.mult$m.out
x.pf<-obj.mult$x.pf.out
plot(x)
points(y,col="red")
lines(m)
lines(x.pf[,1])
