#--COMPARE RESAMPLING SCHEMES--
#plot effective sample size
# additionally, counts number of times resample happens

#load library
source("filter_mod.R")
library("ggplot2")

#set model variables
x0<-0
tau<-100
sigma<-1
sigma.meas<-3
N<-500

#generate states and measurements
set.seed(123)
x<-process.pop(tau=tau,x0=x0,sigma=sigma)
y<-meas.pop(x=x,sigma.meas=sigma.meas)

#run filters
obj.mult<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="multinomial",N.thr.per=0.5)
obj.strat<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="stratified",N.thr.per=0.5)
obj.sys<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="systematic",N.thr.per=0.5)
obj.res<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="residual",N.thr.per=0.5) 


#----------------------------------------Effective Sample Size--------------------------------------------

neff.mult<-obj.mult$N.eff.out
neff.strat<-obj.strat$N.eff.out
neff.sys<-obj.sys$N.eff.out
neff.res<-obj.res$N.eff.out

plot(neff.mult,type="l",main="Effective Sample Size")
lines(neff.strat,col="red")
lines(neff.sys,col="blue")
lines(neff.res,col="green")
legend(75, 100, legend=c("Multinomial", "Stratified","Systematic","Residual"),
       col=c( "black","red","blue","green"),lty=c(1,1,1,1), cex=0.8)

#----------------------------------------ACCURACY & DEGNERACY--------------------------------------------


#unique particles(after resampling)
unique.mult<-rep(NA,tau)
unique.strat<-rep(NA,tau)
unique.sys<-rep(NA,tau)
unique.res<-rep(NA,tau)
for(i in 1:tau){
  unique.mult[i]<-length(unique(obj.mult$x.pf.out[i,]))
  unique.strat[i]<-length(unique(obj.strat$x.pf.out[i,]))
  unique.sys[i]<-length(unique(obj.sys$x.pf.out[i,]))
  unique.res[i]<-length(unique(obj.res$x.pf.out[i,]))
}
plot(unique.mult,type="l")
lines(unique.strat,col="red",type="l")
lines(unique.sys,col="blue")
lines(unique.res,col="green")

#run filters(resampling 0.5%)
obj.mult<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="multinomial",N.thr.per=0.5)
obj.strat<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="stratified",N.thr.per=0.5)
obj.sys<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="systematic",N.thr.per=0.5)
obj.res<-particle.filter.pop(N,x,y,x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="residual",N.thr.per=0.5) 

#assign objects
mult.res<-count.resamples(obj.mult$res.point.out)
strat.res<-count.resamples(obj.strat$res.point.out)
sys.res<-count.resamples(obj.sys$res.point.out)
res.res<-count.resamples(obj.res$res.point.out)

#no of resamples 
mult.res$resamples.out
strat.res$resamples.out
sys.res$resamples.out
res.res$resamples.out

#\\not much difference

#plot histogram comparing time step between resamples
name<-c(rep("multinomial",length(mult.res$tt.out)),rep("stratified",length(strat.res$tt.out))
        ,rep("systematic",length(sys.res$tt.out)),rep("residual",length(res.res$tt.out)))
dat<-c(mult.res$tt.out,strat.res$tt.out,sys.res$tt.out,res.res$tt.out)
df<-data.frame(name,dat)

main.plot<-ggplot(df, aes(dat)) + geom_bar(aes(fill = name), position = "dodge")
main.plot+ylab("Frequency") + xlab("Time Steps Between Resamples")+scale_fill_manual(values=c("black", "green", "red","blue"))+labs(fill = "Scheme")









