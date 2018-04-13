#--Compare performance for different N particles--
#\\ bootstrap for random walk
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

#compute ARMSE for each N
N.list<-seq(10,500,by=10)
ARMSE.list<-rep(NA,length(N.list))
M<-10
MSE<-matrix(NA,tau,M)
for(i in 1:length(N.list)){
  Ni<-N.list[i]
  for(k in 1:M){
    obj.randw<-particle.filter(N=Ni,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
    MSE[,k]<-obj.randw$MSE.k.out
  }
  ARMSE<-mean(sqrt(colSums(MSE)/tau))
  ARMSE.list[i]<-ARMSE
}
plot(N.list,ARMSE.list,main=paste("Average RMSE(", M ,"MC runs)"),xlab="N",ylab="ARMSE")


#plot Effective Sample Size (average percentage effective sample size drops)
#\\varince of weights drop
N.list<-seq(10,500,by=10)
AN.eff.list<-rep(NA,length(N.list))
N.eff<-matrix(NA,tau,M)
M<-10
for(i in 1:length(N.list)){
  Ni<-N.list[i]
  for(k in 1:M){
    obj.randw<-particle.filter.pop(N=Ni,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
    N.eff[,k]<-obj.randw$N.eff.out
  }
  AN.eff<-colSums(N.eff/Ni)/tau
  AN.eff.list[i]<-mean(AN.eff)
}
plot(N.list,AN.eff.list,main=paste("% Effective Sample Size(", M ,"MC runs)"),xlab="N",ylab="Average %N.eff")


