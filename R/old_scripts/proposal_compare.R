#--COMPARISON OF PROPOSALS--
# random walk example 
# all methods done with 50% resampling
# for each case change N.thr.per to 0 or 1

#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-100
sigma<-2
sigma.meas<-1
N<-500 
M<-100

#generate states and measurements
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)

#generate simulated data y
y.mat<-matrix(NA,nrow=tau,ncol=M)
for(k in 1:M){
  y.mat[,k]<-rand.y.1D(x,sigma.meas=sigma.meas)
}


#-------------------------------------PRIOR -----------------------------------

#set up variables
N.eff.prior<-matrix(NA,tau,M)
MSE.prior<-matrix(NA,tau,M)
means.prior<-matrix(NA,tau,M)

#run pf with prior proposals
for(k in 1:M){
  obj.prior<-particle.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)
  MSE.prior[,k]<-obj.prior$MSE.k.out
  N.eff.prior[,k]<-obj.prior$N.eff.out
  means.prior[,k]<-obj.prior$m.out
}

#extract info
MSE.prior.time<-rowMeans(MSE.prior)
AMSE.prior<-mean(MSE.prior)
var.prior<-var(means.prior[40,])
N.eff.prior<-rowMeans(N.eff.prior)


#-------------------------------------OPTIMUM -----------------------------------

#set up variables
N.eff.optimum<-matrix(NA,tau,M)
MSE.optimum<-matrix(NA,tau,M)
means.optimum<-matrix(NA,tau,M)

# run pf with optimum proposal 
for(k in 1:M){
  obj.optimum<-optimum.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)
  MSE.optimum[,k]<-obj.optimum$MSE.k.out
  N.eff.optimum[,k]<-obj.optimum$N.eff.out
  means.optimum[,k]<-obj.optimum$m.out
}

#extract info
MSE.optimum.time<-rowMeans(MSE.optimum)
AMSE.optimum<-mean(MSE.optimum)
var.optimum<-var(means.optimum[40,])
N.eff.optimum<-rowMeans(N.eff.optimum)


#-------------------------------------FIXED-----------------------------------


#set up variables
N.eff.fixed<-matrix(NA,tau,M)
MSE.fixed<-matrix(NA,tau,M)
means.fixed<-matrix(NA,tau,M)

# run pf  with fixed proposal
for(k in 1:M){
  obj.fixed<-fixed.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)
  MSE.fixed[,k]<-obj.fixed$MSE.k.out
  N.eff.fixed[,k]<-obj.fixed$N.eff.out
  means.fixed[,k]<-obj.fixed$m.out
}

#extract info
MSE.fixed.time<-rowMeans(MSE.fixed)
AMSE.fixed<-mean(MSE.fixed)
var.fixed<-var(means.fixed[40,])
N.eff.fixed<-rowMeans(N.eff.fixed)


#-------------------------------------LIKELIHOOD-----------------------------------

#set up variables
N.eff.like<-matrix(NA,tau,M)
MSE.like<-matrix(NA,tau,M)
means.like<-matrix(NA,tau,M)

#run pf with likelihood proposal
for(k in 1:M){
  obj.like<-likelihood.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)
  MSE.like[,k]<-obj.like$MSE.k.out
  N.eff.like[,k]<-obj.like$N.eff.out
  means.like[,k]<-obj.like$m.out
}
MSE.like.time<-rowMeans(MSE.like)
AMSE.like<-mean(MSE.like)
var.like<-var(means.like[40,])
N.eff.like<-rowMeans(N.eff.like)

#-------------------------------------RESULTS -----------------------------------

#-------- N.thr.per= 0 --------
#Scroll down for plotting setup N.thr.per=1

#AMSE results
AMSE.prior 
AMSE.optimum 
AMSE.fixed
AMSE.like 

#plot means at t=40
plot(x,main="Proposal Comparison",xlab="Time,t")
lines(means.prior[,40])
lines(means.like[,40],col="orange")
lines(means.optimum[,40],col="blue")
lines(means.fixed[,40],col="red")
legend(20, 20, legend=c("PF Prior", "Optimum","Fixed","Likelihood"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)

#plot MSE
plot(MSE.prior.time,type="l",main="MSE",ylab="MSE",xlab="Time,t")
lines(MSE.optimum.time,col="blue")
lines(MSE.like.time,col="orange")
lines(MSE.fixed.time,col="red")
legend(60, 60, legend=c("PF Prior", "Optimum","Fixed","Likelihood"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)


#plot comparison effective  sample size
plot(N.eff.prior,type="l",main="Effective Sample Size",xlab="Time,t",ylab="Effective Sample Size",ylim=c(0,500))
lines(N.eff.fixed,col="red")
lines(N.eff.like,col="orange")
lines(N.eff.optimum,col="blue")
legend(70, 400, legend=c("PF Prior", "Optimum","Fixed","Likelihood"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)


#-------- N.thr.per= 1 --------


#AMSE results
AMSE.prior #0.8131669
AMSE.optimum # 0.8099655
AMSE.fixed #10.69816
AMSE.like #0.8186595

#plot means at t=40
plot(x,main="Proposal Comparison",xlab="Time,t")
lines(means.prior[,40])
lines(means.like[,40],col="orange")
lines(means.optimum[,40],col="blue")
lines(means.fixed[,40],col="red")
legend(20, 20, legend=c("PF Prior", "Optimum","Fixed","Likelihood"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)

#plot MSE
plot(MSE.fixed.time,col="red",type="l",main="MSE",ylab="MSE",xlab="Time,t",ylim=c(0,40))
lines(MSE.optimum.time,col="blue")
lines(MSE.like.time,col="orange")
lines(MSE.prior.time)
legend(60, 30, legend=c("PF Prior", "Optimum","Fixed","Likelihood"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)


#plot comparison effective  sample size
plot(N.eff.prior,type="l",main="Effective Sample Size",xlab="Time,t",ylab="Effective Sample Size",ylim=c(0,500))
lines(N.eff.fixed,col="red")
lines(N.eff.like,col="orange")
lines(N.eff.optimum,col="blue")
legend(70, 100, legend=c("PF Prior", "Optimum","Fixed","Likelihood"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)


