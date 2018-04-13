#--COMPARISON OF DIFFERENT FILTERS(RANDOM WALK)--
#two cases: (sigma, sigma.meas)=(3,1) or (1,3)
# all methods resample every time


#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-100
sigma<-1 #1
sigma.meas<-3 #3
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

#------------------------------------BOOTSTRAP---------------------------------------
# resamples at every step

# run filter M times
MSE.boot.mat<-matrix(NA,tau,M)
N.eff.boot<-matrix(NA,tau,M)
for(k in 1:M){
  obj.boot<-particle.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
  MSE.boot.mat[,k]<-obj.boot$MSE.k.out
  N.eff.boot[,k]<-obj.boot$N.eff.out
}

#compute MSE  and AMSE 
MSE.boot<-colMeans(MSE.boot.mat)
AMSE.boot<-mean(MSE.boot)
MSE.boot<-rowMeans(MSE.boot.mat)
N.eff.boot<-rowMeans(N.eff.boot)

#assign mean
m.boot<-obj.boot$m.out

#------------------------------------Exact---------------------------------------
# no resampling for this filter

# run filter M times
MSE.exact.mat<-matrix(NA,tau,M)
N.eff.exact<-matrix(NA,tau,M)
for(k in 1:M){
  obj.exact<-exact.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard")
  MSE.exact.mat[,k]<-obj.exact$MSE.k.out
  N.eff.exact[,k]<-obj.exact$N.eff.out
}

#compute MSE  and AMSE 
MSE.exact<-colMeans(MSE.exact.mat)
AMSE.exact<-mean(MSE.exact)
MSE.exact<-rowMeans(MSE.exact.mat)
N.eff.exact<-rowMeans(N.eff.exact)


#assign mean
m.exact<-obj.exact$m.out


#------------------------------------Optimum---------------------------------------
# resample 50% of the time 


#compute MSE
MSE.opt.mat<-matrix(NA,tau,M)
N.eff.opt<-matrix(NA,tau,M)
for(k in 1:M){
  obj.opt<-optimum.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=1)
  MSE.opt.mat[,k]<-obj.opt$MSE.k.out
  N.eff.opt[,k]<-obj.opt$N.eff.out
}


#compute MSE  and AMSE
MSE.opt<-colMeans(MSE.opt.mat)
AMSE.opt<-mean(MSE.opt)
MSE.opt<-rowMeans(MSE.opt.mat)
N.eff.opt<-rowMeans(N.eff.opt)

#assign mean
m.opt<-obj.opt$m.out


#------------------------------------Likelihood---------------------------------------
# resample 50% of the time


#compute MSE
MSE.like.mat<-matrix(NA,tau,M)
N.eff.like<-matrix(NA,tau,M)
for(k in 1:M){
  obj.like<-likelihood.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
  MSE.like.mat[,k]<-obj.like$MSE.k.out
  N.eff.like[,k]<-obj.like$N.eff.out
}


#compute MSE  and AMSE
MSE.like<-colMeans(MSE.like.mat)
AMSE.like<-mean(MSE.like)
AMSE.like 
MSE.like<-rowMeans(MSE.like.mat)
N.eff.like<-rowMeans(N.eff.like)

#assign mean
m.like<-obj.like$m.out


#------------------------------------RESULTS---------------------------------------

#------(sigma, sigma.meas)=(3,1)-------

#display AMSE
AMSE.like #2.79889(1,3) 0.8977392(3,1) 
AMSE.opt #2.611359 (1,3) 0.8954476 (3,1)  
AMSE.boot #2.625013(1,3) 0.8964912 (3,1) 
AMSE.exact #2.626761 (1,3) 0.894724 (3,1) 


#plot means(There is hardly any difference)
plot(x,main=paste("Random Walk (",N, "particles)"),xlab="t")
lines(m.boot)
lines(m.opt,col="blue")
lines(m.like, col="red")
lines(m.exact, col="orange")
legend(20, 25, legend=c("Bootstrap", "PF Optimum","Likelihood","Exact"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)


#plot effective sample size
#obj.boot$N.eff.out
plot(N.eff.boot,ylab="N.eff",xlab="Time,t", main=" Effective Sample Size",type="l",ylim=c(1,500))
lines(N.eff.like,col="red")
lines(N.eff.opt,col="blue")
legend(70, 100, legend=c("Bootstrap", "PF Optimum","Likelihood"),
       col=c("black", "blue","red"),lty=c(1,1,1,1), cex=0.8)

#plot MSE(There is hardly any difference)
plot(MSE.boot,main="MSE",xlab="Time,t",type="l",ylim=c(0,5))
lines(MSE.like,col="red")
lines(MSE.opt,col="blue")
lines(MSE.exact,col="orange")

#------(sigma, sigma.meas)=(1,3)-------


#display AMSE
AMSE.like #2.745584 (1,3) 0.8984577(3,1)
AMSE.opt #2.629193 (1,3) 0.8950922(3,1)
AMSE.boot #2.625013(1,3) 0.8964912(3,1)
AMSE.exact #2.626761(1,3) 0.894724(3,1)

#old results(keep just in case)
#AMSE.like #2.79889(1,3) 0.8977392(3,1) 0.6573997(1,1)
#AMSE.opt #2.611359 (1,3) 0.8954476 (3,1)  0.623635(1,1)
#AMSE.boot #2.625013(1,3) 0.8964912 (3,1) 0.6381266(1,1)
#AMSE.exact #2.626761 (1,3) 0.894724 (3,1) 0.6228132(1,1)

#plot means(There is hardly any difference)
plot(x,main=paste("Random Walk (",N, "particles)"),xlab="t")
lines(m.boot)
lines(m.opt,col="blue")
lines(m.like, col="red")
lines(m.exact, col="orange")

#plot effective sample size
plot(N.eff.boot,ylab="N.eff",xlab="Time,t", main=" Effective Sample Size",type="l",ylim=c(1,500))
lines(N.eff.like,col="red")
lines(N.eff.opt,col="blue")
legend(70, 100, legend=c("Bootstrap", "PF Optimum","Likelihood"),
       col=c("black", "blue","red"),lty=c(1,1,1), cex=0.8)

#plot MSE(There is hardly any difference)
plot(MSE.boot,main="MSE",xlab="Time,t",type="l",ylim=c(0,5))
lines(MSE.like,col="red")
lines(MSE.opt,col="blue")
lines(MSE.exact,col="orange")

