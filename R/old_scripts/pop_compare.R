#--COMPARISON OF DIFFERENT FILTERS(NON LINEAR EXAMPLE)--
#two cases: (sigma, sigma.meas)=(3,1) or (1,3)
# all methods resample when neff drops below 50% of N

#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-100
sigma<-3 #1
sigma.meas<-1 #3
N<-500
M<-100

#generate states and measurements
set.seed(123)
x<-process.pop(tau=tau,x0=x0,sigma=sigma)

#generate simulated data y
y.mat<-matrix(NA,nrow=tau,ncol=M)
for(k in 1:M){
  y.mat[,k]<-meas.pop(x=x,sigma.meas=sigma.meas)
}

#------------------------------------Particle Filter ---------------------------------------
#prior proposal


# run filter M times
MSE.pf.mat<-matrix(NA,tau,M)
N.eff.pf<-matrix(NA,tau,M)
for(k in 1:M){
  obj.pf<-particle.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)
  MSE.pf.mat[,k]<-obj.pf$MSE.k.out
  N.eff.pf[,k]<-obj.pf$N.eff.out
}

#compute MSE and AMSE
MSE.pf<-colMeans(MSE.pf.mat)
AMSE.pf<-mean(MSE.pf)
MSE.pf<-rowMeans(MSE.pf.mat)
N.eff.pf<-rowMeans(N.eff.pf)

#assign mean
m.pf<-obj.pf$m.out




#-----------------------------------Regularized---------------------------------------

#run filter M times
MSE.reg.mat<-matrix(NA,tau,M)
N.eff.reg<-matrix(NA,tau,M)
for(k in 1:M){
  obj.reg<-regularized.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=0.5)
  MSE.reg.mat[,k]<-obj.reg$MSE.k.out
  N.eff.reg[,k]<-obj.reg$N.eff.out
}

#compute MSE and AMSE
MSE.reg<-colMeans(MSE.reg.mat)
AMSE.reg<-mean(MSE.reg)
MSE.reg<-rowMeans(MSE.reg.mat)
N.eff.reg<-rowMeans(N.eff.reg)

#assign mean
m.reg<-obj.reg$m.out



#------------------------------------Rejection---------------------------------------
#extremely slow(need to leave over night if results need to be checked)

#run filter M times
MSE.rej.mat<-matrix(NA,tau,M)
for(k in 1:M){
  obj.rej<-rejection.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,Mnorm=dnorm(0,0,sigma.meas))
  MSE.rej.mat[,k]<-obj.rej$MSE.k.out
  print(paste("k",k)) #print it to see progress
}

#compute MSE and AMSE
MSE.rej<-colMeans(MSE.rej.mat)
AMSE.rej<-mean(MSE.rej)
MSE.rej<-rowMeans(MSE.rej.mat)

#assign mean
m.rej<-obj.rej$m.out



#------------------------------------Auxiliary---------------------------------------


#run filter M times
MSE.aux.mat<-matrix(NA,tau,M)
N.eff.aux<-matrix(NA,tau,M)
for(k in 1:M){
  obj.aux<-auxiliary.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas)
  MSE.aux.mat[,k]<-obj.aux$MSE.k.out
  N.eff.aux[,k]<-obj.aux$N.eff.out
}

#compute MSE and AMSE
MSE.aux<-colMeans(MSE.aux.mat)
AMSE.aux<-mean(MSE.aux)
MSE.aux<-rowMeans(MSE.aux.mat)
N.eff.aux<-rowMeans(N.eff.aux)

#assign mean
m.aux<-obj.aux$m.out


#------------------------------------RESULTS---------------------------------------

#------(sigma, sigma.meas)=(1,3)-------

plot(x,main=paste("Non-Linear Model(",N, "particles)"),xlab="t")
lines(m.pf)
lines(m.reg,col="blue")
lines(m.aux, col="red")
legend(76, -12, legend=c("PF Prior","Regularised","Auxiliary"),
       col=c("black", "blue","red"),lty=c(1,1,1), cex=0.8)


#plot effective sample sie
plot(N.eff.pf,ylab="N.eff",xlab="Time,t", main="Effective Sample Size",type="l")
lines(N.eff.reg,col="blue")
lines(N.eff.aux,col="red")
legend(70, 100, legend=c("PF Prior","Regularised","Auxiliary"),
       col=c("black", "blue","red"),lty=c(1,1,1,1), cex=0.8)

#plot MSE(hardly any difference)
plot(MSE.pf,main="Average MSE",xlab="Time,t",type="l")
lines(MSE.reg,col="blue")
lines(MSE.aux,col="red")

AMSE.pf #24.38401 (1,3)
AMSE.aux #21.86791 (1,3)
AMSE.reg #24.24915 (1,3) 


#------(sigma, sigma.meas)=(3,1)-------

plot(x,main=paste("Non-Linear Model(",N, "particles)"),xlab="t")
lines(m.pf)
lines(m.reg,col="blue")
lines(m.aux, col="red")
lines(m.rej, col="orange")
legend(76, -12, legend=c("PF Prior","Regularised","Auxiliary","Rejection"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)


#plot effective sample sie
plot(N.eff.pf,ylab="N.eff",xlab="Time,t", main="Effective Sample Size",type="l")
lines(N.eff.reg,col="blue")
lines(N.eff.aux,col="red")
legend(70, 450, legend=c("PF Prior","Regularised","Auxiliary"),
       col=c("black", "blue","red"),lty=c(1,1,1,1), cex=0.8)

#plot MSE(hardly any difference)
plot(MSE.pf,main="Average MSE",xlab="Time,t",type="l")
lines(MSE.reg,col="blue")
lines(MSE.aux,col="red")
lines(MSE.rej,col="orange")

AMSE.pf #25.68344 (3,1)  
AMSE.aux #33.08328 (3,1) 
AMSE.reg #25.29179 (3,1) 
AMSE.rej #22.85898(3,1) 





