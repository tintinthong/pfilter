
#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-100
sigma<-1
sigma.meas<-3
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
#compute RMSE
MSE.pf<-matrix(NA,tau,M)
for(k in 1:M){
  obj.pf<-particle.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,
                              sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)
  MSE.pf[,k]<-obj.pf$MSE.k.out
  print(k)
}
RMSE.pf<-sqrt(colMeans(MSE.pf))
ARMSE.pf<-mean(RMSE.pf)
ARMSE.pf #4.638(100)

#mc error
sd(RMSE.pf)

#compute MSE 
MSE.pf<-rowMeans(MSE.pf)
plot(MSE.pf)

#assign mean
m.pf<-obj.pf$m.out

plot(x)
lines(m.pf)
#------------------------------------Extended---------------------------------------
#\\does not work for large sigma.meas 
#\\use optimum proposal

#generate simulated data y
set.seed(123)
x<-process.pop(tau=tau,x0=x0,sigma=sigma)
y.mat.ext<-matrix(NA,nrow=tau,ncol=M)
for(k in 1:M){
  y.mat.ext[,k]<-extended.meas.1D(x0,x,sigma.meas=sigma.meas)
}

#compute RMSE
MSE.ext<-matrix(NA,tau,M)
for(k in 1:M){
  obj.ext<-extended.filter(N=N,x=x,y=y.mat.ext[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard")
  MSE.ext[,k]<-obj.ext$MSE.k.out
}
RMSE.ext<-sqrt(colMeans(MSE.ext))
ARMSE.ext<-mean(RMSE.ext)
ARMSE.ext 

#mc error
sd(RMSE.ext)

#compute MSE 
MSE.ext<-rowMeans(MSE.ext)
plot(MSE.ext)

#assign mean
m.ext<-obj.ext$m.out

#comparison of the approximation plots
plot(x)
points(y.mat.ext[,1],col="red")
points(y.mat[,1],col="blue ")


#-----------------------------------Regularized---------------------------------------

#compute RMSE
MSE.reg<-matrix(NA,tau,M)
for(k in 1:M){
  obj.reg<-regularized.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=0.5)
  MSE.reg[,k]<-obj.reg$MSE.k.out
}
RMSE.reg<-sqrt(colMeans(MSE.reg))
ARMSE.reg<-mean(RMSE.reg)
ARMSE.reg #  4.914272(100)

#mc error
sd(RMSE.reg) 

#compute MSE 
MSE.reg<-rowMeans(MSE.reg)
plot(MSE.reg)

#assign mean
m.reg<-obj.reg$m.out

plot(x)
lines(m.reg)

#------------------------------------Rejection---------------------------------------


#compute RMSE
MSE.rej<-matrix(NA,tau,M)
for(k in 1:M){
  obj.rej<-rejection.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas)
  MSE.rej[,k]<-obj.rej$MSE.k.out
  print(paste("k",k))
}
RMSE.rej<-sqrt(colMeans(MSE.rej))
ARMSE.rej<-mean(RMSE.rej)
ARMSE.rej  #took a long time:4.733422 
RMSE.rej[1] #5.236246
hist(RMSE.rej)

#mc error
sd(RMSE.rej) #0.3220664

#compute MSE 
MSE.rej<-rowMeans(MSE.rej)
plot(MSE.rej)

#assign mean
m.rej<-obj.rej$m.out

plot(x)
lines(m.rej)


#------------------------------------Auxiliary---------------------------------------


#compute RMSE
MSE.aux<-matrix(NA,tau,M)
for(k in 1:M){
  obj.aux<-auxiliary.filter.pop(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas)
  MSE.aux[,k]<-obj.aux$MSE.k.out
}
RMSE.aux<-sqrt(colMeans(MSE.aux))
ARMSE.aux<-mean(RMSE.aux)
ARMSE.aux 

#mc error
sd(RMSE.aux)

#compute MSE 
MSE.aux<-rowMeans(MSE.aux)
plot(MSE.aux)

#assign mean
m.aux<-obj.aux$m.out




#------------------------------------Final Plot---------------------------------------

#not yet done yet

plot(x,main=paste("Random Walk (",N, "particles)"),xlab="t")
lines(m.pf)
lines(m.reg,col="blue")
lines(m.ext, col="red")
lines(m.exact, col="orange")
legend(20, 10, legend=c("PF Prior", "Extended","Rejection","Auxiliary"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)

#\\exact performs the worse



