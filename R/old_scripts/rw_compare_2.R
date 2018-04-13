
#load library
source("filter_mod.R")

#set model variables
x0<-0
tau<-100
sigma<-1
sigma.meas<-0.5
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
#compute RMSE
MSE.boot<-matrix(NA,tau,M)
for(k in 1:M){
  obj.boot<-particle.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
  MSE.boot[,k]<-obj.boot$MSE.k.out
  print(k)
}
RMSE.boot<-sqrt(colMeans(MSE.boot))
ARMSE.boot<-mean(RMSE.boot)
ARMSE.boot  #0.4794538

#mc error
sd(RMSE.boot) #0.03934531(100) #0.03517375(1000)

#compute MSE 
MSE.boot<-rowMeans(MSE.boot)
plot(MSE.boot)

#assign mean
m.boot<-obj.boot$m.out

#------------------------------------Exact---------------------------------------

#compute RMSE
MSE.exact<-matrix(NA,tau,M)
for(k in 1:M){
  obj.exact<-exact.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard")
  MSE.exact[,k]<-obj.exact$MSE.k.out
}
RMSE.exact<-sqrt(colMeans(MSE.exact))
ARMSE.exact<-mean(RMSE.exact)
ARMSE.exact #0.4675411(100)

#mc error
sd(RMSE.exact)

#compute MSE 
MSE.exact<-rowMeans(MSE.exact)
plot(MSE.exact)

#assign mean
m.exact<-obj.exact$m.out



#------------------------------------Optimum---------------------------------------

#compute RMSE
MSE.opt<-matrix(NA,tau,M)
for(k in 1:M){
  obj.opt<-optimum.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr=0.5)
  MSE.opt[,k]<-obj.opt$MSE.k.out
}
RMSE.opt<-sqrt(colMeans(MSE.opt))
ARMSE.opt<-mean(RMSE.opt)
ARMSE.opt 

#mc error
sd(RMSE.opt)

#compute MSE 
MSE.opt<-rowMeans(MSE.opt)
plot(MSE.opt)

#assign mean
m.opt<-obj.opt$m.out


#------------------------------------Likelihood---------------------------------------

#compute RMSE
MSE.like<-matrix(NA,tau,M)
for(k in 1:M){
  obj.like<-likelihood.filter(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard")
  MSE.like[,k]<-obj.like$MSE.k.out
}
RMSE.like<-sqrt(colMeans(MSE.like))
ARMSE.like<-mean(RMSE.like)
ARMSE.like 

#mc error
sd(RMSE.like)

#compute MSE 
MSE.like<-rowMeans(MSE.like)
plot(MSE.like)

#assign mean
m.like<-obj.like$m.out

#------------------------------------Final Plot---------------------------------------
plot(x,main=paste("Random Walk (",N, "particles)"),xlab="t")
lines(m.boot)
lines(m.opt,col="blue")
lines(m.like, col="red")
lines(m.exact, col="orange")
legend(20, 10, legend=c("Bootstrap", "PF Optimum","Likelihood","Exact"),
       col=c("black", "blue","red","orange"),lty=c(1,1,1,1), cex=0.8)

#\\exact performs the worse



