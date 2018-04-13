
#\\all of this is wrong. Asymptotic variance is calculated conditional on previous draws. It can be shown 
#in creal and doucet and johansen that aysmpt


#--Compare performance for different N particles--
#\\ bootstrap for random walk
#\\2)bootstrap for population model

#\\always set working directory to file location
#\\ modules meant to be sourced ie functions are dependent on source("module.R")- Do not attempt to use just one function 

#load library
source("filter_mod.R")

M<-100
N<-100
test<-rep(NA,length(x))
test.norm<-rep(NA,length(x))
test.wi<-rep(NA,length(M))
for(k in 1:M){
  x<-rnorm(N,5)
  w<-dnorm(x,3,1,log=TRUE)
  w<-w/sum(w)
  test.norm[k]<-mean(x*w)
  test.wi[k]<-w[1]
  test[k]<-mean((1/N)*sample(x,length(x),prob=w,replace=TRUE))
}

sd(test.wi) #0.005274663(sd.y=1) 0.0006609808(sd.y=4)#\\weights of the ith particle degnerate when the likelihood becomes peaked
sd(test)  # 0.1543296(sd.y=1) 0.1248949(sd.y=4)#\\more degeneracy the more monte carlo error in resampled estimates
sd(test) #0.1534726(N=100) 0.04931775(N=1000)#\\more particles less mc error
sd(test.norm) #0.09766345(N=100) 0.03193381(N=1000) #\\more particles less mc error
mean((test.norm-5)^2)#2445.185(sd.y=1) 2449.725(sd.y=4)\\more degneracy, mse can either go up or down
sd(test.norm)# 0.001105887(sd.y=1) 0.0009983367(sd.y=4) \\more degneracy, higher mc error

sd(test.norm) #0.001046114
sd(test) #0.1436265
#\\mc error increases with resampling
#\\resampled estimates have higher variance

mean(test.norm)-5
mean(test)-5
#\\resampled estimates have maybe slightly higher bias
#\\normal to have a bias because these measurements not generated around true mean


mean((test.norm-5)^2) #\\MSE
(mean(test.norm)-5)#\\bias
var(test.norm)
(mean(test.norm)-5)^2+var(test.norm)


mean(test-5)^2 #\\MSE  
(mean(test)-5)#\\bias
var(test)
(mean(test)-5)^2+var(test) 


#MSE of resampled is lower

#\\questions: how is the bias computed?
#\\why is the MSE lower?





#\\M=100, N=100, sd.y=1 is default




#\\mse is a measure of bias and variance

#\\note: weights do not degenerate more when the mean of the likelihood is not the same 
#\\as the mean of th proposal

#---------------------------------------RANDOM WALK EXAMPLE------------------------------------------

#set model variables
x0<-0
tau<-100
sigma<-1
sigma.meas<-0.5
N<-100

#generate states and measurements
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)


#generate simulated data y
y.mat<-matrix(NA,nrow=tau,ncol=M)
for(k in 1:M){
  y.mat[,k]<-rand.y.1D(x,sigma.meas=sigma.meas)
}



#compute variance for each N
#\\only compute variance for filtered estimate of particular time t 
N.thr.list<-seq(0,1,0.1)
M<-100
MC.error.list<-rep(NA,length(N.thr.list) )
mean.list<-rep(NA,M)
#try for different amount of resamples
for(i in 1:length(N.thr.list)){
  N.thri<-N.thr.list[i]
  #compute mc error
  for(k in 1:M){
    #y<-rand.y.1D(x,sigma.meas=sigma.meas)
    obj.randw<-particle.filter.path(N=N,x=x,y=y.mat[,k],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=N.thri)
    x.pf<-obj.randw$x.pf.out
    mean.list[k]<-mean(x.pf[tau,])
  }
  MC.error.list[i]<-sd(mean.list)
}
print(mean.list)

plot(N.thr.list,MC.error.list,main=paste("MC.error(", M ,"MC runs)"),xlab="N.thr",ylab="MC.error")


#\\more we resample, the higher the monte carlo error
#\\after a particular threshold range, you'' find that the mc error actually decreases the monte you resample
#\\ the sample does not impoverish 

#\\degenerate case does not gradually impoversih because after resampling, they are dispersed 
#through propoagation step. They only face one impoverishment problem






















MC.error.list[1]




#compute variance for each N
#\\only compute variance for filtered estimate of particular time t 
N.thr.list<-seq(0,1,0.1)
M<-100
MC.error.list<-rep(NA,length(N.thr.list) )
mean.list<-matrix(NA,nrow=tau,ncol=M)
#try for different amount of resamples
for(i in 1:length(N.thr.list)){
  N.thri<-N.thr.list[i]
  #compute mc error
  for(k in 1:M){
    obj.randw<-particle.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=N.thri)
    x.pf<-obj.randw$x.pf.out
    mean.list[,k]<-rowMeans(x.pf)
  }
  MC.error.list[i]<-sd(mean.list)
}

plot(N.thr.list,MC.error.list,main=paste("MC.error(", M ,"MC runs)"),xlab="N.thr",ylab="MC.error")



length(rowMeans(x.pf))
length(mean.list[,1])



N.thr.list<-seq(0,1,by=0.1)
var.list<-rep(NA,M)
MC.error.list<-rep(NA,length(N.thr.list) )
#try for different amount of resamples
for(i in 1:length(N.thr.list)){
  N.thri<-N.thr.list[i]
  #compute mc error
  obj.randw<-particle.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=N.thri)
  x.pf<-obj.randw$x.pf.out
  MC.error.list[i]<-sd(x.pf[8,])/sqrt(N)
}

plot(N.thr.list,MC.error.list,main=paste("MC.error(", M ,"MC runs)"),xlab="N.thr",ylab="MC.error")





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



#mc error
M<-500
x<-rep(NA,M)
for(i in 1:M){
  N<-10
  obj<-replicate(i,rnorm(N))
  x[i]<-sd(colMeans(obj))
}
plot(x)
1/sqrt(10)
#1/N is the variance

N<-1000
x<-rep(NA,N)
for(i in 1:N){
  x[i]<-var(rnorm(i))
}



