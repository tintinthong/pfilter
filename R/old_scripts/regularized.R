#--REGULARIZED WORKING FILE--
#\\1)regularized for random walk
#\\2)regularized for population model


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

#run particle filter
obj.randw<-regularized.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,
                           resample.type="standard",N.thr.per=1)

#extract info
m<-obj.randw$m.out
x.pf<-obj.randw$x.pf.out
conf<-obj.randw$conf.out
N.eff<-obj.randw$N.eff.out
res.point<-obj.randw$res.point.out

#main plot
plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
points(y,col="red")
legend(20, 10, legend=c("States, x", "Measurements, y"),
       col=c("black", "red"), pch=c(1,1), cex=0.8)

#plot mean over states
plot(x,main=paste("Random Walk 1D (",N, "particles)"),xlab="t")
lines(m,col="blue")
for(i in 1:N){
  lines(x.pf[,i])
}

#plot credible intervals

#ggplot
#library(ggplot2)
#df <- data.frame(t=1:tau,mu=x,L=conf[,1],U=conf[,2])
#ggplot(df, aes(x = t, y = mu)) +geom_point(size = 1) +geom_errorbar(aes(ymax = U, ymin = L))

#base package
plot(x,main=paste("Random Walk 1D (","N", "particles)"),xlab="t")
for(i in 1:tau){
  segments(i,conf[i,1],i,conf[i,2])
}

#plot Effective
plot(N.eff,type="b")

#plot resampled points
plot(res.point,type="b")

#compute RMSE
M<-50
MSE<-matrix(NA,tau,M)
for(k in 1:M){
  obj.randw<-regularized.filter(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
  MSE[,k]<-obj.randw$MSE.k.out
}
RMSE<-sqrt(colSums(MSE)/tau)
ARMSE<-mean(RMSE)
ARMSE

#plot RMSE
plot(RMSE)



#---------------------------------------POPULATION EXAMPLE------------------------------------------

#set model variables
x0<-0
tau<-100
sigma<-1
sigma.meas<-1
N<-500

#generate states and measurements
set.seed(123)
x<-process.pop(tau=tau,x0=x0,sigma=sigma)
y<-meas.pop(x=x,sigma.meas=sigma.meas)


#run particle filter
obj.pop<-regularized.filter.pop(N=N,x=x,y=y,x0=x0,sigma=sigma,
                             sigma.meas=sigma.meas,resample.type="standard",N.thr.per=0.5)

#extract info
m<-obj.pop$m.out
x.pf<-obj.pop$x.pf.out
conf<-obj.pop$conf.out
N.eff<-obj.pop$N.eff.out
res.point<-obj.pop$res.point.out

#diagnostic plots

#main plot
plot(x,main=paste("Population (",N, "particles)"),xlab="t")
points(y,col="red")
legend(40, -10, legend=c("States, x", "Measurements, y"),
       col=c("black", "red"),pch=c(1,1), cex=0.8)

#plot mean over states
plot(x,main=paste("Population (",N, "particles)"),xlab="t")
lines(m,col="blue")

#plot confidence intervals

#ggplot
#library(ggplot2)
#df <- data.frame(t=1:tau,mu=x,L=conf[,1],U=conf[,2])
#ggplot(df, aes(x = t, y = mu)) +geom_point(size = 1) +geom_errorbar(aes(ymax = U, ymin = L))

#base package
plot(x,main=paste("Population (","N", "particles)"),xlab="t")
for(i in 1:tau){
  segments(i,conf[i,1],i,conf[i,2])
}
plot(density(x.pf[4,]))
hdi(x.pf[4,])

#plot Effective
plot(N.eff,type="b")

#plot resampled points
plot(res.point,type="b")

#compute RMSE
M<-50
MSE<-matrix(NA,tau,M)
for(k in 1:M){
  obj.pop<-regularized.filter.pop(N=N,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
  MSE[,k]<-obj.pop$MSE.k.out
}
RMSE<-sqrt(colSums(MSE)/tau)
ARMSE<-mean(RMSE)
ARMSE


















#------------------------------mess about regularized-----------------------------------------------
# Original distribution is exp(rate = 5)
N = 1000
x <- rexp(N, rate = 5)

par(mfrow=c(1,1))
hist(x, prob = TRUE)
lines(density(x))
lines(density(x,kernel="epanechnikov"),col="red")

# Store the bandwith of the estimated KDE
bw <- density(x)$bw

# Draw from the sample and then from the kernel
means <- sample(x, N, replace = TRUE)
hist(rnorm(N, mean = means, sd = bw), prob = TRUE)


n <- 100000
x <- apply(matrix(runif(3*n, -1, 1), 3), 2, median)
plot(density(x))

xx <- faithful$eruptions
fit <- density(xx)
N <- 1e6
x.new <- rnorm(N, sample(xx, size = N, replace = TRUE), fit$bw)
plot(fit)
lines(density(x.new), col = "blue")
