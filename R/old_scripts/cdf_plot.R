#--PLOT EMPIRICAL CDF VS TRUE CDF--

#set up 
mu<-0
sigma<-1
N<-500 #50 particles or 500 particles
x = seq(-3,3,.1)
x.p<-rnorm(N,mu,sigma)

#plot cdf
cdf<-pnorm(x,mu,sigma)
plot(x,cdf,type="l",main=paste("True Cdf vs Empricial Cdf","(",N,"Particles)"),ylab="Cdf",xlab="x",col="red")

#plot ecdf
breaks <-seq(min(x.p), max(x.p), length.out=N)
count<-cut(x.p, breaks, right=FALSE)
freq<-table(count)
cum.freq<-cumsum(freq)
cum.freq<-c(0,cum.freq)
cum.prob<-cum.freq/N
points(cum.prob~breaks, type = "s")
