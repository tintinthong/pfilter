#particles

#set up 
mu<-2
sigma<-5
N<-200
x = seq(-5,15,.1)
pdf<-dnorm(x,mu,sigma)
x.p<-rnorm(N,mu,sigma)

#Plotting the true density as a sanity check
truth = .3*dnorm(x,0,1) + .5*dnorm(x,4,1)
truth2=.3*dnorm(x.p,0,1) + .5*dnorm(x.p,4,1)
#plot(density(rand.samples),main="Density Estimate of the Mixture Model",ylim=c(0,.2),lwd=2)


plot(x,truth,col="red",type="l",main="Filtered Density(200 particles)",ylab="Density")
segments(x0=x.p, y0=rep(0,N), y1=0.1*rep(0.2,N)) 
lines(x,pdf)
legend(6,0.15, legend=c("Proposal", "True posterior"),
       col=c( "black","red"),lty=c(1,1), cex=0.8)


#reweight 
w<-truth2/dnorm(x.p,mu,sigma)
plot(x,truth,col="red",type="l",main="Filtered Density(200 particles)",ylab="Density")
lines(x,pdf)
segments(x0=x.p, y0=rep(0,N), y1=w*0.3*rep(0.2,N)) 
#legend(6,0.15, legend=c("Proposal", "True posterior"),
#       col=c( "black","red"),lty=c(1,1), cex=0.8)



#set up 
mu<-2
sigma<-5
N<-50
x = seq(-5,15,.1)
#pdf<-dnorm(x,mu,sigma)
x.p<-rnorm(N,mu,sigma)


#plot cdf
cdf<-pnorm(x,mu,sigma)
plot(x,cdf,type="l",main=paste("True Cdf vs Empricial Cdf","(",N,"Particles)"),ylab="Cdf",xlab="x",col="red")

#pdf.p<-dnorm(sort(x.p),mu,sigma)

#plot ecdf
breaks <-seq(min(x.p), max(x.p), length.out=N)
count<-cut(x.p, breaks, right=FALSE)
freq<-table(count)
cum.freq<-cumsum(freq)
#cbind(cum.freq)
cum.freq<-c(0,cum.freq)
cum.prob<-cum.freq/N
points(cum.prob~breaks, type = "s")




