#-- PLOT PEAKED LIKELIHOOD & PLOT HIST OF WEIGHTS--
# Prior sample in tails of likelihood
# Prior centred at previous particle ; Likelihood centred at measurement
# Illustrate variance of weights(conditional) is high when likelihood is peaked

#set up
sigma<-2
sigma.y<-0.5 
mux<-3
muy<-0
N<-500
wpast<-0.5

#Generate Likelihood and Priors
x<-seq(from=-3,to=10,length.out=1000)
pdf.x<-dnorm(x,mux,sigma)
pdf.y<-dnorm(muy,x,sigma.y)

#Plot Likelihood and Prior
plot(x,pdf.y,type="l",main="Likelihood vs Prior( time= t)",col="red",ylab="probability")
lines(x,pdf.x)
legend(4, 0.6, legend=c("Prior Density", "Likelihood Density"),
       col=c( "black","red"),lty=c(1,1), cex=0.8)

#Plot histogram of weights
rand.x<-rnorm(N,mux,sigma)
like<-dnorm(muy,rand.x,sigma.y)
w<-like*wpast/sum(like*wpast)
my.bin.width<-0.05
hist(w,main=paste("Histogram of weights(",N,"particles)"),xlab="Weights")


