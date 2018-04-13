norm<-function (n, alpha) 
{
  vec <- vector("numeric", n)
  x <- 0
  vec[1] <- x
  for (i in 2:n) {
    innov <- runif(1, -alpha, alpha)
    can <- x + innov
    aprob <- min(1, dnorm(can)/dnorm(x))
    u <- runif(1)
    if (u < aprob) 
      x <- can
    vec[i] <- x
  }
  vec
}


normvec<-norm(10000,1)
par(mfrow=c(2,1))
plot(ts(normvec))
hist(normvec,30)
par(mfrow=c(1,1))

gamm<-function (n, a, b) 
{
  mu <- a/b
  sig <- sqrt(a/(b * b))
  vec <- vector("numeric", n)
  x <- a/b
  vec[1] <- x
  for (i in 2:n) {
    can <- rnorm(1, mu, sig)
    aprob <- min(1, (dgamma(can, a, b)/dgamma(x, 
                                              a, b))/(dnorm(can, mu, sig)/dnorm(x, 
                                                                                mu, sig)))
    u <- runif(1)
    if (u < aprob) 
      x <- can
    vec[i] <- x
  }
  vec
}


vec<-gamm(10000,2.3,2.7)
par(mfrow=c(2,1))
plot(ts(vec))
hist(vec,30)
par(mfrow=c(1,1))

#met hastings in uscented particle filter link in browser
#smooths particles 

for(i in 1:N){
  for(t in 2:tau){
    u<-runif(1,0,1)
    xi<-rnorm(1,x.pf[t-1,i],sigma)
    if(xi<min(1,dnorm(y[t],xi,sigma.meas)/dnorm(y[t],x.pf[t-1,i],sigma.meas))){
      x.pf[t-1,i]<-xi
    }
  }
}

plot(x)
for(i in 1:tau){
  lines(x.pf[,i])
}


