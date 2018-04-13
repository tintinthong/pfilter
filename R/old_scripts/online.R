# Simulate some fake data
set.seed(123)
tau <- 100
x <- cumsum(rnorm(tau))
y <- x + rnorm(tau)
M<-1

# Begin particle filter
N <- 100
x.pf <- matrix(rep(NA,(tau+1)*N),nrow=tau+1)
list.N.thr<-seq(10,500,10)
list.mean<-rep(NA,M)
list.mcerror<-rep(NA,length(list.N.thr))


# 1. Initialize
x.pf[1, ] <- rnorm(N)
m <- rep(NA,tau)
w<-1

for(j in 1:length(list.N.thr)){
  N.thri<-list.N.thr[j]

for(k in 1:M){
  for (t in 2:(tau+1)) {
    
    # 2. Importance sampling step
    x.pf[t, ] <- x.pf[t-1,] + rnorm(N)
    
    #Likelihood
    w.tilde <- dnorm(y[t-1], mean=x.pf[t, ])*w
    
    #Normalize
    w <- w.tilde/sum(w.tilde)
    
    # NOTE: This step isn't part of your description of the algorithm, but I'm going to compute the mean
    m[t-1] <- sum(w*x.pf[t,])
    N.eff<-1/sum(w^2)
    
    # 3. Resampling step
    if(N.eff<N.thri){
      s <- sample(1:N, size=N, replace=TRUE, prob=w)
      x.pf <- x.pf[, s]
      w<-rep(1/N,N)
    }
    
    
    # Note: resample WHOLE path, not just x.pf[t, ]
  }

  list.mean[k]<-m[tau]
  
}

list.mcerror[j]<-var(list.mean)
}

plot(list.mcerror)

plot(x)
for(i in 1:N){
  lines(x.pf[,i])
}
