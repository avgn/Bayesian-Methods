####Gibbs sampler for Galaxy data, with finite k and non updating (fixed) Variances ####


###Set Up###


##Open data and tools, plot data, which we will call X.
library(MASS)
library(label.switching)
gal<-galaxies
n<-length(gal)
hist(gal, breaks=40, freq=FALSE,xlab="Galaxy Speed", main="Histogram of Galaxy Data")


##Set non-updating parameters
k<-7                     #Number of normal kernels
a<-1                     #Sum of Dirichlet distribution concentration parameters, assumed symmetric and therefore all equal to a/k.
sig<-(0.8/k)*sd(gal)     #Standard Deviation of the normal kernels. NB: Notice that the variance of kernels is not updating.
sigstr<-(0.8/k)*sd(gal)  #Standard Deviation of Prior (normal) distribution on the means of the kernels, which we will call "centroids". NB: Notice that the variance of the prior on the centroids is not updating 
center<-mean(gal)        #Mean of the Prior (normal) distribution on the centroids
M<-10000                 #Number of iteratons of the Gibb Sampler


##Creation of storage space for updating parameters
PI<-array(0, c(M,k))        #For the k possible j, Weight of the kernel j in joint mixed distribution of the data.
PiXandZ<-array(0, c(M,k))   #For the k possible j(i), Joint probability of label Z(i) belonging to the kernel j(i) and data X(i) being extracted, given PI. PiXandZ is equal to probability of Z(i)belonging to kernel j(i) given PI times probability of data X(i)being extracted given Z(i) belonging to J(i) (and given PI).
PiX<-array(0, c(M,1))       #Marginal probability of data X(i) being extracted given PI, which is the sum over all kernels j of PiXandZ
PiZ<-array(0, c(M,k))       #For the k possible j(i), Posterior Probability of observation (i) belonging to the kernel J(i). Equal to the ratio of PiXandZ to PiX.
Z<-array(0, c(M,n))         #Labels denoting the kernel to whom the observation belongs. Extracted from multinomial with probabilities PiZ, or initially PI.
Nj<-array(0, c(M,k))        #Numbers of observations in kernels
theta<-array(0, c(M, k))    #Centroids of the k normal kernels. Each is extracted from a normal whose mean is the average of the prior theta and of the updating theta, weighted by the respective precision. The variance of the normal is the inverse of the sum of the precisions.


###Gibbs Sampler###


##Generate initial values from Priors

# Generate initial values PI for weights of kernels from Dirichlet prior
y<-rgamma(k,a/k,1)                                   
PI[1, ]<-y/sum(y)

# Generate initial values Z for labels from multinomial distribution given the initial PI weights
for (i in 1:n) {
  Z[1, i]<-which(rmultinom(1, 1, PI[1,])==1)
}

# Compute initial values Nj counting numbers of observations' labels Z in the kernels
for (j in 1:k) {
  Nj[1, j]<-sum(Z[1,]==j)
}

#Generate initial values THETA of centroids of kernels from normal prior.
for (j in 1:k) {
  theta[1,j]<-rnorm(1, mean(gal), sigstr)              
}


## Begin Iterative Updating
for (m in 2:M) {

#Generate Posterior Weights PI from Dirichlet Distribution whose concentration parameters are updated by Nj
y<-rep(0,k)
for(j in 1:k) {
  y[j]<-rgamma(1,(a/k)+Nj[m-1,j],1)
}
PI[m,]<-y/sum(y)
  
#Generate Posterior Labels Z from full conditional of Z(i) given all other Z, all X, all PI and all THETA, which is Multinomial with Posterior weights PiZ. 

for(i in 1:n) {
  for(j in 1:k) {
    PiXandZ[m,j]<-(PI[m,j]*dnorm(gal[i],theta[m-1,j],sig))
  }
  PiX[m]<-sum(PiXandZ[m,])  
  PiZ[m,]<-PiXandZ[m,]/PiX[m]
  Z[m,i]<-which(rmultinom(1,1,PiZ[m,])==1)
}
 
# Compute Nj from the Z
for(j in 1:k) {
  Nj[m,j]<-length(which(Z[m,]==j))
} 

#Generate Posterior Centroids
for(j in 1:k) {
  mean_Nj<-sum(gal[Z[m,]==j])/Nj[m,j]
  if (is.na(mean_Nj)) {
    theta[m,j]<-rnorm(1, center, sigstr)
  } else {
    theta[m, j]<-rnorm(1,((mean_Nj*(sig^2/Nj[m,j])^(-1)+theta[m-1,j]*sigstr^(-2))/((sig^2/Nj[m,j])^(-1)+sigstr^(-2))),((sig^2/Nj[m,j])^(-1)+sigstr^(-2))^(-1/2))
  }
}
  
# end Gibbs
}



### Label Switching Adjustment###

#Generate label reordering array
ls <- dataBased(gal, k, Z)
LS <- array(unlist(ls), dim=c(M,k))

#Generate storage space for sorted output
theta.sor <- array(0, dim=c(M,k))
PI.sor <- array(0, dim=c(M,k))
PiZ.sor <- array(0, dim=c(M,k))
Nj.sor <- array(0, dim=c(M,k))

#Reorder results according to reordering array
for (m in 1:M) {
  for (j in 1:k) {
    theta.sor[m,j] <- theta[m,LS[m,j]]
    PI.sor[m,j] <- PI[m,LS[m,j]]
    PiZ.sor[m,j] <- PiZ[m,LS[m,j]]
    Nj.sor[m,j] <- Nj[m,LS[m,j]]
  }
}

##Convergence, with and without LS adjustment##

# plot of convergence in the first iterations, before label switching correction##
windows()
par(mfrow=c(2,1))
matplot(theta[1:300,], type="l", main="Convergence of Thetas", ylab="theta")
matplot(PI[1:300,], type="l", main="Convergence of PIs ", ylab="pi")

# plot of convergence in the first iterations, after label switching correction##
windows()
par(mfrow=c(2,1))
matplot(theta.sor[1:300,], type="l", main="Convergence of Thetas sorted with LS", ylab="theta")
matplot(PI.sor[1:300,], type="l", main="Convergence of PIs sorted with LS", ylab="pi")

# plot of last iterations before label switching correction
windows()
par(mfrow=c(2,1))
matplot(theta[(M-1000):M,], type="l", main="Ergodic Thetas", ylab="theta")
matplot(PI[(M-1000):M,], type="l", main="Ergodic PIs", ylab="pi")

# plot of last iterations after label switching correction
windows()
par(mfrow=c(2,1))
matplot(theta.sor[(M-1000):M,], type="l", main="Ergodic Thetas sorted with LS", ylab="theta")
matplot(PI.sor[(M-1000):M,], type="l", main="Ergodic PIs sorted with LS", ylab="pi")

##Means, with and without LS adjustment##

# means of theta with LS after burnt-in (2000 obs)
theta.sor.means<-rep(0,k)
for (j in 1:k) {
  theta.sor.means[j]<-mean(theta.sor[2000:M,j])
}
theta.sor.means

# means of PI with LS after burnt-in (2000 obs)
PI.sor.means<-rep(0,k)
for (j in 1:k) {
  PI.sor.means[j]<-mean(PI.sor[2000:M,j])
}
PI.sor.means

# means of theta without LS after burnt-in (2000 obs)
theta.means<-rep(0,k)
for (j in 1:k) {
  theta.means[j]<-mean(theta[2000:M,j])
}
theta.means

# means of PI without LS after burnt-in (2000 obs)
PI.means<-rep(0,k)
for (j in 1:k) {
  PI.means[j]<-mean(PI[2000:M,j])
}
PI.means

## Likelihood, AIC and BIC
Lik <- rep(0,82)
for (j in 1:k) {
  Lik <- Lik+PI.sor.means[j]*dnorm(gal,theta.sor.means[j],sig)
}
logLik <- sum(log(Lik))
AIC <- -2*logLik+2*(2*k)
BIC <- -2*logLik+(2*k)*log(n)
HQC <- -2*logLik+2*(2*k)*log(log(n))

## Plot comparing data and fitted distribution with LS adjustment##

# With LS adjustment
windows()
x.axis<-seq(min(gal),max(gal),1)
hx<-PI.sor.means[1]*dnorm(x.axis,theta.sor.means[1],sig)
for (j in 2:k) {
  hx<-hx+PI.sor.means[j]*dnorm(x.axis,theta.sor.means[j],sig)
}
hist(gal, breaks=40, freq=FALSE, xlab="Galaxy Speed", ylim=c(0,max(hx)), main="Galaxy Data vs Fitted Mixture Distribution, with LS adj.")
lines(x.axis,hx,col=4)

# Without LS adjustment
windows()
x.axis<-seq(min(gal),max(gal),1)
hx<-PI.means[1]*dnorm(x.axis,theta.means[1],sig)
for (j in 2:k) {
  hx<-hx+PI.means[j]*dnorm(x.axis,theta.means[j],sig)
}
hist(gal, breaks=40, freq=FALSE, xlab="Galaxy Speed", ylim=c(0,max(hx)), main="Galaxy Data vs Fitted Mixture Distribution, without LS adj.")
lines(x.axis,hx,col=4)


### Detail of Output Without Label Switching Adjustment###

#Set filtering farameter
f<-250
Filtf <- rep(1/(f+1), (f+1))

## Theta

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.theta<-rep(0,k)
for(j in 1:k) {
  plot(theta[,j],ylim=c(min(theta),max(theta)),xlab="Iteration")
  mean.theta[j]<-mean(theta[,j])
}
plot(mean.theta,ylim=c(min(theta),max(theta)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
theta.filt<-array(0, c(M,k))
for(j in 1:k) {
  theta.filt[,j]<-filter(theta[,j],Filtf)
  plot(theta.filt[,j],ylim=c(min(theta),max(theta)),xlab="Iteration")
}


## PI

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.PI<-rep(0,k)
for(j in 1:k) {
  plot(PI[,j],ylim=c(min(PI),max(PI)),xlab="Iteration")
  mean.PI[j]<-mean(PI[,j])
}
plot(mean.PI,ylim=c(min(PI),max(PI)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
PI.filt<-array(0, c(M,k))
for(j in 1:k) {
  PI.filt[,j]<-filter(PI[,j],Filtf)
  plot(PI.filt[,j],ylim=c(min(PI),max(PI)),xlab="Iteration")
}


## PiZ

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.PiZ<-rep(0,k)
for(j in 1:k) {
  plot(PiZ[,j],ylim=c(min(PiZ),max(PiZ)),xlab="Iteration")
  mean.PiZ[j]<-mean(PiZ[,j])
}
plot(mean.PiZ,ylim=c(min(PiZ),max(PiZ)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
PiZ.filt<-array(0, c(M,k))
for(j in 1:k) {
  PiZ.filt[,j]<-filter(PiZ[,j],Filtf)
  plot(PiZ.filt[,j],ylim=c(min(PiZ),max(PiZ)),xlab="Iteration")
}

## Nj

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.Nj<-rep(0,k)
for(j in 1:k) {
  plot(Nj[,j],ylim=c(min(Nj),max(Nj)),xlab="Iteration")
  mean.Nj[j]<-mean(Nj[,j])
}
plot(mean.Nj,ylim=c(min(Nj),max(Nj)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
Nj.filt<-array(0, c(M,k))
for(j in 1:k) {
  Nj.filt[,j]<-filter(Nj[,j],Filtf)
  plot(Nj.filt[,j],ylim=c(min(Nj),max(Nj)),xlab="Iteration")
}



###Detail of Output With Label Switching Adjustment###


#Set filtering farameter
f<-250
Filtf<-rep(1/(f+1), (f+1))


## Theta

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.theta.sor<-rep(0,k)
for(j in 1:k) {
  plot(theta.sor[,j],ylim=c(min(theta.sor),max(theta.sor)),xlab="Iteration")
  mean.theta.sor[j]<-mean(theta.sor[,j])
}
plot(mean.theta.sor,ylim=c(min(theta.sor),max(theta.sor)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
theta.sor.filt<-array(0, c(M,k))
for(j in 1:k) {
  theta.sor.filt[,j]<-filter(theta.sor[,j],Filtf)
  plot(theta.sor.filt[,j],ylim=c(min(theta.sor),max(theta.sor)),xlab="Iteration")
}


## PI

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.PI.sor<-rep(0,k)
for(j in 1:k) {
  plot(PI.sor[,j],ylim=c(min(PI.sor),max(PI.sor)),xlab="Iteration")
  mean.PI.sor[j]<-mean(PI.sor[,j])
}
plot(mean.PI.sor,ylim=c(min(PI.sor),max(PI.sor)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
PI.sor.filt<-array(0, c(M,k))
for(j in 1:k) {
  PI.sor.filt[,j]<-filter(PI.sor[,j],Filtf)
  plot(PI.sor.filt[,j],ylim=c(min(PI.sor),max(PI.sor)),xlab="Iteration")
}


## PiZ

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.PiZ.sor<-rep(0,k)
for(j in 1:k) {
  plot(PiZ.sor[,j],ylim=c(min(PiZ.sor),max(PiZ.sor)),xlab="Iteration")
  mean.PiZ.sor[j]<-mean(PiZ.sor[,j])
}
plot(mean.PiZ.sor,ylim=c(min(PiZ.sor),max(PiZ.sor)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
PiZ.sor.filt<-array(0, c(M,k))
for(j in 1:k) {
  PiZ.sor.filt[,j]<-filter(PiZ.sor[,j],Filtf)
  plot(PiZ.sor.filt[,j],ylim=c(min(PiZ.sor),max(PiZ.sor)),xlab="Iteration")
}

## Nj

#Unfiltered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
mean.Nj.sor<-rep(0,k)
for(j in 1:k) {
  plot(Nj.sor[,j],ylim=c(min(Nj.sor),max(Nj.sor)),xlab="Iteration")
  mean.Nj.sor[j]<-mean(Nj.sor[,j])
}
plot(mean.Nj.sor,ylim=c(min(Nj.sor),max(Nj.sor)))

#Filtered
windows()
par(mfrow=c(2,ceiling((k+1)/2)))
Nj.sor.filt<-array(0, c(M,k))
for(j in 1:k) {
  Nj.sor.filt[,j]<-filter(Nj.sor[,j],Filtf)
  plot(Nj.sor.filt[,j],ylim=c(min(Nj.sor),max(Nj.sor)),xlab="Iteration")
}

