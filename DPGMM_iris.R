library(datasets)
library(DPpackage)
library(ggplot2)
library(gridExtra)

# load data
data(iris)
attach(iris)

# visualize data
hist(Petal.Length, breaks = 30, freq = FALSE)
hist(Petal.Width, breaks = 30, freq = FALSE)

qplot(Petal.Length, Petal.Width, ylab = "Petal Width", main = "Petal Length vs. Petal Width in Fisher's Iris data")
qplot(Petal.Length, Petal.Width, color = Species, xlab = "Petal Length", ylab = "Petal Width", main = "Petal Length vs. Petal Width in Fisher's Iris data")

mean(Petal.Length)
mean(Petal.Width)
var(Petal.Length)
var(Petal.Width)

# set parameters
m1 <- c(7, 1)
s1 <- matrix(c(0.3,0,0,0.28),ncol=2)

prior <- list(alpha=1,m1=m1,k0=1/5,nu1=4,psiinv1=s1)

state <- NULL

nburn <- 5000
nsave <- 10000
nskip <- 10
ndisplay <- 1000
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

# fit the model
fit1 <- DPdensity(y=cbind(Petal.Length, Petal.Width),prior=prior,mcmc=mcmc,state=state,status=TRUE)
clusters <- fit1$state$ss

fit1
# plot clusters' density
plot(fit1, output = "param", ask=FALSE)
clusters
DPrandom(fit1)
# plot estimated parameters' densities
plot(DPrandom(fit1,predictive=TRUE),ask=FALSE)

# plot clustered data
qplot(Petal.Length, Petal.Width, color = factor(clusters), xlab = "Petal Length", ylab = "Petal Width", main = "DPGMM Clustering")

# plot clustered and real data
p1 <- qplot(Petal.Length, Petal.Width, color = factor(clusters), xlab = "Petal Length", ylab = "Petal Width", main = "DPGMM Clustering")
p2 <- qplot(Petal.Length, Petal.Width, color = Species, xlab = "Petal Length", ylab = "Petal Width", main = "Petal Length vs. Petal Width in Fisher's Iris data")
grid.arrange(p1, p2, ncol = 2)
