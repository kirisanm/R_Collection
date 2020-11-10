rm(list = ls())
set.seed(1)
# This is the R code soultion for 2020 STA510 mandatory assignment 2
#
prob1 <- FALSE
prob2 <- FALSE
prob3 <- TRUE
prob4 <- FALSE

if(prob1){
#
# 1,c)
#
# simulate many replica of the monte carlo estimator
Nsim <- 10000
n <-1000
theta <- numeric(Nsim)
for(i in 1:Nsim){
  X <- runif(n,min=-1.0,max=1.0)
  Y <- runif(n,min=-1.0,max=1.0)
  theta[i] <- 4.0*mean(X^2+Y^2<=1.0)
}

# check against N(pi,pi*(4-pi)/n) density.
# Formally, the distribution of the Monte Carlo estimator is a shifted Binomial
# distribution, but since we are working with a sum of iid random variables, the
# central limit theorem tells us that the estimator will be well-approximated 
# by a normal distribution.

hist(theta,probability = TRUE)
xg <- seq(from=min(theta),to=max(theta),length.out = 1000)
lines(xg,dnorm(xg,mean=pi,sd=sqrt(pi*(4-pi)/n)))
# looks good...

#
# 1,d)
#


# the sought probability translates to P(3.14000000<=theta<3.15000000)
print(paste0("probability: ",mean(theta>=3.14 & theta<3.15)))
# alternatively, by converting to strings, don't worry if you 
# don't see what's going on here:
print(paste0("probability: ",mean(grepl("3.14",format(x=theta,digits=6)))))


#
# 1,f)
#
# We do not expect any improvement as the g-function lacks the monotone 
# property. In fact g(X,Y) = g(-X,-Y), and therefore introducing these
# antithetic variables simply amounts to reducing the effective sample
# size by a factor 0.5.
# Check by simulating regular antithetic variables:

Nsim <- 10000
n <-1000
thetaA <- numeric(Nsim)
for(i in 1:Nsim){
  X1 <- runif(n/2,min=-1.0,max=1.0)
  X <- c(X1,-X1)
  Y1 <- runif(n/2,min=-1.0,max=1.0)
  Y <- c(Y1,-Y1)
  thetaA[i] <- 4.0*mean(X^2+Y^2<=1.0)
}

print(paste0("SD of antithetic MC integral, n = ",n," : ",sd(thetaA)))
print(paste0("Theoretical SD of regular MC integral, n = ",n,", : ",sqrt(pi*(4-pi)/n)))
print(paste0("Theoretical SD of regular MC integral, n = ",n/2,", : ",sqrt(pi*(4-pi)/(0.5*n))))
# Indeed, introducing regular antithetics here reduced efficiency corresponding to reducing n
# by a factor 0.5.


#
#  1,g)
#


# from assignment
shift <- function(u){return(((u+2.0)%%2.0)-1.0)}
# obtain integration weights 
X1 <- X1 <- runif(n/2,min=-1.0,max=1.0)
X <- c(X1,shift(X1))
Y1 <- runif(n/2,min=-1.0,max=1.0)
Y <- c(Y1,shift(Y1))
wts <- 4.0*as.numeric(X^2+Y^2<=1.0)
thetaA2 <- mean(wts)
print(paste0("Shift antithetics; estimate: ",thetaA2)) # not too bad 
print(paste0("correlation of g(X,Y) and g(V,W) : ",cor(wts[1:(n/2)],wts[(n/2+1):n])))
# negative correlation here is good, but not very strong

thetaA2 <- numeric(Nsim)
for(i in 1:Nsim){
  X1 <- runif(n/2,min=-1.0,max=1.0)
  X <- c(X1,shift(X1))
  Y1 <- runif(n/2,min=-1.0,max=1.0)
  Y <- c(Y1,shift(Y1))
  thetaA2[i] <- 4.0*mean(X^2+Y^2<=1.0)
}
print(paste0("SD of shift-based antithetic MC integral, n = ",n," : ",sd(thetaA2)))
print(paste0("Theoretical SD of regular MC integral, n = ",n,", : ",sqrt(pi*(4-pi)/n)))
# i.e. a rather modest improvement, in line with the rather weak negative correlation


#
# 1.h)
#
n.sim <- 10000
n <- 1000
ImpSest <- function(sigma,n,n.sim){
  ests <- numeric(n.sim)
  for(i in 1:n.sim){
    X <- rnorm(n=n,sd=sigma)
    Y <- rnorm(n=n,sd=sigma)
    wts <- (X^2+Y^2<=1.0)*(1.0/(dnorm(X,sd=sigma)*dnorm(Y,sd=sigma)))
    ests[i] = mean(wts)
  }
  return(ests)
}
est1 <- ImpSest(0.3,n,n.sim)
print(paste0("Importance sampling mean, sigma = 0.3 : ",mean(est1)))
print(paste0("Importance sampling SD, sigma = 0.3 : ",sd(est1)))

est2 <- ImpSest(0.624127,n,n.sim)
print(paste0("Importance sampling mean, sigma = 0.624127 : ",mean(est2)))
print(paste0("Importance sampling SD, sigma = 0.624127 : ",sd(est2)))

est3 <- ImpSest(1.0,n,n.sim)
print(paste0("Importance sampling mean, sigma = 1.0 : ",mean(est3)))
print(paste0("Importance sampling SD, sigma = 1.0 : ",sd(est3)))

print(paste0("Crude MC would give SD : ",sqrt(pi*(4-pi)/n)))
# in this case, the family of importance densities considered are not competitive
# with the crude Monte Carlo estimator.
} # end problem 1

#
# problem 2
#

if(prob2){  
  # the intensity function
  lambda <- function(t){
    return((29.7*(1.0+cos(2*pi*(t+0.1)))*(1.0-0.5*exp(-0.1*t)))+0.6)
  }
  # checks in footnote
  #lambda(0.5)
  #lambda(1.0)
  
  
  #
  # point 2.a
  #
  message("2a")
  #the distribution is Poisson(m), where m = \int_0^1 \lambda(t) dt
  m <- integrate(f=lambda,lower=0.0,upper=1.0)$value
  print(paste0("X has a Poisson distribution with parameter : ",m))
  m2 <- integrate(f=lambda,lower=5.0,upper=6.0)$value
  print(paste0("expected # storms in 2025 : ",m2 ))
  m3 <- integrate(f=lambda,lower=0.0,upper=2.0)$value
  print(paste0("expected # storms in 2020 through 2021 : ", m3))
  print(paste0("SD of # storms in 2020 through 2021 : ", sqrt(m3)))
  
  message("2c")
  
  # function simulating a single trajectory
  simStormTimes <- function(a,b){
    if(a<0.0 || b<= a) error("bad arguments")
    lambdaMax <- 60.0 # bound from previous point
    Ns <- rpois(1,lambda=lambdaMax*(b-a)) # number of events for the homogenous process
    s <- sort(runif(n=Ns,min=a,max=b)) # the events in the homogenous process
    u <- runif(n=Ns) # accept probs
    return(s[lambda(s)/lambdaMax > u])
  }
  
  # simulations to check results
  N.sim <- 10000
  Nstorms <- numeric(N.sim)
  
  # case 1: the year 2020
  a <- 0.0
  b <- 1.0
  for(i in 1:N.sim){
    Nstorms[i] <- length(simStormTimes(a,b))
  }
  hist(Nstorms,probability = TRUE,breaks=-0.5 + 0:(max(Nstorms)+1))
  xg <- 0:max(Nstorms)
  points(xg,dpois(xg,lambda = m),col="red") # fairly good accordance
  
  
  # case 2: the year 2025
  a <- 5.0
  b <- 6.0
  for(i in 1:N.sim){
    Nstorms[i] <- length(simStormTimes(a,b))
  }
  print(paste0("mean # storms, 2025 : ",mean(Nstorms)," , true : ",m2))
  
  # case 3: 2020 and 2021 combined
  a <- 0.0
  b <- 2.0
  for(i in 1:N.sim){
    Nstorms[i] <- length(simStormTimes(a,b))
  }
  print(paste0("mean # storms, 2020 and 2021 : ",mean(Nstorms)," , true : ",m3))
  print(paste0("SD # storms, 2020 and 2021 : ",sd(Nstorms)," , true : ",sqrt(m3)))
  
  
  #
  # 2.d
  #
  message("2d")
  
  # mean size of claims
  cfun <- function(t){return(10.0*exp(0.05*t))}
  
  # the simulation code
  claimSize <- function(a,b){
    stimes <- simStormTimes(a,b)
    nstorms <- length(stimes)
    if(nstorms>0){
      csize <- sum(rexp(nstorms, rate=1.0/cfun(stimes)))
    } else {
      csize <- 0.0
    }
    return(csize)
  }
  
  set.seed(1)
  n.sim <- 100000
  claims <- numeric(n.sim)
  a <- 0.0
  b <- 1.0
  for(i in 1:n.sim){
    claims[i] <- claimSize(a,b)
  }
  
  print(paste0("Claim size expectation: ",mean(claims)))
  print(paste0("Claim size standard deviation: ",sd(claims)))
  cil <- mean(claims) - sd(claims)/sqrt(n.sim)
  cir <- mean(claims) + sd(claims)/sqrt(n.sim)
  print(paste0("Claim size expectation CI : [",cil,"  ",cir,"]"))
  # the required fund to cover claims with 97.5% probability is the 0.975 quantile
  print(paste0("Claim size 0.975 quantile: ",quantile(claims,probs=c(0.975))))
  
  #
  # 2.e
  #
  message("2e")
  # Rao Blackwellized estimator of mean
  cX <- numeric(n.sim)
  for(i in 1:n.sim){
    ts <- simStormTimes(a,b)
    N <- length(ts)
    if(N>0){
      cX[i] <- sum(cfun(ts)) # E(X|N,storm times)
    } else {
      cX[i] <- 0.0 # no storms      
    }
  }
  print(paste0("claim size expectation based on RB estimator : ",mean(cX)))
  print(paste0("SD of RB estimator : ",sd(cX)/sqrt(n.sim)))
  print(paste0("SD of plain estimator : ",sd(claims)/sqrt(n.sim)))
  # I.e. a reduction of standard deviation, but not a huge difference
  
  #
  # 2 f)
  #
  message("2f)")
  #
  # estimator of variance based on law of total variance
  # see theory solution for background
  #
  cV <- numeric(n.sim) # conditional variances
  cE <- numeric(n.sim) # conditional expectations
  
  for(i in 1:n.sim){
    ts <- simStormTimes(a,b)
    N <- length(ts)
    if(N>0){
      cE[i] <- sum(cfun(ts)) # E(X|N,storm times)
      cV[i] <- sum(cfun(ts)^2) # Var(X|N, storm times)
    } else {
      cE[i] <- 0.0 # no storms
      cV[i] <- 0.0 # no storms
    }
  }
  V.estimate <- mean(cV) + var(cE) 
  print(paste0("RB estimate of variance : ",V.estimate))
  print(paste0("plain estimate of variance : ",var(claims)))
  # fairly similar; now repeat several times to compare:
  nrep <- 100 # number of replications
  n.sim <- 100 # reduce number of simulations to make it go faster
  cE <- numeric(n.sim)
  cV <- numeric(n.sim)
  claims <- numeric(n.sim)
  V.plain <- numeric(nrep)
  V.RB <- numeric(nrep)
  for(rep in 1:nrep){
    for(i in 1:n.sim){
      ts <- simStormTimes(a,b)
      N <- length(ts)
      if(N>0){
        claims[i] <- sum(rexp(N, rate=1.0/cfun(ts))) # plain
        cE[i] <- sum(cfun(ts)) # E(X|N,storm times)
        cV[i] <- sum(cfun(ts)^2) # Var(X|N, storm times)
      } else {
        claims[i] <- 0.0 # no storms
        cE[i] <- 0.0 # no storms
        cV[i] <- 0.0 # no storms
      }
    }
    V.plain[rep] <- var(claims) # regular MC estimate
    V.RB[rep] <- mean(cV) + var(cE) # RB estimate
  }
  
  
  print(paste0("SD of plain Variance estimator : ",sd(V.plain)))
  print(paste0("SD of RB variance estimator : ",sd(V.RB)))
  # quite large improvement!
  print(paste0("mean of plain Variance estimator : ",mean(V.plain)))
  print(paste0("mean of RB Variance estimator : ",mean(V.RB)))  
  # rather similar, no worries about correctness.
} # end problem 2

#
# problem 3
#
if(prob3){
  message("3a")
  
  # function for obtaining B bootstrap samples of R-squared
  # based on data set data and the provided formula
  boot.Rsq <- function(data,formula,B){
    Rsqs <- numeric(B) # output vector
    n <- dim(data)[1] # number of observations in data set
    for(i in 1:B){
      data.star <- data[sample.int(n,replace=TRUE),]
      lm.obj <- lm(formula = formula,data=data.star)
      Rsqs[i] <- summary(lm.obj)$r.squared
    }
    return(Rsqs)
  }
  
  load("prob23.dat") # creates data-frame df
  # obtain bootstrap samples from function above
  Rsqs <- boot.Rsq(data=df,formula=y~x1+x2+x3+x4+x5,B=5000)
  # make histogram
  hist(Rsqs,probability = TRUE,breaks=30,
       main="histogram of bootstrap samples for Rsquared")
  # histogram has a heavier left hand tail, not well-approximated
  # by a Normal distribution
  
  message("3b")
  # first compute R^2 based on the original data
  lm.out <- lm(y~x1+x2+x3+x4+x5,data=df)
  Rsq.hat <- summary(lm.out)$r.squared
  # bootstrap bias
  bs.bias <- mean(Rsqs) - Rsq.hat   
  print(paste0("Bootstrap-based bias : ",bs.bias))
  
  # bootstrap standard deviation
  print(paste0("Bootstrap-based SD : ",sd(Rsqs)))
  
  # standard normal bootstrap interval for R^2
  zah <- qnorm(p=0.5*0.01,lower.tail = FALSE)
  nci <- c(Rsq.hat - zah*sd(Rsqs),Rsq.hat + zah*sd(Rsqs))
  print("standard normal bootstrap interval for R.squared ")
  print(nci)
  # precentile bootstrap interval
  pci <- quantile(Rsqs,probs = c(0.005,1-0.005))
  print("percentile bootstrap interval for R.squared")
  print(pci)
  
  # comparison
  # lenght of the intervals
  print(paste0("length of standard normal interval : ",nci[2]-nci[1]))
  print(paste0("length of standard normal interval : ",pci[2]-pci[1]))
  # percentile interval slightly shorter
  
  
  print("interval difference (normal-percentile")
  print(nci-pci)
  # normal interval somewhat shifted to the right relative to percentile interval.
  # this difference is due to the non-Gaussian features of the bootstrap
  # distribution, where we see from the histogram a heavier left tail, and
  # thin right tail.
  
} # end problem 3

#
# problem 4
#
if(prob4){
  message("4b")
  #
  # inverse transform algorithm
  #
  rltrunc_normal <- function(n,l,mu,sigma){
    Fl <- pnorm(q=l,mean=mu,sd=sigma)
    return(qnorm(p=(Fl + runif(n)*(1.0-Fl)),mean=mu,sd=sigma))
  }
  #
  # density function
  #
  dltrunc_normal <- function(r,l,mu,sigma){
    Fl <- pnorm(q=l,mean=mu,sd=sigma)
    return(dnorm(x=r,mean=mu,sd=sigma)*(r>=l)/(1.0-Fl))
  }
  
  l <- 1.8
  mu <- 2.2
  sigma <- 1.5
  
  x <- rltrunc_normal(100000,l,mu,sigma) # random draws
  xg <- seq(from=1.5, to=max(x)+1, length.out = 1000)
  hist(x,xlim=c(0,max(x)+1),probability = TRUE,
       breaks = seq(from=l,to=max(x)+1,by=0.25))
  lines(xg,dltrunc_normal(xg,l,mu,sigma))
  
  #
  # functions for calulating mean and variance
  # notice parameters are taken from global scope above
  fm <- function(r){
    return(r*dltrunc_normal(r,l,mu,sigma))
  }
  fm2 <- function(r){
    return(r^2*dltrunc_normal(r,l,mu,sigma))
  }
  Rmean <- integrate(fm,lower=l,upper=15.0)
  Rm2 <- integrate(fm2,lower=l,upper=15.0)
  RSD <- sqrt(Rm2$value - Rmean$value^2)
  
  print(paste0("R mean : ",Rmean$value))
  print(paste0("R mean, simulations : ",mean(x)))
  print(paste0("R SD : ",RSD))
  print(paste0("R SD, simulations : ",sd(x)))
  # i.e. in good accordance
  
  message("4c")
  #
  # simulate random variates from joint distribution
  #
  rRS <- function(n){
    R <- rltrunc_normal(n,1.8,2.2,1.5)
    S <- runif(n,min=0.5/R,0.8/R)
    return(cbind(R,S))
  }
  #
  # joint probability density function (needed later)
  #
  dRS <- function(r,s){
    return(dltrunc_normal(r,1.8,2.2,1.5)*dunif(s,min=0.5/r,max=0.8/r))
  }
  
  simRS <- rRS(1000)
  par(mfrow=c(1,2))
  plot(simRS[,"R"],simRS[,"S"],pch=19,cex=0.05,xlab="R",ylab="S",
       main="scatterplot")
  simRS <- rRS(100000)
  print("mean vector")
  print(colMeans(simRS))
  print("covariance matrix")
  print(cov(simRS))
  print("correlation matrix")
  print(cor(simRS))
  
  #
  # 4.d)
  #
  message("4d)")
  hist(rRS(100000)[,"S"],probability = TRUE,xlab="S",main="S-marginal")
  # from histogram; it is immediatly clear that the marginal
  # distribution of S is not uniform
  
  # function that evaluates the marginal density of S 
  dSmarg <- Vectorize(function(s){return(integrate(f=dRS,lower=1.8,upper=15,s)$value)})
  sgrid <- seq(from=0.0,to=0.5,length.out = 100)
  lines(sgrid,dSmarg(sgrid),col="red")
  
  message("4e")
  source("seir_model.R")
  
  nsim <- 500 # total number of simulations
  l <- 1.8
  mu <- 2.2
  sigma <- 1.5
  
  # construct the antithetics
  Fl <- pnorm(q=l,mean=mu,sd=sigma)
  U.R <- runif(nsim/2)
  R <- c(qnorm(p=(Fl + U.R*(1.0-Fl)),mean=mu,sd=sigma),
        qnorm(p=(Fl + (1-U.R)*(1.0-Fl)),mean=mu,sd=sigma))
  S.left <- 0.5/R
  S.right <- 0.8/R
  
  uu <- runif(nsim/2)
  U.S <- c(uu,1-uu) 
  # notice a + (b-a)*U is uniform(a,b) when U is uniform(0,1)
  S <- S.left + (S.right-S.left)*U.S
  
  # run simulations
  M <- numeric(nsim)
  H <- matrix(0.0,nsim,365) # prepare for next point
  for(i in 1:nsim){
    ff <- seir_sim(R0max = R[i],soc_dist_Rfac = S[i])
    t1 <- which(ff$dates=="2021-01-01")
    t2 <- which(ff$dates=="2021-12-31")
    H[i,] <- rowSums(1000000*ff$ode[t1:t2,c("HH","HC","CC")])
    M[i] <- max(H[i,])
  }
  
  M1 <- M[1:(nsim/2)]
  M2 <- M[(nsim/2+1):nsim]
  print(paste0("correlation of antithetic function values for M: ",cor(M1,M2)))
  print(paste0("correlation of antithetic function values for M^2: ",cor(M1^2,M2^2)))
  # both are slightly negative, indicating that the antithetics are 
  # a somewhat good idea in this situation, but will not make a large difference
  # 
  
  #
  # 4.f)
  #
  message("4f")
  pdf("example_traject.pdf")
  par(mfrow=c(1,1))
  # make a data frame to facilitate plotting with dates
  df <- data.frame(t(H),dates=ff$dates[t1:t2])
  
  # in the upper plot, simply show some of the trajectories:
  plot(df$dates,df$X1,type="l",
       xlab="dates in 2021",
       ylab="required hospital beds",
       main="example trajetories") 
  for(i in 2:50) lines(df$dates,df[,i],col=i)
  dev.off()
  # compute daily statistics for plotting
  daily.stats <- matrix(0.0,365,6)
  for(t in 1:365){
    daily.stats[t,1] <- mean(H[,t])
    daily.stats[t,2] <- median(H[,t])
    daily.stats[t,3:6] <- quantile(H[,t],probs=c(0.005,0.1,0.9,0.995))
  }
  
  df2 <- data.frame(daily.stats,ff$dates[t1:t2])
  colnames(df2) <- c("mean","median","q05","q10","q90","q995","dates")
  pdf("uncertainty.pdf")
  plot(df2$dates,df2$mean,ylim=c(0,300),col=0,
       xlab="dates in 2021",
       ylab="required hospital beds",
       main="distribution representation") # simply sets up plot window
  polygon(c(df2$dates,rev(df2$dates)),c(df2$q05,rev(df2$q995)),
          col="red",density=100)
  polygon(c(df2$dates,rev(df2$dates)),c(df2$q10,rev(df2$q90)),
          col="blue",density=100)
  lines(df2$dates,df2$mean,lwd=2)
  lines(df2$dates,df2$median,col="green")
  legend("bottomright",lty=c(1,1,1,1),lwd=c(10,10,2,1),
         legend=c("99%","90%","mean","median"),
         col=c("red","blue","black","green"),)
  dev.off()
} # end problem 4

