rm(list=ls())

#IMPORTANT MESSAGE:
#Some problems has to be ran/executed individually,
#as some of the functions/variables has the same name

#Problem 1A
func_1a <- function(x){
  density <- ((2^(1/4)*gamma(3/4))/pi)*exp(-(x^4/2))
}

#MC
nsim <- 10000
x <- rnorm(nsim)
func_1a <- func_1a(x)
plot(x, func_1a)

###################################################################

#problem 1B
func_1b <- function(x){
  density <- (2*x*exp(-x^2))
}
nsim <- 10000
x <- rnorm(nsim)
func_1b <- func_1b(x)
hist(func_1b)

mean_1b <- mean(func_1b)
print(mean_1b)

#integrate(func_1b, -Inf, Inf)

###################################################################
#problem 1C

#Problem 2A

###################################################################

#Problem 3A
#load the data set
df <- data.frame(read.table(file="logistic_regression_data.txt"))
x <- df$x
y <- df$y

#Function returning a log-posterior kernel for theta = (alpha, beta)
logistic.lp <- function(theta){
  alpha <- theta[1]
  beta <- theta[2]
  #log-likelihood
  Eeta <- exp(alpha + beta*x)
  #print(Eeta)
  p <- Eeta/(1.0+Eeta)
  log.like <- sum(dbinom(y, size = 1, prob = p, log = TRUE))
  
  #priors
  log.prior <- dnorm(alpha, sd = 10, log = TRUE) + dnorm(beta, sd = 10, log = TRUE)
  
  #log-posterior kernel
  return(log.like + log.prior)
}

theta_hat <- c(-0.102, 1.993)
sigma <- matrix(c(0.00653, -0.00058, -0.00058, 0.01689), ncol = 2)

#>install.packagees("mvtnorm") is necessary
library(mvtnorm)

#independent 2D MH algorithms
two_D_MH <- function(lp, sigma, theta_, n){
  out <- matrix(0, n, 2)
  #print(out)
  out[1, ] <- theta_
  k <- 0
  
  #old/previous
  old <- lp(out[ , 1]) - dmvnorm(theta_, mean = theta_, sigma = sigma, log = TRUE)
  
  for(i in 2:n){
    theta_star <- rmvnorm(1, mean = theta_, sigma = sigma)
    new <- lp(theta_star) - dmvnorm(theta_star, mean = theta_, sigma = sigma, log = TRUE)
    
    alpha <- exp(min(0, new - old))
    if (runif(1) < alpha && is.finite(alpha)){
      out[i, ] <- theta_star
      new <- old
      #print(new)
      k <- k + 1
    } else{
      out[i, ] <- out[i - 1, ]
    }
  }
  print(k/n)
  return(out)
}
output <- two_D_MH(logistic.lp, sigma = (1/2)*sigma, theta_ = theta_hat, n =10000)

# Acceptance rate is 0.99. 
# Acceptance rate for independent MH should be as high as possible

ts.plot(output)

print("Covariance matrix: ")
(1/2)*sigma

print("Estimated covariance matrix: ")
cov(output)

#>install.packages("coda") is necessary

library(coda)
ESS <- function(x){ return(as.numeric(coda::effectiveSize(x))) }
ESS(output)

chain <- mcmc(data = output)
varnames(chain) <- c("alpha", "beta")   

geweke.diag(chain)   #Geweke's Covergence Diagnostic
geweke.plot(chain)   #Geweke-Brooks Plot

###################################################################

#Problem 3B
alpha <- chain[ ,1]
beta <- chain[ ,2]

m <- function(x){
  exp(alpha + beta*x)/(1 + exp(alpha + beta*x))
}

x <- seq(-5, 5, length.out = 1000)
val <- 1000  #value

medi <- numeric(val)   #Median variable
q5 <- numeric(val)     #0.05 quantiles
q5_val <- 0.05
q95 <- numeric(val)    #0.95 quantiles
q95_val <- 0.95

for(i in 1:val){
  medi[i] <- median(m(x[i]))
  q5[i] <- quantile(m(x[i]), q5_val)
  q95[i] <- quantile(m(x[i]), q95_val)
}

plot(x, medi, type = "l", lwd = 2)
lines(x, q5, col = "red", lwd = 2)
lines(x, q95, col = "blue", lwd = 2)
legend(-5, 1, legend = c("95%", "50%", "5%"), col = c("blue", "black", "red"), lty = 1, cex = 0.8)

###################################################################

#Problem 3C
for(i in 1:1000){
  sum_ <- sum((m(x[i]) > 0.8))
  certainty <- sum_/length(alpha)
  if(certainty >= 0.99){
    value <- x[i]
    break
  }
}

print(value)

###################################################################

#Problem 4A

###################################################################
#Problem 4B
df <- data.frame(read.table(file="linear_regression_data.txt"))
x <- df$x
y <- df$y


x_sum <- sum(x)
y_sum <- sum(y)
xy_sum <- sum(x*y)
xx_sum <- sum(x^2)
n <- length(x)

gibbs <- function(theta, n.iter){
  alpha <- theta[1]
  beta <- theta[2]
  tau <- theta[3]
  
  out <- matrix(0, n.iter, 3)
  out[1, ] <- c(alpha, beta, tau)
  
  
  x_sum <- sum(x)
  y_sum <- sum(y)
  xy_sum <- sum(x*y)
  xx_sum <- sum(x^2)
  n <- length(x)
  
  for(i in 2:n.iter){
    alpha <- rnorm(1, mean = (tau*(y_sum - beta*x_sum))/(tau*n + 0.01), sd = sqrt(1/(tau*n + 0.01)))
    beta <- rnorm(1, mean = (tau*(xy_sum - alpha*x_sum))/(tau*xx_sum + 0.01), sd = sqrt(1/(tau*xx_sum + 0.01)))
    var <- sum((y - alpha - beta*x)^2)
    tau <- rgamma(1, shape = (n/2 + 1), scale = (1/(var/2 + 1)))
    
    out[i, ] <- c(alpha, beta, tau)
  }
  return(out)
}

###################################################################

#Problem 4C
result <- gibbs(theta = c(1, -1, 0.5), n.iter = 10000)
print(result)

ts.plot(result)

#>install.packages("coda") is necessary

chain <- mcmc(result)
varnames(chain) <- c("alpha", "beta", "tau")
geweke.diag(chain)
geweke.plot(chain)

library(coda)
ESS <- function(x){ return(as.numeric(coda::effectiveSize(x))) }
ESS(result)

lm.result <- lm(y~x, data = df)
mean(result[ ,1])

lm.result$coeffincients[1]
mean(result[ ,2])

lm.result$coefficients[2]
colMeans(result)

###################################################################

#Problem 4D

###################################################################

#Problem 4E
mean(result[ ,1]) + mean(out[ ,2])

###################################################################

#Problem 5A

#Problem 6A

#Problem 6B

#Problem 6C

#Problem 6D
