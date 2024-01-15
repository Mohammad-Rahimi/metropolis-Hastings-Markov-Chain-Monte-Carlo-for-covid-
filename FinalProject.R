
rm(list=ls())
library(Bolstad2)

setwd('G:/My Drive/wvu/2nd semester/601/codes/homework11/project')

'''

polandcases<-read.csv("./first_peak_poland.csv")
newdaily<-polandcases$new_cases_smoothed_per_million
#date<-as.Date(polandcases$date,format="%b %d %Y")
date<-seq(1,length(newdaily))
plot(date,newdaily,type="l", main="poland first peak")
#axis(1,date,format(date, "%b %d"), cex.axis=0.8)
####
polandcases2<-read.csv("./second_peak_poland.csv")
newdaily2<-polandcases2$new_cases_smoothed_per_million
#date<-as.Date(polandcases$date,format="%b %d %Y")
date2<-seq(1,length(newdaily2))
plot(date2,newdaily2,type="l", main="poland second peak")
#axis(1,date,format(date, "%b %d"), cex.axis=0.8)

'''
###
polandpeak1<-read.csv("./first_peak_poland_clean.csv")
polandpeak2<-read.csv("./second_peak_poland_clean.csv")
Peak1 <- as.matrix(polandpeak1)
Peak2 <- as.matrix(polandpeak2)
count <- 0

#Define Metropolis-Hastings algorithm

MHmcmc <- function(sigma, posterior, data, steps = 1000, target = 0.2,startValue = NULL, randomSeed = NULL ) #changed here
{
  if (steps < 100) {
    warning("Function should take at least 100 steps")
  }
  #determine number of parameter dimensions
  np <- length(sigma)
  if (any(sigma <= 0)) 
    stop("All standard deviations must be strictly non-zero and positive")
  
  targetSample <- matrix(rep(0, (np+1)*steps), nrow = steps, byrow = TRUE)
  
  if (!is.null(randomSeed)) 
    set.seed(randomSeed)
  z <- rnorm(steps, 0, sigma[1])
  for (n in 2:np){
    z <- cbind(z, rnorm(steps, 0, sigma[n]))
  }
  u <- runif(steps)
  if (is.null(startValue)) 
    startValue <- z[1,]
  
  i1 <- 1
  nstep = 1
  accept = 1
  af <- accept/nstep
  
  targetSample[1,] <- c(startValue, af)
  g <- rep(0, steps)
  proposal <- matrix(rep(0, np*steps), nrow = steps, byrow = TRUE)
  alpha <- rep(0, steps)
  
  g[1] <- posterior(targetSample[1,c(1:np)], data)
  print(g[1])
  for (n in 2:steps) {
    proposal[n,] <- targetSample[i1,c(1:np)] + z[n,]
    g[n] <- posterior(proposal[n,], data)
    k3 <- g[n]
    k4 <- g[i1]
    alpha[n] <- ifelse(k3/k4 > 1, 1, k3/k4)
    
    print( is.na(u[n]) )
    print( is.na(alpha[n]) )
    
    if (u[n] >= alpha[n]) {
      targetSample[n,] <- targetSample[i1,]
    }
    else {
      targetSample[n,] <- c(proposal[n,], af)
      i1 <- n
      accept <- accept + 1
    }
    if (nstep >= 200){
      af <- accept/nstep
      if (af > target){
        z <- z * 1.1
      } else if (af < target){
        z <- z * 0.9
      }
      nstep = 0
      accept = 0
    } else {
      nstep = nstep + 1
    }
  }
  
  oldPar <- par(mfrow = c(2, 2), pty = "s")
  h <- density(targetSample[,1])
  ymax <- max(c(h$y)) * 1.05
  #plot(h$x, h$y, type = "l", col = "light blue", xlim = range(targetSample[,1]), 
     #  ylim = c(0, ymax), main = "Posterior Parameter 1", 
     #  xlab = "x", ylab = "Density")
  
  #h <- hist(targetSample[,2], plot = FALSE)
  #ymax <- max(c(h$density)) * 1.05
  #hist(targetSample[,2], prob = TRUE, col = "light blue", xlim = range(targetSample[,2]), 
     #  ylim = c(0, ymax), main = "Posterior Parameter 2", 
      # xlab = "x", ylab = "Density")
  
  #box()
  #plot(targetSample[,1], type = "l", main = "", ylab = "Target 1 Sample")
  #plot(targetSample[,2], type = "l", main = "", ylab = "Target 2 Sample")
  
  par(oldPar)
  invisible(targetSample)
}

LLH <- function(theta, data){
  print("In LLH")
  
  Yhat <- SIRsimRHS(c(theta[1], theta[2],1000000, 23.577, 0), nrow(data))
  # Yhat <- SIRsimRHS(c(theta[1], theta[2], 37797000, 23.577, 0), nrow(data))
  
  #Yhat <- res
  #Yhat <- theta[1] * data[,2] + theta[2]
  SSE <- sum((log10(data[,2]) - log10(Yhat))^2)
  
  #post <- SSE^(-1/2) + .Machine$double.eps
  post<-exp(-SSE)+ .Machine$double.eps
  
  
  return(post)
  
}

SIRsimRHS<- function(param, T) {
  # Simulate an SIR epidemic
  # beta is infection rate constant, eta is recovery rate constant
  # s0 initial susceptibles, i0 initial infected, r0 initial recovered, 
  # simulation length T
  # returns a matrix size (T+1)*3 with columns s, i, r respectively
  # implicitly deltaT is 1 day
  
  beta <- 10^param[1]
  eta <- 10^param[2]
  #beta <-param[1]
  #eta <-param[2]
  s0 <- param[3]
  i0 <- param[4]
  r0 <- param[5]
  s <- rep(0, T)
  i <- rep(0, T)
  r <- rep(0, T)
  RHS <- rep(0,T)
  
  s[1] <- s0
  i[1] <- i0
  r[1] <- r0
  for (j in 1:(T-1)) { #121
    s[j+1] <- s[j] - beta * s[j] * i[j]
    i[j+1] <- i[j] + beta * s[j] * i[j] - eta * i[j]
    r[j+1] <- r[j] + eta * i[j]
    ############Question???? 
    RHS[j+1] <- beta * s[j] * i[j]  #2-122  # THIS SHOULD NEVER BE NEGATIVE .Machine$double.eps. 1000-->3
    if(RHS[j+1]<0){
       RHS[j+1]<-.Machine$double.eps
     }
    
  }
  
  return(RHS)
}

############################################################################################

#Define data points
#NC <- matrix(c(2,  0.988850026,  2,	1.021569338, 2,	0.989580636, 4,	0.853398211, 4,	0.853895625, 4,	1.28523392, 4,	1.007472244, 
               #10,	1.262924495, 10,	0.997538414, 10,	1.105481178, 10,	0.731030243, 10,	0.90302567), nrow = 12, byrow = TRUE)

#P75 <- matrix(c(2,  2.256871955, 2,	2.45309514, 2,	2.560123986, 4,	2.949653866, 4,	3.384441969, 4,	2.475052574, 4,	2.892547513, 4,	1.904991892, 
                #10,	3.785547475, 10,	5.079586603, 10,	4.495100205, 10,	4.084163797, 10,	3.413358253, 10,	3.922278872), nrow = 14, byrow = TRUE)

#plot(NC[,1], NC[,2], ylim = c(0.5, 6), ty = "p", col = "red")
#points(P75[,1], P75[,2], col = "blue")


#Likelihood evaluation


#Set up standard deviation of proposal distribution
#step length value
Sig <- c(0.001, 0.001)

TNC1 <- MHmcmc(Sig, LLH, Peak1, steps = 4000, target = 0.2, startValue = c(-8.56,0.33))# starting values
#c(-6,-10))
j <- seq(1,nrow(TNC1)-1,by = 10)
plot(j, TNC1[j,1], ty = "l", col = "red")
plot(j, TNC1[j,2], ty = "l", col = "blue")



res <- SIRsimRHS(c(-8.56,0.33, 1e6, 54, 0), 122) 
plot(seq(1, length(res)), res, ty = "l", col = "orange")


#lines(seq(1, nrow(res)), res[,2], col = "blue")
#lines(seq(1, nrow(res)), res[,3], col = "black")

TNC2 <- MHmcmc(Sig, LLH, Peak1, steps = 4000, target = 0.2, startValue = c(runif(1),runif(1)))
TNC3 <- MHmcmc(Sig, LLH, Peak1, steps = 4000, target = 0.2, startValue = c(runif(1),runif(1)))
TNC4 <- MHmcmc(Sig, LLH, Peak1, steps = 4000, target = 0.2, startValue = c(runif(1),runif(1)))



plot(TNC1[,1], TNC1[,2], xlim = c(-0.5, 1), ylim = c(-0.5,1.5), ty = "l", col = "red")
lines(TNC2[,1], TNC2[,2], col = "blue")
lines(TNC3[,1], TNC3[,2], col = "green")
lines(TNC4[,1], TNC4[,2], col = "orange")

#TP751 <- MHmcmc(Sig, LLH, P75, steps = 4000, target = 0.2, startValue = c(0,1))
#TP752 <- MHmcmc(Sig, LLH, P75, steps = 4000, target = 0.2, startValue = c(runif(1),runif(1)))
#TP753 <- MHmcmc(Sig, LLH, P75, steps = 4000, target = 0.2, startValue = c(runif(1),runif(1)))
#TP754 <- MHmcmc(Sig, LLH, P75, steps = 4000, target = 0.2, startValue = c(runif(1),runif(1)))

# Use Gelman-Rubin potential improvement statistic
# Ratio of variance between chains / variance within chain
#

slope <- cbind(TNC1[,1], TNC2[,1], TNC3[,1], TNC4[,1])
Xval <- seq(100, nrow(TNC1), by = 100)
GRx <- rep(0, length(Xval))
for (i in 1:length(Xval)){
  tmp <- GelmanRubin(slope[1:Xval[i],])
  GRx[i] <- tmp$R
}

plot(Xval, GRx, ty = "l", xlim = c(0,4000))

par(mfrow = c(1, 2), pty = "s")
h1 <- density(c(TNC1[2000:4000,1], TNC2[2000:4000,1], TNC3[2000:4000,1], TNC4[2000:4000,1]))
ymax <- max(c(h1$y)) * 1.05
plot(h1$x, h1$y, type = "l", col = "black", xlim = range(-0.5,0.5), 
     ylim = c(0, ymax), main = "Posterior Slope", 
     xlab = "x", ylab = "Density")

h1 <- density(c(TNC1[2000:4000,2], TNC2[2000:4000,2], TNC3[2000:4000,2], TNC4[2000:4000,2]))
ymax <- max(c(h1$y)) * 1.05
plot(h1$x, h1$y, type = "l", col = "black", xlim = range(0.5,1.5), 
     ylim = c(0, ymax), main = "Posterior Intercept", 
     xlab = "x", ylab = "Density")


#################################
'''
rm(list = ls())

SIRsim2 <- function(param, T) {
  # Simulate an SIR epidemic
  # beta is infection rate constant, eta is recovery rate constant
  # s0 initial susceptibles, i0 initial infected, r0 initial recovered, 
  # simulation length T
  # returns a matrix size (T+1)*3 with columns s, i, r respectively
  # implicitly deltaT is 1 day
  beta <- param[1]
  eta <- param[2]
  s0 <- param[3]
  i0 <- param[4]
  r0 <- param[5]
  s <- rep(0, T+1)
  i <- rep(0, T+1)
  r <- rep(0, T+1)
  s[1] <- s0
  i[1] <- i0
  r[1] <- r0
  for (j in 1:T) {
    s[j+1] <- s[j] - beta * s[j] * i[j]
    i[j+1] <- i[j] + beta * s[j] * i[j] - eta * i[j]
    r[j+1] <- r[j] + eta * i[j]
  }
  return(matrix(c(s, i, r), ncol = 3))
}


res <- SIRsimRHS(c(0.5, 0.1, 1, 1e-5, 0), 122)

plot(seq(1, nrow(res)), res[,1], ylim = c(0,1), ty = "l", col = "red")
lines(seq(1, nrow(res)), res[,2], col = "blue")
lines(seq(1, nrow(res)), res[,3], col = "black")
'''