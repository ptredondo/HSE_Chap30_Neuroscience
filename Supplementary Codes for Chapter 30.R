##############################################
##### Supplementary Codes for Chapter 30 #####
##############################################

### Setting working directory
setwd("set/working/directory/")

##### Loading needed packages
library(qgam)
library(evgam)
library(VGAM)
library(signal)
library(LSTS)
library(optimParallel)
library(ggplot2)

### Creating clusters for parallel computations
cl <- makeCluster(detectCores()-1); setDefaultCluster(cl = cl)

##################################
##### User-defined functions #####
##################################

### Extracting the five standard frequency bands.
# INPUTS:
# x - numeric vector containing series to be filtered
# fs - sampling rate
# rhythms - character vector of frequency bands
# method - "butter" for Butterworth filter, "fir" for finite impulse response filter
# twosided - logical value (TRUE/FALSE)
# OUTPUT: list of filtered frequency band signals

extract_rhythms <- function(x,fs,rhythms = c("delta","theta","alpha","beta","gamma"),
                            method = "butter",twosided = TRUE){
  
  y <- vector("list",length(rhythms))
  
  for(i in 1:length(rhythms)){
    
    if(rhythms[i] == "delta"){
      
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(0.5,4)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(0.5,4)/(0.5 * fs), type="pass")
      }
      
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
      
    } else if(rhythms[i] == "theta"){
      
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(4,8)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(4,8)/(0.5 * fs), type="pass")
      }
      
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
      
    } else if(rhythms[i] == "alpha"){
      
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(8,12)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(8,12)/(0.5 * fs), type="pass")
      }
      
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
      
    } else if(rhythms[i] == "beta"){
      
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(12,30)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(12,30)/(0.5 * fs), type="pass")
      }
      
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
      
    } else if(rhythms[i] == "gamma"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(30,45)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(30,45)/(0.5 * fs), type="pass")
      }
      
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
      
    } 
    
  }
  
  return(y)
  
}

### Calculating time-varying relative power of the five standard frequency bands.
# INPUTS:
# y - numeric vector containing series of interest
# blength - block length of the moving segments for periodogram calculations
# samprate - sampling rate
# smooth - logical value (TRUE/FALSE)
# OUTPUT: matrix containing the time-varying relative spectral power (row) for 
#         each frequency band (column)
tv.relpo <- function(y,blength = 500,samprate = 100,smooth = TRUE){
  
  bnum <- length(y)/blength
  relpo <- array(NA,dim=c(bnum,5))
  
  for(b in 1:bnum){
    sub <- y[blength*(b-1) + 1:blength]
    sm.perio <- smooth.periodogram(sub,plot = TRUE,spar = 0.5)$smooth.periodogram
    omega <- (1:(blength/2))/blength
    fband.limit <- c(0,4,8,12,30,samprate/2)/samprate
    for(f in 1:5){
      f_ind <- (1:(blength/2))[omega > fband.limit[f] & omega <= fband.limit[f+1]]
      relpo[b,f] <- sum(sm.perio[f_ind])/sum(sm.perio)
    }
  }
  
  if(smooth){
    for(f in 1:5){
      relpo[,f] <- smooth.spline(1:bnum, relpo[,f], spar=0.5)$y
    }
  }
  
  return(relpo)
}

### Negative log-likelihood of H&T model without covariates
# INPUTS:
# par - parameters of the H&T model
# Y0 - numeric vector of the conditioning variable (Y1 given Y0)
# Y1 - numeric vector of the conditioned variable (Y1 given Y0)
# OUTPUT: value of the negative log-likelihood
nll0 <- function(par,Y0,Y1){
  
  alpha <- par[1]
  beta <- par[2]
  sig <- par[3]
  mu <- par[4]
  
  #Add in any constraints here
  if(alpha < -1 | alpha > 1 | beta > 1 | beta < 0 | sig < 0){return(1e10)}
  
  nll <- sum(0.5*log(2*pi) + log(sig*Y0^beta) + 0.5*((Y1 - (alpha*Y0 + mu*Y0^beta))/(sig*Y0^beta))^2)
  
  if(is.finite(nll)){
    return(nll)
  } else{return(1e10)}
}

### Negative log-likelihood of H&T model with covariates
# INPUTS:
# par - parameters of the H&T model
# Y0 - numeric vector of the conditioning variable (Y1 given Y0)
# Y1 - numeric vector of the conditioned variable (Y1 given Y0)
# X - matrix of covariates containing splines variables
# X0 - matrix of covariates containing non-splines variables
# OUTPUT: value of the negative log-likelihood
nllx <- function(par,Y0,Y1,X,X0,lambda = 0){
  
  n.cov <- ncol(X)
  n.cov2 <- ncol(X0)
  alpha.vec <- par[1:n.cov]
  beta.vec <- par[n.cov + 1:n.cov2]
  sig.vec <- par[n.cov + n.cov2 + 1:n.cov2]
  mu.vec <- par[n.cov + 2*n.cov2 + 1:n.cov2]
  
  alpha_t <- X %*% matrix(alpha.vec,ncol = 1)
  beta_t <- X0 %*% matrix(beta.vec,ncol = 1)
  sig_t <- X0 %*% matrix(sig.vec,ncol = 1)
  mu_t <- X0 %*% matrix(mu.vec,ncol = 1)
  
  #Add in any constraints here
  if(any(alpha_t < -1) | any(alpha_t > 1 ) | any(beta_t > 1) | any(beta_t < 0) | any(sig_t < 0)){return(1e10)}
  
  nll <- sum(0.5*log(2*pi) + log(sig_t*Y0^beta_t) + 
               0.5*((Y1 - (alpha_t*Y0 + mu_t*Y0^beta_t))/(sig_t*Y0^beta_t))^2) +
    lambda*sum((alpha.vec[1:(n.cov-2)] - 2*alpha.vec[2:(n.cov-1)] + alpha.vec[3:n.cov])^2)
  
  if(is.finite(nll)){
    return(nll)
  } else{return(1e10)}
}

### Fitting Heffernen and Tawn model with b-splines as covariates
# INPUTS:
# Y0 - numeric vector of the conditioning variable (Y1 given Y0)
# Y1 - numeric vector of the conditioned variable (Y1 given Y0)
# X - matrix of covariates containing splines variables
# X0 - matrix of covariates containing non-splines variables
# init.par - initial parameters for estimation
# lambda_grid - numeric vector of lambda values for AIC selection
# u_thres - numeric value in (0,1) which serves as the quantile threshold
# OUTPUT: list containing "param.est" - parameter estimates,
#         "param.t" - parameter estimates per time t, "AIC" - AIC of fitter model
htgam <- function(Y0,Y1,X,X0,init.par = NULL,
                  lambda_grid = seq(0,2,by = 0.05),u_thres = 0.9){
  
  Tlength <- length(Y0)
  n.cov <- ncol(X)
  n.cov2 <- ncol(X0)
  mydata <- data.frame(time = 1:Tlength,Y0 = Y0,Y1 = Y1,X,X0)
  y_thres <- qlaplace(u_thres)
  mydata_ht <- subset(mydata,Y0 > y_thres)
  Y0.excess <- matrix(mydata_ht$Y0,ncol = 1)
  Y1.excess <- matrix(mydata_ht$Y1,ncol = 1)
  X.excess <- as.matrix(mydata_ht[,c(1:n.cov + 3)])
  X0.excess <- as.matrix(mydata_ht[,c(n.cov + 3 + 1:n.cov2)])
  
  if(is.null(init.par)){
    init.par0 <- runif(4)
    fit0 <- optimParallel(par=init.par0,fn=nll0,Y0=Y0.excess,Y1=Y1.excess,lower = c(-1,0))
    init.par <- c(fit0$par[1],rep(0.5,n.cov-1),fit0$par[2],0,fit0$par[3],0,fit0$par[4],0)
  }
  
  AIC_all <- array(NA,dim = length(lambda_grid))
  par_all <- array(NA,dim = c(n.cov+3*n.cov2,length(lambda_grid)))
  
  for(j in 1:length(lambda_grid)){
    fit1<-optimParallel(par=init.par,fn=nllx,Y0=Y0.excess,Y1=Y1.excess,X=X.excess,X0=X0.excess,
                        lambda=lambda_grid[j],lower = c(-1,0))
    par_all[,j] <- fit1$par
    AIC_all[j] <- 2*(length(fit1$par)+2) + 2*fit1$value
  }
  
  final.ind <- which.min(AIC_all)
  AIC <- AIC_all[final.ind]
  par.est <- par_all[,final.ind]
  
  alpha.est <- par.est[1:n.cov]
  beta.est <- par.est[n.cov + 1:n.cov2]
  sig.est <- par.est[n.cov + n.cov2 + 1:n.cov2]
  mu.est <- par.est[n.cov + 2*n.cov2 + 1:n.cov2]
  param.est <- list(alpha.est = alpha.est, beta.est = beta.est,
                    sig.est = sig.est, mu.est = mu.est)
  
  alpha_t <- X %*% matrix(alpha.est,ncol = 1)
  beta_t <- X0 %*% matrix(beta.est,ncol = 1)
  sig_t <- X0 %*% matrix(sig.est,ncol = 1)
  mu_t <- X0 %*% matrix(mu.est,ncol = 1)
  param.t <- list(alpha_t = alpha_t, beta_t = beta_t, sig_t = sig_t, mu_t = mu_t)
  
  output <- list(param.est = param.est, param.t = param.t, AIC = AIC)    
  
  return(output)
  
}

##############################################################################
##### The following codes are used to implement the procedures discussed #####
##### in Chapter 30 (in the same order as they appear in the chapter).   #####
##############################################################################

##### Importing dataset and pre-loaded objects
load("Epilepsy Data.RData")

#########################
##### Section 1.3.1 #####
#########################

### Fitting 0.9-quantile regression for X_t (amplitudes from channel T3) in terms
### of a smooth function of time t.
mydata <- data.frame(X_t = abs(ep[,"T3"]),time = 1:Tlength,
                     seizure = c(rep(0,35000),rep(1,15000)))
qgam.fit <- qgam(X_t~s(time, k=10, bs="bs"),data = mydata,qu = 0.90)
mydata$thres <- predict(qgam.fit, newdata = mydata)
mydata$excess <- mydata$X_t - mydata$thres
mydata$excess <- ifelse(mydata$excess > 0,mydata$excess,0)
mydata$ecdf <- rank(mydata$X_t)/(length(mydata$X_t)+1)

### Figure 1.4
par(mfrow = c(1,2),mar = c(5,5,2,2))
plot(mydata$time/samprate,mydata$X_t,type = "l",
     xlab = "Time (in seconds)", ylab = "Amplitude",
     main = "Time-varying Threshold (90-th Quantile)",
     lwd = 2)
lines(mydata$time/samprate,mydata$thres,col = "blue",lwd = 3)
abline(v = 350,col = "orange",lwd = 2,lty = 2)
plot(mydata$time/samprate,mydata$excess,type = "l",
     xlab = "Time (in seconds)", ylab = "Excess",
     main = "Exceedances above the Threshold",
     lwd = 2)
abline(v = 350,col = "orange",lwd = 2,lty = 2)

### Fitting GPD to the excess of X_t over the time-varying 0.9-quantile threshold.
mydata_gpd <- subset(mydata,excess > 0)
fit_gpd <- evgam(list(excess ~ s(time, k=10,bs = "bs"), ~ seizure),
                 data=mydata_gpd, family="gpd")

### Q-Q plot with standard Laplace margins
param.est <- predict(fit_gpd,mydata_gpd,type = "response")
scale.est <- param.est[,1]; shape.est <- param.est[,2]
u_t <- VGAM::pgpd(mydata_gpd$excess,scale = scale.est,shape = shape.est)
x_emp <- sort(u_t); y_emp <- qlaplace(x_emp)
x_the <- (1:length(u_t))/(length(u_t)+1); y_the <- qlaplace(x_the)

### Generating the 0.995-quantile over time based on the 0.9-quantile regression
### and fitted GPD for the excesses over the time-varying threshold.
param.est <- predict(fit_gpd,mydata,type = "response",se.fit = TRUE)
scale.est <- param.est$fitted[,1]
shape.est <- param.est$fitted[,2]
mydata$q_0.995 <- mydata$thres + VGAM::qgpd(0.95,scale = scale.est,shape = shape.est)

### Confidence bounds for the parameter estimates
scale.up <- param.est$fitted[,1] + 2*param.est$se.fit[,1]
scale.lo <- param.est$fitted[,1] - 2*param.est$se.fit[,1]
shape.up <- param.est$fitted[,2] + 2*param.est$se.fit[,2]
shape.lo <- param.est$fitted[,2] - 2*param.est$se.fit[,2]

### Figure 1.5
par(mfrow = c(1,2),mar = c(5,5,2,2))
plot_min <- min(y_the,y_emp,na.rm = TRUE) # plot limits
plot_max <- max(y_the,y_emp,na.rm = TRUE) # plot limits
plot(y_the,y_emp,ylim = c(plot_min,plot_max),xlim = c(plot_min,plot_max),
     xlab = "Theoretical", ylab = "Empirical",
     main = "Q-Q Plot for the Fitted GPD Model")
abline(a = 0,b = 1,lwd = 2.5)
mean(mydata$X_t > mydata$q_0.995) # Prob(X_t > q_0.995) = 0.00422

plot(mydata$time/samprate,mydata$X_t,type = "l",lwd = 2,
     xlab = "Time (in seconds)", ylab = "Amplitude",
     main = "Estimated 99.5-th Quantile over Time")
abline(v = 350,lwd = 2,lty = 2,col = "orange")
lines(mydata$time/samprate,mydata$thres,col = "blue",lty = 2,lwd = 3)
lines(mydata$time/samprate,mydata$q_0.995,col = "purple",lty = 2,lwd = 3)
text(150,10,expression(P(Y(t) > q[0.995](t)) == 0.0042))

### Figure 1.6
par(mfrow = c(1,2),mar = c(5,5,2,2))
plot(mydata$time/samprate,scale.est,type = "l",lwd = 3,ylim = c(0,1),
     xlab = "Time (in Seconds)", ylab = "",
     main = "Time-varying Scale Parameter")
polygon(x = c(mydata$time/samprate, rev(mydata$time/samprate)),
        y = c(scale.up, rev(scale.lo)),
        col = adjustcolor("gray",alpha.f = 0.5), border = NA)
lines(mydata$time/samprate,scale.est,type = "l",lwd = 3)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
abline(h = 0,lty = 2,lwd = 1.25)

plot(mydata$time/samprate,shape.est,type = "l",lwd = 3,ylim = c(-0.1,0.35),
     xlab = "Time (in Seconds)", ylab = "",
     main = "Time-varying Shape Parameter")
polygon(x = c(mydata$time/samprate, rev(mydata$time/samprate)),
        y = c(shape.up, rev(shape.lo)),
        col = adjustcolor("gray",alpha.f = 0.5), border = NA)
lines(mydata$time/samprate,shape.est,type = "l",lwd = 3)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
abline(h = 0,lty = 2,lwd = 1.25)

### Extracting the five standard frequency bands of channel T3 over locally
### stationary blocks of length 500
blength <- 500
bnum <- Tlength/blength
X <- ep$T3
X.f <- array(NA,dim = c(Tlength,5))
for(b in 1:bnum){
  sub <- X[blength*(b-1) + 1:blength]
  f.sub <- extract_rhythms(sub,fs = samprate)
  for(f in 1:5){
    X.f[blength*(b-1) + 1:blength,f] <- f.sub[[f]]
  }
}

### Estimating the time-varying scale parameter for each frequency band.
### Results of the following codes are pre-loaded with the .RData file
### Uncomment if needed
# scale.f <- array(NA,dim = c(Tlength,5))
# for(f in 1:5){
#   mydata <- data.frame(time = 1:Tlength,X_t = abs(X.f[,f]),
#                        seizure = c(rep(0,35000),rep(1,15000)))
#   qgam.fit <- qgam(X_t~s(time, k=10, bs="bs"),data = mydata,qu = 0.9)
#   mydata$thres <- predict(qgam.fit, newdata = mydata)
#   mydata$excess <- mydata$X_t - mydata$thres
#   mydata$excess <- ifelse(mydata$excess > 0,mydata$excess,0)
# 
#   mydata_gpd <- subset(mydata,excess > 0)
#   fit_gpd <- evgam(list(excess ~ s(time, k=10,bs = "bs"), ~ seizure)
#                    , data=mydata_gpd, family="gpd")
# 
#   param.est <- predict(fit_gpd,mydata,type = "response")
#   scale.f[,f] <- param.est[,1]
# }

### Calculating relative spectral power for channel T3 
relpo <- tv.relpo(y = ep$T3,blength = 500,samprate = samprate)

### Figure 1.7
par(mfrow = c(1,2),mar = c(5,5,2,2))
plot((1:Tlength)/samprate,scale.f[,1],type = "l",ylim = c(0,1),lwd = 2,
     main = "Scale Parameter per Frequency Band",
     ylab = "Estimated Scale",xlab = "Time (in Seconds)")
lines((1:Tlength)/samprate,scale.f[,2],col = "red",lwd = 2)
lines((1:Tlength)/samprate,scale.f[,3],col = "blue",lwd = 2)
lines((1:Tlength)/samprate,scale.f[,4],col = "green",lwd = 2)
lines((1:Tlength)/samprate,scale.f[,5],col = "purple",lwd = 2)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
legend("top", legend=c("Delta", "Theta", "Alpha", "Beta", "Gamma"),
       col=c("black","red", "blue","green","purple"),
       lty = 1,cex = 0.5,horiz = TRUE,lwd = 2)

plot(1:bnum*blength/samprate,relpo[,1],type = "l",ylim = c(0,1.05),lwd = 2,
     main = "Relative Spectral Band Power (Smooth)",
     ylab = "Relative Spectral Power",xlab = "Time (in seconds)")
lines(1:bnum*blength/samprate,relpo[,2],col = "red",lwd = 2)
lines(1:bnum*blength/samprate,relpo[,3],col = "blue",lwd = 2)
lines(1:bnum*blength/samprate,relpo[,4],col = "green",lwd = 2)
lines(1:bnum*blength/samprate,relpo[,5],col = "purple",lwd = 2)
abline(v = 35,lty = 2,col = "orange",lwd = 2)
legend("top", legend=c("Delta", "Theta", "Alpha", "Beta", "Gamma"),
       col=c("black","red", "blue","green","purple"),
       lty = 1,cex = 0.5,horiz = TRUE,lwd = 2)

#########################
##### Section 1.3.2 #####
#########################

### Transforming signal amplidues of each channel to the standard Laplace scale
blength <- 1250
bnum <- Tlength/blength
channels <- c("F7","F8","T3")
mydata_lap <- array(NA,dim = c(Tlength,3))
for(ch in 1:3){
  mydata <- data.frame(time = 1:Tlength,X_t = abs(ep[,channels[ch]]),
                       seizure = c(rep(0,35000),rep(1,15000)))
  qgam.fit <- qgam(X_t~s(time, k=10, bs="bs"),data = mydata,qu = 0.9)
  mydata$thres <- predict(qgam.fit, newdata = mydata)
  mydata$excess <- mydata$X_t - mydata$thres
  mydata$excess <- ifelse(mydata$excess > 0,mydata$excess,0)
  
  mydata$ecdf <- rep(NA,Tlength)
  for(b in 1:bnum){
    mydata$ecdf[blength*(b-1) + 1:blength] <- rank(mydata$X_t[blength*(b-1) + 1:blength])/(blength+1)
  }
  
  mydata_gpd <- subset(mydata,excess > 0)
  fit_gpd <- evgam(list(excess ~ s(time, k=10,bs = "bs"), ~ seizure),
                   data=mydata_gpd, family="gpd")
  param.est <- predict(fit_gpd,mydata,type = "response")
  scale.est <- param.est[,1]
  shape.est <- param.est[,2]
  mydata$gpdcdf <- VGAM::pgpd(mydata$excess,scale = scale.est,shape = shape.est)
  u_t <- ifelse(mydata$excess > 0,
                0.9 + (1 - 0.9)*mydata$gpdcdf,
                0.9*mydata$ecdf)
  mydata_lap[,ch] <- qlaplace(u_t)
}

### Estimating cross-extremograms
y_thres <- qlaplace(0.95)
mydata_pre <- mydata_lap[1:35000,] > y_thres
mydata_post <- mydata_lap[35001:50000,] > y_thres
l_max <- 50
time.pre <- (l_max+1):(35000-l_max)
time.post <- (l_max+1):(15000-l_max)
cxgram <- array(NA,dim = c(2,2,2*l_max+1))
for(ch in 1:2){
  for(l in (-l_max):l_max){
    cxgram[ch,1,l_max + 1 + l] <- sum(mydata_pre[time.pre,ch] & mydata_pre[time.pre-l,3])/sum(mydata_pre[time.pre-l,3])
    cxgram[ch,2,l_max + 1 + l] <- sum(mydata_post[time.post,ch] & mydata_post[time.post-l,3])/sum(mydata_post[time.post-l,3])
  }
}

### Figure 1.8a
par(mfrow = c(1,2),mar = c(5,5,2,2),oma = c(0,0,2,0))
plot((-l_max):l_max,cxgram[1,1,],type = "h",lwd=2.5,ylim = c(0,1),
     main = "Before Seizure Onset",ylab = expression(chi[u](h)),xlab = "Lag h")
abline(v = 0,col = "gray",lty = 2)
plot((-l_max):l_max,cxgram[1,2,],type = "h",lwd=2.5,ylim = c(0,1),
     main = "After Seizure Onset",ylab = expression(chi[u](h)),xlab = "Lag h")
abline(v = 0,col = "gray",lty = 2)
mtext("(a) Extremogram between F7 and T3",side = 3,font = 2,
      outer = TRUE,adj = 0.1)

### Figure 1.8b
par(mfrow = c(1,2),mar = c(5,5,2,2),oma = c(0,0,2,0))
plot((-l_max):l_max,cxgram[2,1,],type = "h",lwd=2.5,ylim = c(0,1),
     main = "Before Seizure Onset",ylab = expression(chi[u](h)),xlab = "Lag h")
abline(v = 0,col = "gray",lty = 2)
plot((-l_max):l_max,cxgram[2,2,],type = "h",lwd=2.5,ylim = c(0,1),
     main = "After Seizure Onset",ylab = expression(chi[u](h)),xlab = "Lag h")
abline(v = 0,col = "gray",lty = 2)
mtext("(b) Extremogram between F8 and T3",side = 3,font = 2,
      outer = TRUE,adj = 0.1)

### Extracting the five standard frequency bands of channels F7, F8 and T3 over 
### locally stationary blocks
blength <- 500
bnum <- Tlength/blength
X <- ep[,c("F7","F8","T3")]
X.f <- array(NA,dim = c(Tlength,5,3))
for(ch in 1:3){
  for(b in 1:bnum){
    sub <- X[blength*(b-1) + 1:blength,ch]
    f.sub <- extract_rhythms(sub,fs = samprate)
    for(f in 1:5){
      X.f[blength*(b-1) + 1:blength,f,ch] <- f.sub[[f]]
    }
  }
}

### Estimating band coherence for each frequency band
l_max <- 50
coh_comp <- array(NA,dim = c(2,5,bnum,2*l_max+1))
for(f in 1:5){
  for(b in 1:bnum){
    X.sub <- X.f[blength*(b-1) + 1:blength,f,]
    time.sub <- (l_max+1):(blength-l_max)
    for(l in (-l_max):l_max){
      coh_comp[1,f,b,l + l_max + 1] <- (cor(X.sub[time.sub,1],X.sub[time.sub - l,3]))^2
      coh_comp[2,f,b,l + l_max + 1] <- (cor(X.sub[time.sub,2],X.sub[time.sub - l,3]))^2
    }
    
  }
}
coh_est <- apply(coh_comp,MARGIN = c(1,2,3),FUN = max,na.rm = TRUE)

### Smoothing band coherence estimates for plotting
sm.coh_est <- array(NA,dim=c(2,5,bnum))
for(f in 1:5){
  sm.coh_est[1,f,] <- smooth.spline(1:bnum, coh_est[1,f,], spar=0.5)$y
  sm.coh_est[2,f,] <- smooth.spline(1:bnum, coh_est[2,f,], spar=0.5)$y
}

### Figure 1.9a
par(mfrow = c(1,2),mar = c(5,5,2,2))
plot(1:bnum*blength/samprate,sm.coh_est[1,1,],type = "l",ylim = c(0,1.05),lwd = 2,
     main = "Band Coherence Between F7 and T3 (Smooth)",
     ylab = "Band Coherence",xlab = "Time (in seconds)")
lines(1:bnum*blength/samprate,sm.coh_est[1,2,],col = "red",lwd = 2)
lines(1:bnum*blength/samprate,sm.coh_est[1,3,],col = "blue",lwd = 2)
lines(1:bnum*blength/samprate,sm.coh_est[1,4,],col = "green",lwd = 2)
lines(1:bnum*blength/samprate,sm.coh_est[1,5,],col = "purple",lwd = 2)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
legend("top", legend=c("Delta", "Theta", "Alpha", "Beta", "Gamma"),
       col=c("black","red", "blue","green","purple"),
       lty = 1,cex = 0.5,horiz = TRUE,lwd = 2)

plot(1:bnum*blength/samprate,sm.coh_est[2,1,],type = "l",ylim = c(0,1.05),lwd = 2,
     main = "Band Coherence Between F8 and T3 (Smooth)",
     ylab = "Band Coherence",xlab = "Time (in seconds)")
lines(1:bnum*blength/samprate,sm.coh_est[2,2,],col = "red",lwd = 2)
lines(1:bnum*blength/samprate,sm.coh_est[2,3,],col = "blue",lwd = 2)
lines(1:bnum*blength/samprate,sm.coh_est[2,4,],col = "green",lwd = 2)
lines(1:bnum*blength/samprate,sm.coh_est[2,5,],col = "purple",lwd = 2)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
legend("top", legend=c("Delta", "Theta", "Alpha", "Beta", "Gamma"),
       col=c("black","red", "blue","green","purple"),
       lty = 1,cex = 0.5,horiz = TRUE,lwd = 2)

### Transforming the frequeny bands of each channel to the standard Laplace scale
### Results of the following codes are pre-loaded with the .RData file
### Uncomment if needed
# blength <- 1250
# bnum <- Tlength/blength
# channels <- c("F7","F8","T3")
# mydata.f_lap <- array(NA,dim = c(Tlength,5,3))
# for(ch in 1:3){
#   for(f in 1:5){
#     mydata <- data.frame(time = 1:Tlength,X_t = abs(X.f[,f,ch]),
#                          seizure = c(rep(0,35000),rep(1,15000)))
#     qgam.fit <- qgam(X_t~s(time, k=10, bs="bs"),data = mydata,qu = 0.9)
#     mydata$thres <- predict(qgam.fit, newdata = mydata)
#     mydata$excess <- mydata$X_t - mydata$thres
#     mydata$excess <- ifelse(mydata$excess > 0,mydata$excess,0)
#     
#     mydata$ecdf <- rep(NA,Tlength)
#     for(b in 1:bnum){
#       mydata$ecdf[blength*(b-1) + 1:blength] <- rank(mydata$X_t[blength*(b-1) + 1:blength])/(blength+1)
#     }
#     
#     mydata_gpd <- subset(mydata,excess > 0)
#     fit_gpd <- evgam(list(excess ~ s(time, k=10,bs = "bs"), ~ seizure),
#                      data=mydata_gpd, family="gpd")
#     param.est <- predict(fit_gpd,mydata,type = "response")
#     scale.est <- param.est[,1]
#     shape.est <- param.est[,2]
#     mydata$gpdcdf <- VGAM::pgpd(mydata$excess,scale = scale.est,shape = shape.est)
#     u_t <- ifelse(mydata$excess > 0,
#                   0.9 + (1 - 0.9)*mydata$gpdcdf,
#                   0.9*mydata$ecdf)
#     mydata.f_lap[,f,ch] <- qlaplace(u_t)
#     cat(paste("DONE: Freq ",f," Ch ",ch,sep=""),"\r")
#     
#   }
# }

### Estimating band tail coherence for each frequency band
l_max <- 50
blength <- 1250
bnum <- Tlength/blength
y_thres <- qlaplace(0.95)
tailcoh_comp <- array(NA,dim = c(2,5,bnum,2*l_max+1))
for(f in 1:5){
  for(b in 1:bnum){
    X.sub <- mydata.f_lap[blength*(b-1) + 1:blength,f,]
    
    time.sub <- (l_max+1):(blength-l_max)
    for(l in (-l_max):l_max){
      if(sum(X.sub[time.sub-l,3] > y_thres) > 0){
        tailcoh_comp[1,f,b,l + l_max + 1] <- sum(X.sub[time.sub,1] > y_thres & X.sub[time.sub-l,3] > y_thres)/sum(X.sub[time.sub-l,3] > y_thres)
        tailcoh_comp[2,f,b,l + l_max + 1] <- sum(X.sub[time.sub,2] > y_thres & X.sub[time.sub-l,3] > y_thres)/sum(X.sub[time.sub-l,3] > y_thres)
      } else {
        tailcoh_comp[1,f,b,l + l_max + 1] <- 0
        tailcoh_comp[2,f,b,l + l_max + 1] <- 0
      }
    }
    
  }
}
tailcoh_est <- apply(tailcoh_comp,MARGIN = c(1,2,3),FUN = max,na.rm = TRUE)

### Smoothing tail coherence estimates for plotting
sm.tailcoh_est <- array(NA,dim=c(2,5,bnum))
for(f in 1:5){
  sm.tailcoh_est[1,f,] <- smooth.spline(1:bnum, tailcoh_est[1,f,], spar=0.6)$y
  sm.tailcoh_est[2,f,] <- smooth.spline(1:bnum, tailcoh_est[2,f,], spar=0.6)$y
}

### Figure 1.9b
par(mfrow = c(1,2),mar = c(5,5,2,2))
plot(1:bnum*blength/samprate,sm.tailcoh_est[1,1,],type = "l",ylim = c(0,1.05),lwd = 2,
     main = "Band Tail Coherence Between F7 and T3 (Smooth)",
     ylab = "Band Tail Coherence",xlab = "Time (in seconds)")
lines(1:bnum*blength/samprate,sm.tailcoh_est[1,2,],col = "red",lwd = 2)
lines(1:bnum*blength/samprate,sm.tailcoh_est[1,3,],col = "blue",lwd = 2)
lines(1:bnum*blength/samprate,sm.tailcoh_est[1,4,],col = "green",lwd = 2)
lines(1:bnum*blength/samprate,sm.tailcoh_est[1,5,],col = "purple",lwd = 2)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
legend("top", legend=c("Delta", "Theta", "Alpha", "Beta", "Gamma"),
       col=c("black","red", "blue","green","purple"),
       lty = 1,cex = 0.5,horiz = TRUE,lwd = 2)

plot(1:bnum*blength/samprate,sm.tailcoh_est[2,1,],type = "l",ylim = c(0,1.05),lwd = 2,
     main = "Band Tail Coherence Between F8 and T3 (Smooth)",
     ylab = "Band Tail Coherence",xlab = "Time (in seconds)")
lines(1:bnum*blength/samprate,sm.tailcoh_est[2,2,],col = "red",lwd = 2)
lines(1:bnum*blength/samprate,sm.tailcoh_est[2,3,],col = "blue",lwd = 2)
lines(1:bnum*blength/samprate,sm.tailcoh_est[2,4,],col = "green",lwd = 2)
lines(1:bnum*blength/samprate,sm.tailcoh_est[2,5,],col = "purple",lwd = 2)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
legend("top", legend=c("Delta", "Theta", "Alpha", "Beta", "Gamma"),
       col=c("black","red", "blue","green","purple"),
       lty = 1,cex = 0.5,horiz = TRUE,lwd = 2)

#########################
##### Section 1.3.3 #####
#########################

### Preparing bootstrap data based on block bootstrap with 
### varying quantile thresholds for fitting the GPD model
### Results of the following codes are pre-loaded with the .RData file
### Uncomment if needed
# set.seed(3) # seed for reproducibility
# 
# blength <- 1250 # block length for ECDF estimation
# bnum <- 50000/blength
# slength <- 25 # segment length for block bootstrap
# snum <- blength/slength
# 
# mydata_lap_bs <- array(NA,dim = c(50000,3,500))
# channels <- c("F7","F8","T3")
# for(bs in 1:500){
#   
#   # Block bootstrap
#   mydata_bs <- array(NA, dim = c(50000,3))
#   for(b in 1:bnum){
#     sub <- as.matrix(ep[blength*(b-1) + 1:blength,channels])
#     s.ind <- sample(1:snum,snum,replace = TRUE)
#     mydata_bs[blength*(b-1) + 1:blength,] <- sub[rep((s.ind-1)*slength,each = slength) + rep(1:slength,times = snum),]
#   }
#   
#   # Random quantile threshold (Uniform(0.9,0.95))
#   u0_thres <- runif(1,min = 0.9,max = 0.95)
#   
#   # Transforming bootstrap data to the standard Laplace scale
#   for(ch in 1:3){
#     mydata <- data.frame(time = 1:Tlength,X_t = abs(mydata_bs[,ch]),
#                          seizure = c(rep(0,35000),rep(1,15000)))
#     qgam.fit <- qgam(X_t~s(time, k=10, bs="bs"),data = mydata,qu = u0_thres)
#     mydata$thres <- predict(qgam.fit, newdata = mydata)
#     mydata$excess <- mydata$X_t - mydata$thres
#     mydata$excess <- ifelse(mydata$excess > 0,mydata$excess,0)
#     
#     mydata$ecdf <- rep(NA,Tlength)
#     for(b in 1:bnum){
#       mydata$ecdf[blength*(b-1) + 1:blength] <- rank(mydata$X_t[blength*(b-1) + 1:blength])/(blength+1)
#     }
#     
#     mydata_gpd <- subset(mydata,excess > 0)
#     fit_gpd <- evgam(list(excess ~ s(time, k=10,bs = "bs"), ~ seizure),
#                      data=mydata_gpd, family="gpd")
#     param.est <- predict(fit_gpd,mydata,type = "response")
#     scale.est <- param.est[,1]
#     shape.est <- param.est[,2]
#     mydata$gpdcdf <- VGAM::pgpd(mydata$excess,scale = scale.est,shape = shape.est)
#     u_t <- ifelse(mydata$excess > 0,
#                   u0_thres + (1 - u0_thres)*mydata$gpdcdf,
#                   u0_thres*mydata$ecdf)
#     mydata_lap_bs[,ch,bs] <- qlaplace(u_t)
#   }
# }

### Fitting H&T model for F7 given T3
set.seed(1)
X <- bs(1:Tlength,df = 10,intercept = TRUE) # b-splines
X0 <- cbind(1,c(rep(0,35000),rep(1,15000))) # intercept and seizure indicator
Y0 = mydata_lap[,3] # channel T3
Y1 = mydata_lap[,1] # channel F7
ht.fit <- htgam(Y0 = Y0,Y1 = Y1, X = X,X0 = X0,lambda_grid = seq(0,2,by=0.05))

### Bootstrap estimates for H&T model parameters
### Results of the following codes are pre-loaded with the .RData file
### Uncomment if needed
# f7_alpha_t_bs <- array(NA,dim = c(50000,500))
# f7_beta_bs <- array(NA,dim = c(2,500))
# f7_sig_bs <- array(NA,dim = c(2,500))
# f7_mu_bs <- array(NA,dim = c(2,500))
# 
# for(bs in 1:500){
#   
#   Y0 = mydata_lap_bs[,3,bs] # channel T3
#   Y1 = mydata_lap_bs[,1,bs] # channel F7
#   ht.fit_bs <- htgam(Y0 = Y0,Y1 = Y1, X = X,X0 = X0,
#                   lambda_grid = seq(0,2,by=0.25),
#                   u_thres = runif(1,min = 0.88,max = 0.92))
#   
#   f7_alpha_t_bs[,bs] <- ht.fit_bs$param.t$alpha_t
#   f7_beta_bs[,bs] <- c(ht.fit_bs$param.est$beta.est[1],
#                     ht.fit_bs$param.est$beta.est[1]+ht.fit_bs$param.est$beta.est[2])
#   f7_sig_bs[,bs] <- c(ht.fit_bs$param.est$sig.est[1],
#                    ht.fit_bs$param.est$sig.est[1]+ht.fit_bs$param.est$sig.est[2])
#   f7_mu_bs[,bs] <- c(ht.fit_bs$param.est$mu.est[1],
#                   ht.fit_bs$param.est$mu.est[1]+ht.fit_bs$param.est$mu.est[2])
# }
# 
# ### Confidence bounds for the time-varying alpha parameter
# f7_alpha_t_up <- apply(f7_alpha_t_bs,MARGIN = 1,FUN = quantile,probs = 0.975,na.rm = TRUE)
# f7_alpha_t_lo <- apply(f7_alpha_t_bs,MARGIN = 1,FUN = quantile,probs = 0.025,na.rm = TRUE)

### Figure 1.10
par(mfrow = c(2,2),mar = c(5,5,2,2))
plot(NA,NA,type = "l",lwd = 3,xlim = c(0,Tlength/samprate),ylim = c(-1,1),
     xlab = "Time (in Seconds)", ylab = "",
     main = "Time-varying Alpha Parameter")
polygon(x = c(1:Tlength/samprate, rev(1:Tlength/samprate)),
        y = c(f7_alpha_t_up, rev(f7_alpha_t_lo)),
        col = adjustcolor("gray",alpha.f = 0.5), border = NA)
lines(1:Tlength/samprate,ht.fit$param.t$alpha_t,type = "l",lwd = 3)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
abline(h = 0,lty = 2,lwd = 1.25)

hist(f7_beta_bs[1,],col = adjustcolor("blue",alpha = 0.5),
     xlim = c(0,0.8),ylim = c(0,6),freq = FALSE,
     main = "Bootstrap Beta Parameter",xlab = "", ylab = "")
mtext(expression(hat(beta)),side = 1,line = 4,cex = 2)
hist(f7_beta_bs[2,],col = adjustcolor("red",alpha = 0.5),add = TRUE,freq = FALSE)
legend("topright", legend=c("Before Onset", "After Onset"),
       fill=c(adjustcolor("blue",alpha = 0.5),adjustcolor("red",alpha = 0.5)),
       cex = 0.75,horiz = FALSE)

hist(f7_mu_bs[1,],col = adjustcolor("blue",alpha = 0.5),
     xlim = c(-0.6,0.4),ylim = c(0,6),freq = FALSE,
     main = "Bootstrap Mu Parameter",xlab = "", ylab = "")
mtext(expression(hat(mu)),side = 1,line = 4,cex = 2)
hist(f7_mu_bs[2,],col = adjustcolor("red",alpha = 0.5),add = TRUE,freq = FALSE)
legend("topright", legend=c("Before Onset", "After Onset"),
       fill=c(adjustcolor("blue",alpha = 0.5),adjustcolor("red",alpha = 0.5)),
       cex = 0.75,horiz = FALSE)

hist(f7_sig_bs[1,],col = adjustcolor("blue",alpha = 0.5),
     xlim = c(0.3,1.5),ylim = c(0,6),freq = FALSE,
     main = "Bootstrap Sig Parameter",xlab = "", ylab = "")
mtext(expression(hat(sigma)),side = 1,line = 4,cex = 2)
hist(f7_sig_bs[2,],col = adjustcolor("red",alpha = 0.5),add = TRUE,freq = FALSE)
legend("topright", legend=c("Before Onset", "After Onset"),
       fill=c(adjustcolor("blue",alpha = 0.5),adjustcolor("red",alpha = 0.5)),
       cex = 0.75,horiz = FALSE)

### Fitting H&T model for F8 given T3
set.seed(2)
X <- bs(1:Tlength,df = 10,intercept = TRUE) # b-splines
X0 <- cbind(1,c(rep(0,35000),rep(1,15000))) # intercept and seizure indicator
Y0 = mydata_lap[,3] # channel T3
Y1 = mydata_lap[,2] # channel F8
ht.fit <- htgam(Y0 = Y0,Y1 = Y1, X = X,X0 = X0,lambda_grid = seq(0,2,by=0.05))

### Bootstrap estimates for H&T model parameters
### Results of the following codes are pre-loaded with the .RData file
### Uncomment if needed
# f8_alpha_t_bs <- array(NA,dim = c(50000,500))
# f8_beta_bs <- array(NA,dim = c(2,500))
# f8_sig_bs <- array(NA,dim = c(2,500))
# f8_mu_bs <- array(NA,dim = c(2,500))
# 
# for(bs in 1:500){
# 
#   Y0 = mydata_lap_bs[,3,bs] # channel T3
#   Y1 = mydata_lap_bs[,2,bs] # channel F8
#   ht.fit_bs <- htgam(Y0 = Y0,Y1 = Y1, X = X,X0 = X0,
#                      lambda_grid = seq(0,2,by=0.25),
#                      u_thres = runif(1,min = 0.88,max = 0.92))
# 
#   f8_alpha_t_bs[,bs] <- ht.fit_bs$param.t$alpha_t
#   f8_beta_bs[,bs] <- c(ht.fit_bs$param.est$beta.est[1],
#                     ht.fit_bs$param.est$beta.est[1]+ht.fit_bs$param.est$beta.est[2])
#   f8_sig_bs[,bs] <- c(ht.fit_bs$param.est$sig.est[1],
#                    ht.fit_bs$param.est$sig.est[1]+ht.fit_bs$param.est$sig.est[2])
#   f8_mu_bs[,bs] <- c(ht.fit_bs$param.est$mu.est[1],
#                   ht.fit_bs$param.est$mu.est[1]+ht.fit_bs$param.est$mu.est[2])
# }
# 
# ### Confidence bounds for the time-varying alpha parameter
# f8_alpha_t_up <- apply(f8_alpha_t_bs,MARGIN = 1,FUN = quantile,probs = 0.975,na.rm = TRUE)
# f8_alpha_t_lo <- apply(f8_alpha_t_bs,MARGIN = 1,FUN = quantile,probs = 0.025,na.rm = TRUE)

### Figure 1.11
par(mfrow = c(2,2),mar = c(5,5,2,2))
plot(NA,NA,type = "l",lwd = 3,xlim = c(0,Tlength/samprate),ylim = c(-1,1),
     xlab = "Time (in Seconds)", ylab = "",
     main = "Time-varying Alpha Parameter")
polygon(x = c(1:Tlength/samprate, rev(1:Tlength/samprate)),
        y = c(f8_alpha_t_up, rev(f8_alpha_t_lo)),
        col = adjustcolor("gray",alpha.f = 0.5), border = NA)
lines(1:Tlength/samprate,ht.fit$param.t$alpha_t,type = "l",lwd = 3)
abline(v = 350,lty = 2,col = "orange",lwd = 2)
abline(h = 0,lty = 2,lwd = 1.25)

hist(f8_beta_bs[1,],col = adjustcolor("blue",alpha = 0.5),
     xlim = c(0,0.8),ylim = c(0,6),freq = FALSE,
     main = "Bootstrap Beta Parameter",xlab = "", ylab = "")
mtext(expression(hat(beta)),side = 1,line = 4,cex = 2)
hist(f8_beta_bs[2,],col = adjustcolor("red",alpha = 0.5),add = TRUE,freq = FALSE)
legend("topright", legend=c("Before Onset", "After Onset"),
       fill=c(adjustcolor("blue",alpha = 0.5),adjustcolor("red",alpha = 0.5)),
       cex = 0.75,horiz = FALSE)

hist(f8_mu_bs[1,],col = adjustcolor("blue",alpha = 0.5),
     xlim = c(-0.6,0.4),ylim = c(0,6),freq = FALSE,
     main = "Bootstrap Mu Parameter",xlab = "", ylab = "")
mtext(expression(hat(mu)),side = 1,line = 4,cex = 2)
hist(f8_mu_bs[2,],col = adjustcolor("red",alpha = 0.5),add = TRUE,freq = FALSE)
legend("topright", legend=c("Before Onset", "After Onset"),
       fill=c(adjustcolor("blue",alpha = 0.5),adjustcolor("red",alpha = 0.5)),
       cex = 0.75,horiz = FALSE)

hist(f8_sig_bs[1,],col = adjustcolor("blue",alpha = 0.5),
     xlim = c(0.3,1.5),ylim = c(0,6),freq = FALSE,
     main = "Bootstrap Sig Parameter",xlab = "", ylab = "")
mtext(expression(hat(sigma)),side = 1,line = 4,cex = 2)
hist(f8_sig_bs[2,],col = adjustcolor("red",alpha = 0.5),add = TRUE,freq = FALSE)
legend("topright", legend=c("Before Onset", "After Onset"),
       fill=c(adjustcolor("blue",alpha = 0.5),adjustcolor("red",alpha = 0.5)),
       cex = 0.75,horiz = FALSE)

### Fitting H&T model for each frequency band
set.seed(4) # seed for reproducibility
alpha_t <- array(NA,dim = c(50000,2,5))
for(f in 1:5){
  X <- bs(1:Tlength,df = 10,intercept = TRUE) # b-splines
  X0 <- cbind(1,c(rep(0,35000),rep(1,15000))) # intercept and seizure indicator
  Y0 = mydata.f_lap[,f,3] # channel T3
  
  # Fitting H&T model for F7 given T3
  Y1 = mydata.f_lap[,f,1] # channel F7
  ht.fit <- htgam(Y0 = Y0,Y1 = Y1, X = X,X0 = X0,lambda_grid = seq(0,2,by=0.25))
  alpha_t[,1,f] <- ht.fit$param.t$alpha_t
  
  # Fitting H&T model for F8 given T3
  Y1 = mydata.f_lap[,f,2] # channel F8
  ht.fit <- htgam(Y0 = Y0,Y1 = Y1, X = X,X0 = X0,lambda_grid = seq(0,2,by=0.25))
  alpha_t[,2,f] <- ht.fit$param.t$alpha_t
}

### Preparing data for plotting
blength <- 100
bnum <- Tlength/blength
alpha_sub <- alpha_t[(1:bnum)*blength,,] # sub-sampling at every second
data_melt <- NULL
for(f in 5:1){
  for(b in 1:bnum){
    entry <- data.frame(freq = f,time = b,F7 = alpha_sub[b,1,f],F8 = alpha_sub[b,2,f])
    data_melt <- rbind(data_melt,entry)
  }
}

### Color scale for plot
col.pal1 <- colorRampPalette(c("red",
                             "yellow"))
col.pal2 <- colorRampPalette(c("yellow",
                             "white"))
col.pal3 <- colorRampPalette(c("white",
                             "limegreen"))
col.pal4 <- colorRampPalette(c("limegreen",
                             "blue"))

### Plot settings
theme_setting <- theme(axis.text=element_text(size=14),
                       axis.title=element_text(size=16),
                       axis.text.x = element_text(margin = margin(b = 10)),
                       axis.text.y = element_text(margin = margin(l = 20),size=16),
                       plot.title = element_text(face = "bold",size = 20,hjust = -0.125,
                                                 margin = margin(t=20)),
                       legend.title = element_text(face = "bold",size = 20),
                       legend.text=element_text(size=12),
                       legend.key.size = unit(0.75, 'cm'),
                       plot.margin=grid::unit(c(1,1,1,1), "mm"))

### Figure 1.12
data_melt <- subset(data_melt,time > 200)
ggp1 <- ggplot(data_melt, aes(time, freq)) +  xlab("Time (in seconds)") + 
  ylab("Frequency Band") + 
  geom_tile(aes(fill = F7)) + theme_classic(base_line_size = 0) +
  scale_fill_gradientn(colours = c(col.pal1(30),col.pal2(70),col.pal3(70),col.pal4(30)),
                       limits=c(-1,1),name = expression(~hat(alpha))) +
  theme_setting + labs(title="(a) F7 given T3") +
  geom_segment(x = 350,y = 0.5,xend = 350,yend = 5.5, linetype="dashed", color = "orange",lwd=1) + 
  scale_y_continuous(breaks = c(1,2,3,4,5),labels=c(expression(Omega[1]),
                                                    expression(Omega[2]),
                                                    expression(Omega[3]),
                                                    expression(Omega[4]),
                                                    expression(Omega[5])))
ggp1 # time-varying alpha parameter estimates (F7 given T3)

ggp2 <- ggplot(data_melt, aes(time, freq)) +  xlab("Time (in seconds)") + ylab("Frequency Band") + 
  geom_tile(aes(fill = F8)) + theme_classic(base_line_size = 0) +
  scale_fill_gradientn(colours = c(col.pal1(30),col.pal2(70),col.pal3(70),col.pal4(30)),
                       limits=c(-1,1),name = expression(~hat(alpha))) +
  theme_setting + labs(title="(b) F8 given T3") +
  geom_segment(x = 350,y = 0.5,xend = 350,yend = 5.5, linetype="dashed", color = "orange",lwd=1) + 
  scale_y_continuous(breaks = c(1,2,3,4,5),labels=c(expression(Omega[1]),
                                                    expression(Omega[2]),
                                                    expression(Omega[3]),
                                                    expression(Omega[4]),
                                                    expression(Omega[5])))
ggp2 # time-varying alpha parameter estimates (F8 given T3)



