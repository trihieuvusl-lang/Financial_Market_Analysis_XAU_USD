rm(list=ls())
require(zoo)
require(ggplot2)
require(rugarch)
## install.packages("ggplot2)
source("LoadFundData.R")

# Panel of figures for full sample 
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(P, main="Level", xlab="", ylab="")
plot(r, main="returns", xlab="", ylab="")


##########################################################
# Settings for analysis                                  #
##########################################################

T  <- length(r)                        # full sample, nobs
N  <- 1250                             # we shorten the sample to the last 5 years
r  <- tail(r,N)                        # log-returns
R  <- coredata(r)

##########################
# Descriptive statistics #
##########################

mu = sum(R)/N

R0    <- R - mu
M2    <- sum(R0^2)/N
M3    <- sum(R0^3)/N
M4    <- sum(R0^4)/N

sig <- sqrt(M2)       ## volatility
S   <- M3/(sig^3)     ## skewness
K   <- M4/(sig^4)     ## kurtosis

# the same with package moments
require(moments)
M1   = moment(R, order = 1, central = FALSE, na.rm = TRUE)
M2   = moment(R, order = 2, central = TRUE, na.rm = TRUE)
M3   = moment(R, order = 3, central = TRUE, na.rm = TRUE)
M4   = moment(R, order = 4, central = TRUE, na.rm = TRUE)
mu0  = M1
sig0 = sqrt(M2)
S0   = M3/(sig0^3)
K0   = M4/(sig0^4)

# Annualization
Nyear <- 250
muA   <- mu*Nyear
sigA  <- sig*sqrt(Nyear)


###################
# Normality test  #
###################

#D'Agostino test for skewness
agostino.test(R)

#Anscombe-Glynn test of kurtosis for normal samples
anscombe.test(R)

#Jarque-Bera test of normality
jarque.test(R)

#############################################
# t-Student distribution                    #
#############################################

x = seq(-5,5,0.01)
junk <- cbind(x,dnorm(x),ddist("std",x,shape=10),ddist("std",x,shape=5),ddist("std",x,shape=3))
colnames(junk) <- c("x","vInf","v10","v5","v3")

ggplot(data = data.frame(junk), aes(x = x)) + 
  geom_line(aes(y = vInf,colour="v=Inf")) + 
  geom_line(aes(y = v10,colour="v=10")) +
  geom_line(aes(y = v5,colour="v=5")) +
  geom_line(aes(y = v3,colour="v=3")) +
  labs(title="t-Student distribution", y="", x="", caption="")+
  theme_bw()+
  scale_colour_manual("", 
                      breaks = c("v=Inf", "v=10","v=5","v=3"),
                      values = c("black", "red", "blue", "green"))
## recmommend this function 
# ?? dist / ?? rnorm to check mấy cái hàm ở console 

######################
# Empirical density  #
######################

R0        <- (R-mu)/sig #standarize 
bwdth     <- 0.1

ggplot(data.frame(R0), aes(x = R0)) +
  theme_bw() +
  geom_histogram(binwidth = bwdth, colour = "white", fill = "yellow4", size = 0.1) +
  stat_function(fun = function(x) dnorm(x)*N*bwdth, color = "red", size = 1) +                       # normal distribution
  stat_function(fun = function(x) ddist("std",x,shape=5)*N*bwdth, color = "black", size = 1) +     # t-Student distribution
  xlim(-7,-2)


#######################################################
# quantile-quantile plot (QQplot)                     #
#######################################################

R0     <- (R-mu)/sig

# quantile table
q              <- seq(0.01, 0.20, 0.01)
QQtable        <- data.frame(Qemp = quantile(R0,q),Qteo = qnorm(q))   
names(QQtable) <- c("empirical quantile", "theoretical quantile")

# QQ plot
q            <- seq(0.001, 0.999, 0.001)
Qemp          <- quantile(R0,q)                    # empirical quantiles 
#Qemp         <- quantile(rnorm(N),q)              # hypotherical if we had data from normal distribution
#Qemp         <- quantile(R0,q)                    
#Qteo          <- qnorm(q)                          # theoretical quantile (normal)
v=5;   Qteo  <- qt(q,v)*sqrt((v-2)/v)             # t-Student (variance is v/v-2)
lim0    <- c(-5,5)                           
par(mfrow=c(1,1), cex = 0.7, bty="l")
plot(Qemp,Qteo, main="QQplot", col="red", xlim = lim0, ylim = lim0,
     xlab="Qemp", ylab="Qteo") 
abline(a=0,b=1, lwd=2)


##############################################
# Estimation of t-Studenta parameters        #
##############################################

# Method of moments (K - kurtosis)
v0 <- 4 + 6/(K-3)

# Max likelihood method 
require(MASS)

d0 <- fitdistr(R0, "normal")
d0$loglik

d1 <- fitdistr(R0, "t", m = 0, start = list(s=sqrt((v0-2)/v0), df=v0), lower=c(0.001,3))
d1$loglik
v=d1$estimate[[2]]

# Comparison
(d1$loglik-d0$loglik)/N

# QQ plot
q        <- seq(0.001, 0.999, 0.001)
Qemp     <- quantile(R0,q)                     
Qteo     <- qt(q,v)*sqrt((v-2)/v)              
lim0    <- c(-5,5)                             
par(mfrow=c(1,1), cex = 0.7, bty="l")
plot(Qemp,Qteo, main="QQplot", col="red", xlim = lim0, ylim = lim0,
     xlab="Qemp", ylab="Qteo") 
abline(a=0,b=1, lwd=2)

########################################
# Autocorrelations                     #
########################################

require(tseries)

# ACF for returns
z0 <- acf(R,20, plot=TRUE)
plot(z0$acf[2:20], type="h", main="ACF for returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

# ACF for squarred returns 
z1 <- acf(R^2,20, plot=FALSE)
plot(z1$acf[2:20], type="h", main="ACF for sq returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.4), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

# Cross-correlations
z2 <- ccf(R^2, R, 20, type="correlation", plot=TRUE)
plot(z2$acf[22:40], type="h", main="correlation for sq. ret and ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

plot(z2$acf[20:1], type="h", main="correlation for ret and sq. ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

# panel of plots
par(mfrow=c(2,2), cex = 0.7, bty="l")
plot(z0$acf[2:20], type="h", main="ACF for returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z1$acf[2:20], type="h", main="ACF for sq returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.4), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z2$acf[22:40], type="h", main="correlation for sq. ret and ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z2$acf[20:1], type="h", main="correlation for ret and sq. ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))


# Ljunga-Box test

Box.test(R, lag = 4, type = c("Ljung-Box"))
Box.test(R^2, lag = 4, type = c("Ljung-Box"))




