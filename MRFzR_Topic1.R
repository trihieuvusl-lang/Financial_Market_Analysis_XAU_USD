rm(list = ls())

log(16,2)
log(x=16,2)
log(16, base=2)
log(base=2,x=16)


help(log, package=base)
help(log)
?log
args(log)

# Vectors

x <- c(2,1,8,16,4)

length(x)

2*x-2
x^2;    
log(x,2)

x*x
t(x)%*%x

x[c(2,3,4)]
x[2:4]
x[-c(1,5)]

x > 3
all(x>3)
any(x>3)
which(x>3)

names(x) <- c("a","b","c","d","e")
x[c("c","d","e")]

rep(x=1,times=5)
seq(from=2, to=10, by=2)
1:5
rnorm(5)


# Matrices
A <- matrix((1:9)^2, nrow = 3, ncol=3, byrow=TRUE)
B <- matrix(0, nrow = 3, ncol=3)
C <- diag(1,3)

A
t(A)

dim(A)
nrow(A)
ncol(A)

A1 <- A[c(2,3),c(2,3)]
A2 <- A[-2,-2]
A3 <- A[1:2,1:2]

cbind(A1,A2)
rbind(A1,A2)

rownames(A) <- c("r1", "r2", "r3")
colnames(A) <- c("c1", "c2", "c3")
A
A[,"c2"]
A[,"c2",drop=FALSE]

det(A)
eigen(A)
solve(A)
A %*% solve(A)
A * solve(A)

# Lists

x <- seq(from=0, to = 10, by=0.5)
y <- 2 + 3*x + rnorm(21)
z <- lm(y ~ x)
names(z)

MyList <- list(data = rnorm(5), rozklad="normal", par = list(mean=0, sd=1))
MyList
names(MyList)
summary(MyList)

MyList[[1]]
MyList[["data"]]
MyList$data

MyList[[3]][[2]]
MyList[["par"]][["sd"]]
MyList$par$sd

# class
######################

x <- c(2,1,8,16,4)
class(x)
is.numeric(x)
is.character(x)


x.char <- as.character(x)
x.ts   <- as.ts(x)
x.ts   <- ts(x, start=2000, freq=1)
x.date <- as.Date(x, origin="2000-01-01")

summary(x)
summary(x.char)

# conditioning 
##########################

c   <- exp(1)
sig <- 1
if(sig==1) { 
  u <- log(c) 
} else { 
  u <- (c^(1-sig)-1)/(1-sig) 
} 
u

# functions
######################

util <- function(c,sig=1) {
  if(sig==1) { 
    u <- log(c) 
  } else { 
    u <- (c^(1-sig)-1)/(1-sig) 
  } 
  return(u)
}

util(exp(1))
util(exp(1),2)

# loops
############

x <- NULL
for(i in 1:5){
  x <- c(x,i)
}
x

x   <- 2
eps <- 1e-8
while(abs(x^2-x-1) > eps) {
  x <- x - (x^2-x-1)/(2*x-1)
}
x

###########################################
# Time series, packages zoo, Quandl, etc. #
###########################################
rm(list = ls())

require(zoo)

# Set working directory
# setwd(".")


# Quantmod package
require(quantmod)
getSymbols("SPY", src="yahoo")
head(SPY)
class(SPY)

# zoo objects
data <- as.zoo(SPY$SPY.Close)
head(data)
head(index(data))
head(coredata(data))
plot(data)

# Dates
mydate <- as.Date("01-01-2017")
mydate

mydate <- as.Date("2017-01-01")
mydate

mydate <- as.Date("01-01-2017", format="%d-%m-%Y")
mydate

mydate + 3

mydate1 <- as.Date("2017-02-15")
mydate1

mydate < mydate1
mydate >= mydate1
mydate - mydate1

class(mydate - mydate1)

mydate1 - mydate

dates <- index(data)

mydates <- head(dates, n=100)
weekdays(mydates)
months(mydates)
quarters(mydates)

## Sequence of dates
mydates1 <- seq(from=mydate, to=mydate1, by="1 day")
mydates1

# Business days
require(bizdays)
load_quantlib_calendars("Poland", as.Date(head(mydates1,1)), as.Date(tail(mydates1, 1)))
is.bizday(mydates1, "QuantLib/Poland")

## lubridate package
require(lubridate)

mydate <- ymd("2017-02-01")
mydate

mydate1 <- mdy("01/03/2017")
mydate1

yday(mydate)
wday(mydate)

years(5)
mydate + years(5)

mydate + (0:5)*months(2)

## merging zoo objects
data1 <- data
getSymbols("GC=F", src="yahoo")
data2 <- as.zoo(`GC=F`$`GC=F.Close`)


## Sum
data3 <- merge(data1, data2)

head(data3)
tail(data3)
dim(data3)

## Full join
data3 <- merge(data1, data2)
head(data3)
tail(data3)

## Inner join
data4 <- merge(data1, data2, all=FALSE)
head(data4)
tail(data4)

## Windows
data5 <- window(data, start=as.Date("2020-01-01"), end=as.Date("2020-12-31"))

length(data5)
head(data5)

## lags
data6 <- lag(data, -1)
head(merge(data, data6))
tail(merge(data, data6))

## which package is on the top?
data6 <- stats::lag(data, -1)
head(merge(data, data6))

library(dplyr)
data6 <- dplyr::lag(data, 1)
head(merge(data, data6))

## leads
data7 <- lag(data, 1)   # dplyr on the top...
head(merge(data, data7))
data7 <- stats::lag(data, 1)
head(merge(data, data7))
tail(merge(data, data7))

## differences
data8 <- diff(data)
head(data8)

## Log-returns
data9 <- diff(log(data))
head(data9)

## simple returns
data10 <- diff(data)/stats::lag(data, -1)
head(data10)

plot(merge(data9, data10), plot.type="single", col = c("black", "red"))


## apply function
tmp <- proc.time()
rollsd1 <- rep(0, length(data) -9)
for(i in 1:(length(data) -9)){
  rollsd1[i] <- sd(data[i:(i+9)])
}
proc.time() - tmp

head(rollsd1)

tmp <- proc.time()
rollsd2 <- rollapply(data, width=10, sd, by=1)
proc.time() - tmp

head(rollsd2)

## frequency conversion - xts package
require(xts)

rollsd3 <- apply.weekly(data, sd)
head(rollsd3)

weeklydata <- apply.weekly(data, last)
head(weeklydata)

## A file with WIG data
filename <- "https://stooq.pl/q/d/l/?s=wig&i=d"
x        <- read.csv(filename)
y <- x[,c(1,5)]
colnames(y) <- c("Dates","Price")
write.csv(y,"Example.csv", row.names = FALSE)

P       <- zoo(y$Price, order.by = as.Date(y$Dates))
head(P)
plot(P)
r       <- diff(log(P))


