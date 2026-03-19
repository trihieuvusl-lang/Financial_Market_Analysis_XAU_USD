require(zoo)
setwd("C:/Users/Vuxhieeu/Downloads/Rcodes")
# Choose Your fund from the list
# browseURL("https://stooq.pl/q/?s=xauusd")

#fundname    <- "1018.n"   # <-- example 
fundname    <- "xauusd"       # <-- this line needs to be adjusted

filename    <- paste0(fundname,".csv")

# Creating CSV file (set TRUE only at the first time and then return to FALSE)
CreFile = FALSE
if(CreFile){
  filepath    <- paste0("https://stooq.pl/q/d/l/?s=",fundname,"&i=d")
  x           <- read.csv2(file=filepath, head=TRUE, sep=",", dec=".")
  y           <- x[,c(1,5)]
  colnames(y) <- c("Dates","Price")
  ind         <- as.Date(y$Dates)
  N           <- max(as.numeric(diff(tail(ind,2000))))
  # checks if at least 2000 observations and there is no observation gap longer that 7 days
  if(dim(y)[1]>2000 & N<8){
    write.csv(y,filename, row.names = FALSE)
  }
  # check if you have .csv file in your working directory 
}

# Downloading the data
y     <- read.csv2(file=filename , head=TRUE, sep=",", dec=".")
P     <- zoo(y$Price, order.by = as.Date(y$Dates))
r     <- diff(log(P))

