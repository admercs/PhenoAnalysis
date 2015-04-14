# Load libraries
library(RColorBrewer) # Cynthia Brewer's Color tables
library(stats4)       # includes mle
library(MASS)         # includes fitdistr
library(gamlss.dist)  # includes generalized gamma distribution
library(gamlss)       # includes generalized gamma fitting

source("dayofyear.R")
source("monthdayofyear.R")
source("read_phenoanalysis.R")


#treetype <- "bergahorn"
#treetype <- "edelkastanie"
#treetype <- "haengebirke"
#treetype <- "haselstrauch"
#treetype <- "laerche"
#treetype <- "robinie"
#treetype <- "sommerlinde"
#treetype <- "winterlinde"
#treetype <- "vogelbeere"
treetype <- "buche"

alt.min <- 300
#alt.max <- 500
#alt.min <- 500
alt.max <- 700
#alt.min <- 700
#alt.max <- 900
#alt.min <- 900
#alt.max <- 1100
#alt.min <- 1100
#alt.max <- 1300
#alt.min <- 1300
#alt.max <- 1500
#alt.min <- 1500
#alt.max <- 1700

datadir <- "../../data/phenology/dwh/"
plotdir <- "/Users/stockli/publications/pheno_light/figures/"

use.pdf <- F
missing <- "32767"

tsum.thres <- 273.0+4.0

year.start <- 1960
year.end <- 2011
nyears <- year.end - year.start + 1

# open file
con <- file(paste(datadir,treetype,"_",year.start,"-",year.end,".dat",sep=""),"r")

# read header
for (h in 1:5) {
  temp <- readLines(con,n=1)
}

# read station data
id <- c()
name <- c()
alt <- c()
chx <- c()
chy <- c()

while (substr(temp,1,6) != "number") {
  id <- c(id,as.integer(substr(temp,1,7)))
  name <- c(name,substr(temp,8,31))
  alt <- c(alt,as.integer(substr(temp,32,35)))
  chx <- c(chx,as.integer(substr(temp,46,51)))
  chy <- c(chy,as.integer(substr(temp,54,59)))
  temp <- readLines(con,n=1)
}

ns <- length(id)

# read header
for (h in 1:2) {
  temp <- readLines(con,n=1)
}

# read dates
temp <- read.table(con,header=T, stringsAsFactors=F)

# close file
close(con)

# generate doy by station and year
doy <- array(NA,c(nyears,ns))

for (s in 1:ns) { 
  idx <- which(temp$STA == id[s])
  
  for (i in 1:length(idx)) {
    d <- temp[[7]][idx[i]]

    if (d != missing) {
      day <- as.integer(substr(d,1,2))
      month <- as.integer(substr(d,4,5))
      year <- as.integer(substr(d,7,10))
      
      doy[year-year.start+1,s] <- dayofyear(year,month,day)
    }
  }
}

tair.table <- read_phenoanalysis("/Users/stockli/phenoanalysis-output","Prediction","Swiss-Lowland",
                                 "TMEAN",7.5,8.0,47,47.5,0.5,0.5,
                                 date.start=ISOdatetime(year.start,1,1,0,0,0,tz="GMT"),
                                 date.end=ISOdatetime(year.end,12,31,23,59,59,tz="GMT"))
# extract SOS from model
for (year in year.start:year.end) {
  if (year == year.start) {
    temptime <- ISOdatetime(year,1,1,0,0,0,tz="GMT")
  } else {
    temptime <- c(temptime,ISOdatetime(year,1,1,0,0,0,tz="GMT"))
  }
}
nexp <- 1
tair <- vector(mode="list",length=(1+nexp))
tair[[1]] <- temptime
tsum <- array(NA,c(nyears,ns))
tpast <- array(NA,c(nyears,ns))

for (e in 1:nexp) {
  tair[[1+e]] <- rep(NA,nyears)
  for (year in year.start:year.end) {
    tidx <- which((tair.table[[1]]>=ISOdatetime(year,2,1,0,0,0,tz="GMT")) &
                  (tair.table[[1]]<=ISOdatetime(year,4,30,23,59,59,tz="GMT")))
    tair[[1+e]][year-year.start+1] <- mean(tair.table[[1+e]][tidx],na.rm=T)

    for (s in 1:ns) {
      if (!is.na(doy[(year-year.start+1),s])) {
        temp <- monthdayofyear(year,doy[(year-year.start+1),s])      
        tidx <- which((tair.table[[1]]>=ISOdatetime(year,1,1,0,0,0,tz="GMT")) &
                      (tair.table[[1]]<=ISOdatetime(year,temp[2],temp[1],23,59,59,tz="GMT")))
        tsum[(year-year.start+1),s] <- sum(pmax(tair.table[[1+e]][tidx]-tsum.thres,0),na.rm=T)
        tidx <- which((tair.table[[1]]>=(ISOdatetime(year,temp[2],temp[1],23,59,59,tz="GMT")-60.*86400.)) &
                      (tair.table[[1]]<=ISOdatetime(year,temp[2],temp[1],23,59,59,tz="GMT")-10.*86400.))
        tpast[(year-year.start+1),s] <- mean(tair.table[[1+e]][tidx],na.rm=T)
      }
    }
  }
}

trend <- array(NA,c(ns,2))
for (s in 1:ns) {
  if ((length(which(!is.na(doy[,s]))) > 20) && (alt[s] > alt.min) && (alt[s] < alt.max)) {
    doythres <- 0.5*(min(doy[,s],na.rm=T) + max(doy[,s],na.rm=T))
    
    idx <- which(doy[,s] > doythres)
    y <- doy[idx,s]
    x <- tair[[1+e]][idx]
    x <- tpast[idx,s]
    tr <- lm(y ~ x, na.action="na.exclude")
    coef <- coefficients(tr)
    pv <- summary(tr)$coefficients[2,4]
    print(paste("LATE: Observed Trend (days/K) ",format(coef[[2]],digits=3)," p-value: ",format(pv,digits=3),sep=""))

    if (pv < 0.1) {
      trend[s,1] <- coef[[2]]
    }
    
    idx <- which(doy[,s] <= doythres)

    y <- doy[idx,s]
    x <- tair[[1+e]][idx]
    x <- tpast[idx,s]
    tr <- lm(y ~ x, na.action="na.exclude")
    coef <- coefficients(tr)
    pv <- summary(tr)$coefficients[2,4]
    print(paste("EARLY: Observed Trend (days/K) ",format(coef[[2]],digits=3)," p-value: ",format(pv,digits=3),sep=""))
    if (pv < 0.1) {
      trend[s,2] <- coef[[2]]
    }
  }
}
print(paste("Late:  Observed Trend (days/K): ",format(mean(trend[,1],na.rm=T),digits=3)))
print(paste("Early: Observed Trend (days/K): ",format(mean(trend[,2],na.rm=T),digits=3)))

browser()

plot(tpast[,1],doy[ ,1],xlim=c(min(tpast,na.rm=T),max(tpast,na.rm=T)),ylim=c(min(doy,na.rm=T),max(doy,na.rm=T)))
plot(tsum[ ,1],doy[ ,1],xlim=c(min(tsum,na.rm=T),max(tsum,na.rm=T)),ylim=c(min(doy,na.rm=T),max(doy,na.rm=T)))
plot(tair[[2]],doy[ ,1],xlim=c(min(tair[[2]],na.rm=T),max(tair[[2]],na.rm=T)),ylim=c(min(doy,na.rm=T),max(doy,na.rm=T)))


plot(tair[[2]],doy[ ,1],ylim=c(min(doy[,1],na.rm=T),max(doy[,1],na.rm=T)))
for (s in 2:ns) {
  points(tair[[2]],doy[ ,s])
}
#r <- cor(tair[[2]],obs[[2]],method="spearman",use="complete.obs")
#print(paste("OBS: Spearman Correlation: ", format(r,digits=3),sep=""))

browser()


# define colors (black and white, plus 8 colors)
brewercolors <- c("#000000",brewer.pal(8,"Dark2"),"#FFFFFF") #
plotfile <- "sos_histogram"
width <- 6.0
height <- 6.0

plotfile <- paste(treetype,"_sos_histogram_pheno_swiss_",alt.min,"m-",alt.max,"m",sep="")

if (use.pdf) {
  pdf(file=paste(plotdir,plotfile,'.pdf',sep=""), onefile=FALSE,
      height=height, width=width, pointsize=12)
} else {
  X11(height=height,width=width) # choose appropriate aspect ratio on device
}

#capitalize first letter
treetype2 <- gsub(pattern="^(\\w)",replacement="\\U\\1",x=treetype,perl=TRUE)
maintitle <- paste(treetype2," (",alt.min,"m-",alt.max,"m)",sep="")
ytitle <- "Occurence PDF"
xtitle <- "Start of Season"

x0 <- ISOdatetime(2008,1,1,0,0,0,tz="GMT")
dx <- 86400.

xlim <- c(min(doy,na.rm=T)-5,max(doy,na.rm=T)+5)
xlim <- c(70,170)
ylim <- c(0.0,0.1)

plot(x0 + xlim*dx,ylim,main = maintitle, ylab = ytitle, xlab = xtitle, xlim = x0+xlim*dx, ylim = ylim, 
     type='n', col=brewercolors[1],lwd=2, xaxt="n")
axis.POSIXct(1, at=seq(x0+xlim[1]*dx, x0+xlim[2]*dx, by="week"), format="%b %d")

par.save <- par(no.readonly=T)

alt.all <- array(NA,c(nyears,ns))
chx.all <- array(NA,c(nyears,ns))
chy.all <- array(NA,c(nyears,ns))
id.all <- array(NA,c(nyears,ns))
for (s in 1:ns) {
  alt.all[,s] <- alt[s]
  chx.all[,s] <- chx[s]
  chy.all[,s] <- chy[s]
  id.all[,s] <- id[s]
}

# create histogram of modeled SOS

#  z <- doy[which(!is.na(doy) & (alt.all < 600) & (chx.all < 600000) & (chy.all > 200000))]
z <- doy[which(!is.na(doy) & (alt.all >= alt.min) & (alt.all < alt.max))]

n <- length(z)

if (n > 10) {
  
  mids <- seq(min(z),max(z),by=1)
  breaks <- seq(min(z)-0.5,max(z)+0.5,by=1)
  h = hist(z,plot=F,breaks=breaks)
  y <- h$counts
  
  # draw histogram
  for (b in 1:length(mids)) {
    if (y[b] > 0) {
      lines(x0+(c(mids[b],mids[b]))*dx,c(0,y[b]/n),col=paste(brewercolors[2],"50",sep=""),lwd=4,lend=2)
    }
  }
  
  # fit (generalized) gamma curves
  ggfit <- gamlssML(z,family=GG)
  
  # check whether a normal distribution can be used
  stest <- shapiro.test(z)
  
  # calculate Skewness
  skewness = sum((z-mean(z))^3/sqrt(var(z))^3,na.rm=T)/n
  
  print(paste("n: ",n))
  print(paste("Shapiro-Wilks Normality Test p-value: ",stest$p.value))
  print(paste("Skewness: ",skewness))
  print(paste("Stddev:   ",sd(z)))
  print(paste("Gamma (mu, sigma, nu): ",ggfit$mu, ggfit$sigma, ggfit$nu))
  
  # draw fitted curve
  x <- seq(xlim[1],xlim[2],by=0.1)
  y <- max(h$counts)/max(h$density)*dGG(x,mu = ggfit$mu, sigma = ggfit$sigma, nu = ggfit$nu)
  lines(x0+x*dx,y/n,col=paste(brewercolors[2],"50",sep=""),lwd=4)

  text <- c(paste("n = ",formatC(n,format="f",digits=0),sep=""),
            paste("Shapiro-Wilks Normality p = ",formatC(stest$p.value,width=5,format="f",digits=4),sep=""),
            paste("Skewness = ",formatC(skewness,width=5,format="f",digits=3),sep=""),
            paste("Mean = ",formatC(mean(z),width=5,format="f",digits=2),sep=""),
            paste("SD = ",formatC(sd(z),width=5,format="f",digits=2),sep=""))

  legend("topleft",text,bty="n",lty=0)
 if (use.pdf) {
   dev.off()
 }
 
}
