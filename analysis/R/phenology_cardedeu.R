# Load libraries
library(RColorBrewer) # Cynthia Brewer's Color tables
library(stats4)       # includes mle
library(MASS)         # includes fitdistr
library(gamlss.dist)  # includes generalized gamma distribution
library(gamlss)       # includes generalized gamma fitting

source("dayofyear.R")

treetype <- "Quercus"
treetype <- "Acer"
phase <- "fl"

datadir <- "../../data/cardedeu/"
datafile <- "cardedeu_all.dat"
plotdir <- "/Users/stockli/publications/pheno_light/figures/"

use.pdf <- F
missing <- "32767"

year.start <- 1952
year.end <- 2000
nyear <- year.end - year.start + 1

# open file
con <- file(paste(datadir,datafile,sep=""),"r")

# read dates
temp <- read.table(con,header=T, stringsAsFactors=F)

# close file
close(con)

# generate doy by station and year
#doy <- unlist(subset(temp,select=paste(phase,"_",treetype,sep="")))
doy <- unlist(subset(temp,select=colnames(temp)[which(grepl(paste(phase,"_",sep=""),colnames(temp)))]))



# define colors (black and white, plus 8 colors)
brewercolors <- c("#000000",brewer.pal(8,"Dark2"),"#FFFFFF") #
plotfile <- "sos_histogram"
width <- 6.0
height <- 6.0

plotfile <- paste(phase,"_",treetype,"_sos_histogram_pheno_cardedeu",sep="")

if (use.pdf) {
  pdf(file=paste(plotdir,plotfile,'.pdf',sep=""), onefile=FALSE,
      height=height, width=width, pointsize=12)
} else {
  X11(height=height,width=width) # choose appropriate aspect ratio on device
}

#capitalize first letter
#treetype2 <- gsub(pattern="^(\\w)",replacement="\\U\\1",x=treetype,perl=TRUE)
maintitle <- paste(treetype," (",phase,")",sep="")
ytitle <- "Occurence PDF"
xtitle <- "Start of Season"

x0 <- ISOdatetime(2008,1,1,0,0,0,tz="GMT")
dx <- 86400.

xlim <- c(min(doy,na.rm=T)-5,max(doy,na.rm=T)+5)
#xlim <- c(70,170)
ylim <- c(0.0,0.1)

plot(x0 + xlim*dx,ylim,main = maintitle, ylab = ytitle, xlab = xtitle, xlim = x0+xlim*dx, ylim = ylim, 
     type='n', col=brewercolors[1],lwd=2, xaxt="n")
axis.POSIXct(1, at=seq(x0+xlim[1]*dx, x0+xlim[2]*dx, by="week"), format="%b %d")

par.save <- par(no.readonly=T)

# create histogram of modeled SOS

z <- doy[which(!is.na(doy))]
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
