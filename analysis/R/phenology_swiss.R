# This code generates figures for the Swiss phenology climate sensitivity experiments
# The code depends on the following data files to be present in datadir:
# light_fac.dat
# moist_fac.dat
# sos.dat
# temp_fac.dat

# how to run:
# adjust datadir and plotdir (below)
# select kind of plot (plottype = 1..3
# select whether you would like to have screen plots (use.pdf = F) or pdf plots (use.pdf = T)

# kill all zombie variables (test if needed)
# rm(list = ls(all = TRUE))

# Load libraries
library(RColorBrewer) # Cynthia Brewer's Color tables
library(stats4)       # includes mle
library(MASS)         # includes fitdistr
library(gamlss.dist)  # includes generalized gamma distribution
library(gamlss)       # includes generalized gamma fitting

source("read_phenoanalysis.R")

# set X11 instead of Cairo plotting
#X11.options(type="Xlib")

sitetype <- "swiss"
sitenum <- 7
#sitetype <- "harvard"
#sitenum <- 3

experiments <- c("Assimilation")
experiments <- c("Prediction")
#experiments <- c("Swiss-Prediction-0K")

# directories
modeldir <- "/Users/stockli/phenoanalysis-output"
plotdir <- "/Users/stockli/Desktop/"

# flags and parameters
plottype <- 1 # plot: 1) time series 2) histogram 3) limitations
use.pdf <- F  # pdf or window
use.dwh <- T  # read phenology from DWH
use.png <- T  # convert pdf to png (requires use.pdf to be true)
use.qq <- F   # plot QQ instead of histogram in plot 2
stestlim <- 0.05  # limiting p-value of shapiro-wilks test (lower than the limit: not a normal distribution)

if (sitetype == "swiss") {
  year.start <- 1960
  year.end <- 2010
} else {
  year.start <- 1990
  year.end <- 2012
}

# CODE STARTS HERE

nexp <- length(experiments)
nyears <- year.end - year.start + 1
years <- year.start:year.end
date.start <- ISOdatetime(year.start,1,1,0,0,0,tz="GMT")
date.end <- ISOdatetime(year.end,12,31,23,59,59,tz="GMT")

# in case no site name is given: read all sites of given sitetype from table
if (sitetype != "none") {
  coord <- read.table(paste("../../data/sites/",sitetype,".dat",sep=""),
                      header=T,sep=",",stringsAsFactors=FALSE)
  # remove tabs
  coord[[2]] <- gsub("\t","",coord[[2]])
  coord[[3]] <- gsub("\t","",coord[[3]])
  coord[[10]] <- gsub("\t","",coord[[10]])
  coord[[11]] <- gsub("\t","",coord[[11]])
  coord[[12]] <- gsub("\t","",coord[[12]])
  coord[[13]] <- gsub("\t","",coord[[13]])

  # extract chosen site
  idx <- which(coord[[1]] == sitenum)
  sitename <- coord[[3]][idx]
  xmin <- coord[[4]][idx]
  xmax <- coord[[5]][idx]
  ymin <- coord[[6]][idx]
  ymax <- coord[[7]][idx]
  dx <- coord[[8]][idx]
  dy <- coord[[9]][idx]
  projection <- coord[[11]][idx]
  projparam <- coord[[12]][idx]
  hgt <- coord[[13]][idx]
}

# generate observation list
# 1. element: time (POSIXct)
# 2. element: SOS (DoY)
# 3. element: SOS uncertaity (days)
obs <- vector(mode="list",length=(3))
for (year in year.start:year.end) {
  if (year == year.start) {
    temptime <- ISOdatetime(year,1,1,0,0,0,tz="GMT")
  } else {
    temptime <- c(temptime,ISOdatetime(year,1,1,0,0,0,tz="GMT"))
  }
}
      
if (sitetype == "swiss") {
  delta <- 5.4 # observation uncertainty (days) # This Rutishauser, personal communication
  
  # read sos observations (Spring Index by This Rutishauser 2007, extended to 2010)
  obs.table <- read.table("../../data/phenology/Swiss_Lowland/RutishauserEtAl_2007_StatSpringPlantReconstruction_JGR_updated_2010.dat",header=T)
  has.obs <- T
  
  # extract SOS from observations
  for (year in year.start:year.end) {
    tidx <- which(obs.table[[1]]==year)
    if (year == year.start) {
      if (length(tidx) == 1) {
        tempsos <- obs.table[[2]][tidx]
        tempvar <- delta
      } else {
        tempsos <- NA
        tempvar <- NA
      }
    } else {
      if (length(tidx) == 1) {
        tempsos <- c(tempsos,obs.table[[2]][tidx])
        tempvar <- c(tempvar,delta)
      } else {
        tempsos <- c(tempsos,NA)
        tempvar <- c(tempvar,NA)
      }
    }
  }

  laithres = 0.45
  
} else {
  obs.table <- read.csv("../../data/phenology/Harvard_Forest/hf003-03-spring.csv",
                          header=T,sep=",",stringsAsFactors=FALSE)
  has.obs <- T
  thres <- 50
  ntrees <- 4
#  species <- c('QURU','ACRU','BELE','TSCA','QUAL','QUVE')
  species <- c('QURU','ACRU','QUAL')
  phases <- c('BBRK')  
  nspecies <- length(species)
  nphases <- length(phases)
  nval <- length(obs.table$DATE)
  
  val <- array(data=NA,dim=c(nyears,nphases,nspecies,ntrees))
  
  phaseidx <- vector(length=nphases)
  for (p in 1:nphases) {
    phaseidx[p] <- which(names(obs.table) == phases[p])
  }

  # get SOS dates for individual trees, species and phases by year
  for (v in 1:nval) {
    year <- as.numeric(substr(obs.table$DATE[v],1,4))
    spec <- substr(obs.table$TREEID[v],1,4)
    y <- year - year.start + 1
    s <- which(species == spec)
    t <- as.numeric(substr(obs.table$TREEID[v],6,7))
    if ((y >= 1) & (y <= nyears) & (length(s) != 0) & (t <= ntrees)) {
      for (p in 1:nphases) {
        if ((is.na(val[y,p,s,t])) & (!is.na(obs.table[[phaseidx[p]]][v])) &
            (obs.table[[phaseidx[p]]][v] >= thres) ) {
          val[y,p,s,t] <- obs.table$JULIAN[v]
        }
      }
    }
  }

  # average over trees, species and phases
  tempsos <- apply(val,c(1),mean,na.rm=T)
  tempvar <- apply(val,c(1),sd,na.rm=T)

  laithres = 0.55
  
}
  
obs[[1]] <- temptime
if (has.obs) {
  obs[[2]] <- tempsos
  obs[[3]] <- tempvar
} else {
  obs[[2]] <- rep(NA,nyears)
  obs[[3]] <- rep(NA,nyears)
}

# read sos and limitation factors from model
mod.table <- read_phenoanalysis(modeldir,experiments,sitename,"SOS",xmin,xmax,ymin,ymax,dx,dy,
                                date.start=date.start, date.end=date.end)

tair.table <- read_phenoanalysis(modeldir,experiments,sitename,"TMIN",xmin,xmax,ymin,ymax,dx,dy,
                                date.start=date.start, date.end=date.end)

# extract SOS from model
mod <- vector(mode="list",length=(1+nexp))
mod[[1]] <- temptime
tair <- vector(mode="list",length=(1+nexp))
tair[[1]] <- temptime
for (e in 1:nexp) {
  laimin <- min(mod.table[[1+e]])
  laimax <- max(mod.table[[1+e]])
  threshold <- laithres *(laimax - laimin) + laimin
  tmin <- 50

  mod[[1+e]] <- rep(NA,nyears)
  tair[[1+e]] <- rep(NA,nyears)
  for (year in year.start:year.end) {
    tidx <- which((mod.table[[1]]>=ISOdatetime(year,1,1,0,0,0,tz="GMT")) &
                  (mod.table[[1]]<=ISOdatetime(year,12,31,23,59,59,tz="GMT")))
    tval <- which(mod.table[[1+e]][tidx[tmin:length(tidx)]] >= threshold)

    mod[[1+e]][year-year.start+1] <- as.numeric(format(mod.table[[1]][tidx[tval[1]]+tmin],"%j",TZ="GMT"))

    tidx <- which((tair.table[[1]]>=ISOdatetime(year,2,1,0,0,0,tz="GMT")) &
                  (tair.table[[1]]<=ISOdatetime(year,4,30,23,59,59,tz="GMT")))
    tair[[1+e]][year-year.start+1] <- mean(tair.table[[1+e]][tidx],na.rm=T)
  }
}

doythres <- 0.5*(min(mod[[1+e]],na.rm=T) + max(mod[[1+e]],na.rm=T))
    
idx <- which(mod[[1+e]] > doythres)
y <- mod[[1+e]][idx]
x <- tair[[1+e]][idx]
tr <- lm(y ~ x, na.action="na.exclude")
coef <- coefficients(tr)
pv <- summary(tr)$coefficients[2,4]
print(paste("LATE: Observed Trend (days/K) ",format(coef[[2]],digits=3)," p-value: ",format(pv,digits=3),sep=""))
    
idx <- which(mod[[1+e]] <= doythres)
y <- mod[[1+e]][idx]
x <- tair[[1+e]][idx]
tr <- lm(y ~ x, na.action="na.exclude")
coef <- coefficients(tr)
pv <- summary(tr)$coefficients[2,4]
print(paste("EARLY: Observed Trend (days/K) ",format(coef[[2]],digits=3)," p-value: ",format(pv,digits=3),sep=""))

browser()

plot(tair[[2]],obs[[2]])
r <- cor(tair[[2]],obs[[2]],method="spearman",use="complete.obs")
print(paste("OBS: Spearman Correlation: ", format(r,digits=3),sep=""))

plot(tair[[2]],mod[[2]])
r <- cor(tair[[2]],mod[[2]],method="spearman",use="complete.obs")
print(paste("MOD: Spearman Correlation: ", format(r,digits=3),sep=""))

browser()

# define legend text
if (has.obs) {
  legendtext <- c('OBS',experiments)
} else {
  legendtext <- experiments
}

# define colors (black and white, plus 8 colors)
brewercolors <- c("#000000",brewer.pal(8,"Dark2"),"#FFFFFF") #

if (plottype == 1) {
  plotfile <- "sos_timeseries"
  width <- 10.0
  height <- 6.0
} else if (plottype == 2) { 
  if (use.qq) {
    plotfile <- "sos_qqplot"
    rows <- 2
    cols <- 3
    width <- 12.0
    height <- 8.0
  } else {
    plotfile <- "sos_histogram"
    width <- 6.0
    height <- 6.0
  }
} else {
  plotfile <- "sos_limitation"
  width <- 6.0
  height <- 6.0
}

plotfile <- paste(plotfile,'_',sitename,sep='')

if (use.pdf) {
  pdf(file=paste(plotdir,plotfile,'.pdf',sep=""), onefile=FALSE,
      height=height, width=width, pointsize=12)
} else {
  X11(height=height,width=width) # choose appropriate aspect ratio on device
}

# simple time series with start of season includes uncertainty
if (plottype == 1) {
  
  x <- obs[[1]]
  y <- obs[[2]]
  ydelta <- obs[[3]]
  
  ytitle <- "Start of Season"
  xtitle <- "Years"

  y0 <- ISOdatetime(2008,1,1,0,0,0,tz="GMT")
  dy <- 86400.
  
  ylim <- c(min(c(obs[[2]],mod[[2]]),na.rm=T)-5,max(c(obs[[2]],mod[[2]]),na.rm=T)+5)
  
  plot(x,y0 + y*dy,main = sitename, ylab = ytitle, xlab = xtitle, ylim = y0+ylim*dy,
       type='l', col=brewercolors[1],lwd=4, yaxt="n")

  axis.POSIXct(2, at=seq(y0+ylim[1]*dy, y0+ylim[2]*dy, by="week"), format="%b %d")

  # uncertainty from This Rutishauser's time series
  if (has.obs) {
    polygon(c(x,rev(x)),y0+c(y-ydelta,rev(y+ydelta))*dy,
            col=paste(brewercolors[1],"20",sep=""), border=NA)
    
    print(legendtext[1])
    
    trend <- lm(y ~ years, na.action="na.exclude")
    coef <- coefficients(trend)
    pval <- summary(trend)$coefficients[2,4]
    print(paste("Observed Trend (days/year) ",format(coef[[2]],digits=3)," p-value: ",format(pval,digits=3),sep=""))
    legendtext[1] <- paste(legendtext[1]," (trend (d/y) = ",formatC(coef[[2]],digits=2),")",sep=" ")
  }
  
  for (e in 1:nexp) {

    # maximum model SOS uncertainty given a LAI uncertainty of 0.35 m2/m2
    # this is a very conservative number. The real number is around 2-3 days
    delta <- 5.0
    y <- mod[[1+e]]
    points(x,y0+y*dy,type='l',col=brewercolors[e+1],lwd=4)
    polygon(c(x,rev(x)),y0+c(y-delta,rev(y+delta))*dy,
            col=paste(brewercolors[e+1],"20",sep=""), border=NA)

    print(legendtext[e+1])

    if (has.obs) {
      r <- cor(obs[[2]],y,method="spearman",use="complete.obs")
      print(paste("Spearman Correlation: ",
                  format(r,digits=3),sep=""))
#      print(t.test(obs[[2]],y,alternative="two.sided",paired=T,conf.level=0.95))
#      print(wilcox.test(obs[[2]],y,alternative="two.sided",paired=T,conf.level=0.95))
    }

    trend <- lm(y ~ years, na.action="na.exclude")
    coef <- coefficients(trend)
    pval <- summary(trend)$coefficients[2,4]
    print(paste("Modeled Trend (days/year) ",format(coef[[2]],digits=3)," p-value: ",format(pval,digits=3),sep=""))

    legendtext[1+e] <- paste(legendtext[1+e]," (trend (d/y) = ",formatC(coef[[2]],digits=2)," R = ",formatC(r,digits=2),")",sep=" ")

  }

  legend("topright",legendtext,bty="n",text.col=brewercolors[1:(nexp+1)])  
  
}

browser()

# histogram of occurence for a specific start of season date
if (plottype == 2) {

  if (use.qq) par(mfrow=c(rows,cols))
  
  ytitle <- "Number of Occurences"
  xtitle <- "Start of Season"

  x0 <- ISOdatetime(2008,1,1,0,0,0,tz="GMT")
  dx <- 86400.

  xlim <- c(min(c(obs,mod),na.rm=T)-5,max(c(obs,mod),na.rm=T)+5)
#  xlim <- c(60,140)
  
#  ylim <- c(0,30)
  ylim <- c(0,20)

  plot(x0 + xlim*dx,ylim,main = sitename, ylab = ytitle, xlab = xtitle, xlim = x0+xlim*dx, ylim = ylim, 
       type='n', col=brewercolors[1],lwd=2, xaxt="n")
  axis.POSIXct(1, at=seq(x0+xlim[1]*dx, x0+xlim[2]*dx, by="week"), format="%b %d")

  par.save <- par(no.readonly=T)
  
  # create histogram of observed SOS (This Rutishauser's data)
  if (has.obs) {
    mids <- seq(min(obs,na.rm=T),max(obs,na.rm=T),by=1)
    breaks <- seq(min(obs,na.rm=T)-0.5,max(obs,na.rm=T)+0.5,by=1)
    h = hist(obs,plot=F,breaks=breaks)

    # draw histogram
    x <- mids
    y <- h$counts

    for (b in 1:length(mids)) {
      if (y[b] > 0) {
        lines(x0+c(mids[b],mids[b])*dx,c(0,y[b]),col=paste(brewercolors[1],"50",sep=""),lwd=6,lend=2)
      }
    }

    gmin <- min(obs,na.rm=T)-1
    wmin <- gmin
  
    # fit normal, (generalized) gamma and weibull curves
    nfit <- fitdistr(obs,"normal")
    gfit <- fitdistr(obs-gmin,dgamma, list(shape = 1, rate = 0.1), lower = 0.001)
    wfit <- fitdistr(obs-gmin,"weibull", list(shape = 1, scale = 1), lower = 0.001)
    ggfit <- gamlssML(obs,family=GG)
    
    # check whether a normal distribution can be used
    stest <- shapiro.test(obs)

    # calculate Skewness
    skewness = sum((obs-mean(obs,na.rm=T))^3/sqrt(var(obs,na.rm=T))^3,na.rm=T)/length(obs)

    print(legendtext[1])
    print(paste("Shapiro-Wilks Normality Test p-value: ",stest$p.value))
    print(paste("Skewness: ",skewness))
    print(paste("Gamma (mu, sigma, nu): ",ggfit$mu, ggfit$sigma, ggfit$nu))
  
    # draw fitted curve
    x <- seq(xlim[1],xlim[2],by=0.1)

    y <- max(h$counts)/max(h$density)*dGG(x,mu = ggfit$mu, sigma = ggfit$sigma, nu = ggfit$nu)

    lines(x0+x*dx,y,col=brewercolors[1],lwd=4)
    
    # Q-Q plot
    if (use.qq) {
      legend("topleft","a)",bty="n",text.col=brewercolors[1])

      # emprical quantiles
      yq <- obs

      qlim <- c(min(yq,na.rm=T)-5,max(yq,na.rm=T)+5)
    
      par(mfg=c(ceiling((2)/cols),((1) %% cols)+1))

      plot(0,0,xlab="Theoretical Quantiles",ylab="Empirical Quantiles",
           main=paste(legendtext[1]," (Q-Q Plot)",sep=""),type="n",
           xlim=qlim,ylim=qlim)

      # theoretical quantiles for normal fit
      f <- qnorm(ppoints(yq), mean=nfit$estimate[[1]], sd=nfit$estimate[[2]])
      xq <- f[order(order(yq))]
      points(xq,yq+0.25,pch=24,col=brewercolors[1], cex=1.5)

      # theoretical quantlies for gamma fit
#      f <- qgamma(ppoints(yq), shape=gfit$estimate[[1]], rate=gfit$estimate[[2]])+gmin
#      f <- qweibull(ppoints(yq), shape=wfit$estimate[[1]], scale=wfit$estimate[[2]])+wmin
      f <- qGG(ppoints(yq),mu = ggfit$mu, sigma = ggfit$sigma, nu = ggfit$nu)
      xq <- f[order(order(yq))]
      points(xq,yq,pch=21,col=brewercolors[1], cex=1.5)
    
      lines(qlim,qlim,lwd=2)
      legend("bottomright",c("gamma","normal"),bty="n",pch=c(21,24),text.col=brewercolors[1])  

      legend("topleft","b)",bty="n",text.col=brewercolors[1])
    }

  }
  
  for (e in 1:nexp) {
    # create histogram of modeled SOS
    mids <- seq(min(mod[,e],na.rm=T),max(mod[,e],na.rm=T),by=1)
    breaks <- seq(min(mod[,e],na.rm=T)-0.5,max(mod[,e],na.rm=T)+0.5,by=1)
    h = hist(mod[,e],plot=F,breaks=breaks)
    y <- h$counts

    # draw histogram
    par(par.save)
    if (use.qq) par(mfg=c(1,1))
    for (b in 1:length(mids)) {
      if (y[b] > 0) {
        lines(x0+(c(mids[b],mids[b]))*dx,c(0,y[b]),col=paste(brewercolors[1+e],"50",sep=""),lwd=6,lend=2)
      }
    }

    gmin <- min(mod[,e],na.rm=T)-1
    wmin <- gmin
    
    # fit curve
    nfit <- fitdistr(mod[,e],"normal")
    gfit <- fitdistr(mod[,e]-gmin,"gamma", list(shape = 1, rate = 0.1), lower = 0.001)
    wfit <- fitdistr(mod[,e]-wmin,"weibull", list(shape = 1, scale = 1), lower = 0.001)
    ggfit <- gamlssML(mod[,e],family=GG)

    # check whether a normal distribution can be used
    stest <- shapiro.test(mod[,e])

    # calculate Skewness
    skewness = sum((mod[,e]-mean(mod[,e],na.rm=T))^3/sqrt(var(mod[,e],na.rm=T))^3,na.rm=T)/length(mod[,e])
    
    print(legendtext[e+1])
    print(paste("Shapiro-Wilks Normality Test p-value: ",stest$p.value))
    print(paste("Skewness: ",skewness))
    print(paste("Gamma (mu, sigma, nu): ",ggfit$mu, ggfit$sigma, ggfit$nu))

    # evaluate probability of DOY of 1 and 5 April occurrence
    doy1 <- (unclass(ISOdatetime(2007,4,1,0,0,0,tz="GMT")-ISOdatetime(2007,1,1,0,0,0,tz="GMT")))[[1]]
    doy2 <- (unclass(ISOdatetime(2007,4,5,0,0,0,tz="GMT")-ISOdatetime(2007,1,1,0,0,0,tz="GMT")))[[1]]
    
    p1 <- pGG(doy1,mu = ggfit$mu, sigma = ggfit$sigma, nu = ggfit$nu,lower.tail=F)
    p2 <- pGG(doy2,mu = ggfit$mu, sigma = ggfit$sigma, nu = ggfit$nu,lower.tail=F)
    print(paste("p of 1 April :",p1))
    print(paste("p of 5 April :",p2))

    # draw fitted curve
    x <- seq(xlim[1],xlim[2],by=0.1)

    y <- max(h$counts)/max(h$density)*dGG(x,mu = ggfit$mu, sigma = ggfit$sigma, nu = ggfit$nu)

    lines(x0+x*dx,y,col=brewercolors[1+e],lwd=4)

    # Q-Q plot
    if (use.qq) {
      
      # emprical quantiles
      yq <- mod[,e]

      qlim <- c(min(yq,na.rm=T)-5,max(yq,na.rm=T)+5)

      par(mfg=c(ceiling((e+2)/cols),((e+1) %% cols)+1))
      plot(0,0,xlab="Theoretical Quantiles",ylab="Empirical Quantiles",
           main=paste(legendtext[e+1]," (Q-Q Plot)",sep=""),type="n",
           xlim=qlim,ylim=qlim)

      # theoretical quantiles for normal fit
      f <- qnorm(ppoints(yq), mean=nfit$estimate[[1]], sd=nfit$estimate[[2]])
      xq <- f[order(order(yq))]
      points(xq,yq+0.25,pch=24,col=brewercolors[1+e],cex=1.5)

      # theoretical quantlies for gamma fit
#      f <- qgamma(ppoints(yq), shape=gfit$estimate[[1]], rate=gfit$estimate[[2]])+gmin
#      f <- qweibull(ppoints(yq), shape=wfit$estimate[[1]], scale=wfit$estimate[[2]])+gmin
      f <- qGG(ppoints(yq),mu = ggfit$mu, sigma = ggfit$sigma, nu = ggfit$nu)
      xq <- f[order(order(yq))]
      points(xq,yq,pch=21,col=brewercolors[1+e],cex=1.5)
   
      lines(qlim,qlim,lwd=2)
      legend("bottomright",c("gamma","normal"),bty="n",pch=c(21,24),text.col=brewercolors[1])  

      if (e == 1) legend("topleft","c)",bty="n",text.col=brewercolors[1])
      if (e == 2) legend("topleft","d)",bty="n",text.col=brewercolors[1])
      if (e == 3) legend("topleft","e)",bty="n",text.col=brewercolors[1])
      if (e == 4) legend("topleft","f)",bty="n",text.col=brewercolors[1])
    }
      
  }

  if (use.qq) par(mfg=c(1,1))

  if (sitename != "Switzerland") {
    if (has.obs) {
      legendtext[1] <- paste(legendtext[1],':   mean=',format(mean(obs,na.rm=T),digits=4),' sd=',format(sd(obs,na.rm=T),digits=2),sep='')
    }
    for (e in 1:nexp) {
      legendtext[1+e] <- paste(legendtext[1+e],':   mean=',format(mean(mod[,e],na.rm=T),digits=4),' sd=',format(sd(mod[,e],na.rm=T),digits=2),sep='')
    }
  }

  if (!has.obs) legendtext[1] = ""
  legend("topright",legendtext,bty="n",text.col=brewercolors[1:(nexp+1)])  
  
}

# seasonal mean and interannual variability of light and temperature limitation curves
if (plottype == 3) {

  xlim <- c(ISOdatetime(2007,2,1,0,0,0,tz="GMT"),ISOdatetime(2007,10,31,0,0,0,tz="GMT"))
  ylim <- c(0,1.15)
  xtitle <- "Month"
  ytitle <- "Limitation"

  # only use days 1-365 since leap years have not the same sampling
  
  x <- ISOdatetime(2007,1,1,0,0,0,tz="GMT") + seq(0,(ndays-2))*86400
  y <- rowMeans(matrix(light.fac[[3]],c(ndays,nyears)),na.rm=T)[1:(ndays-1)]

  plot(x,y,main = sitename, ylab = ytitle, xlab = xtitle, xlim = xlim, ylim = ylim, 
       type='n', col=brewercolors[1],lwd=2)


  # plot light limitation factor
  y <- rowMeans(matrix(light.fac[[3]],c(ndays,nyears)),na.rm=T)[1:(ndays-1)]
  lines(x,y,col=brewercolors[1],lwd=4,lty=4)
  
  for (e in 1:nexp) {
    # plot temperature limitation factors
    y <- rowMeans(matrix(temp.fac[[2+e]],c(ndays,nyears)),na.rm=T)[1:(ndays-1)]
    delta <- sqrt(rowMeans(matrix(temp.fac[[2+e]]^2,c(ndays,nyears)),na.rm=T) -
                  rowMeans(matrix(temp.fac[[2+e]],c(ndays,nyears)),na.rm=T)^2)[1:(ndays-1)]
    polygon(c(x,rev(x)),c(y-delta,rev(y+delta)),
            col=paste(brewercolors[1+e],"20",sep=""), border=NA)
    lines(x,y,col=brewercolors[1+e],lwd=4, lty=2)
  }

  for (e in 1:nexp) {
    # plot SOS dates
    y <- rowMeans(matrix(temp.fac[[2+e]],c(ndays,nyears)),na.rm=T)[1:(ndays-1)]
    mm <- mean(mod[,e],na.rm=T)
    ms <- sd(mod[,e],na.rm=T)
    points(ISOdatetime(2007,1,1,0,0,0,tz="GMT") + c(mm,mm)*86400.,c(0,y[mm]),type='l', lend=1,
           col=brewercolors[e+1],lwd=2)
    points(ISOdatetime(2007,1,1,0,0,0,tz="GMT") + c(0,mm)*86400.,c(y[mm],y[mm]),type='l', lend=1,
           col=brewercolors[e+1],lwd=2)
#    points(ISOdatetime(2007,1,1,0,0,0,tz="GMT") + mm*86400.,y[mm],pch=19,
#           col=brewercolors[e+1],lwd=2)
#    polygon(c(x,rev(x)),y0+c(y-delta,rev(y+delta))*dy,
#            col=paste(brewercolors[e+1],"20",sep=""), border=NA)

    # print relative photoperiod limitation
    print(legendtext[e+1])
    light <- matrix(light.fac[[2+e]],c(ndays,nyears))
    temp <- matrix(temp.fac[[2+e]],c(ndays,nyears))
    moist <- matrix(moist.fac[[2+e]],c(ndays,nyears))
    doy <- mod[,e]
    avg <- 0
    rellim <- array(dim=nyears)
    for (y in 1:nyears) {
      rellim[y] <- 100.* mean((1.-light[(doy[y]-avg):doy[y],y])) / mean((1.-light[(doy[y]-avg):doy[y],y]) +
                                                              (1.- temp[(doy[y]-avg):doy[y],y]) +
                                                              (1.-moist[(doy[y]-avg):doy[y],y]))
    }
    print(paste("Relative photoperiod limitation: ",mean(rellim)))

    limindex <-  array(dim=nyears)
    for (y in 1:nyears) {
      limindex[y] <- mean(1.-light[(doy[y]-avg):doy[y],y]) / max(c(mean(1.-temp[(doy[y]-avg):doy[y],y]),0.01))
    }
    print(paste("Index of Photoperiod limitation: ",mean(limindex)))
    
  }

  legend(ISOdatetime(2007,7,1,0,0,0,tz="GMT"),0.6,c("Scenario:",legendtext3),
         bty="n",text.col=brewercolors[1:(nexp+1)])
  legend("topleft",legendtext2,bty="n",text.col=brewercolors[1],lty=c(2,4),lwd=2)

  text(ISOdatetime(2007,5,5,0,0,0,tz="GMT"),0,"Mean SOS Dates", adj=c(0,0.5))

  
  
}

if (use.pdf) {
  ## close PDF file
  dev.off()
  
  ## PDF -> PNG
  if (use.png) {
    pngfile <- paste(plotfile,'.png',sep="")
    cmd <- paste("convert -density 300 -depth 8 -resize 1280 -flatten ",plotdir,plotfile,".pdf ",plotdir,pngfile,sep="")
    system(cmd)
  }
}


