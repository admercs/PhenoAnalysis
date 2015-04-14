###################################
# Plot Statistical Spring Plant

###################################
# Log

# 16/2/2007: "Index.2005.JuergUncert.dat" gespeichert als "PhenoSpringIndex_RutishauserEtAl_inReview.dat" und an Annette geschickt
# 10/7/2008: Updated Plot bis 2007 gemacht und StandAloneVersion zum weitergeben aufbereitet
# 27/1/2012: Anpassungen an Update bis 2010 für Retos Paper

# Daten







###################################
# Figure 6
# Recon 1702-2005

# Directory

indir <- "/Users/stockli/modelfarm/tower_data/selected/phenology/Swiss_Lowland/"
output <- "StatSpringPlant_RutishauserEtAl_1702_2007"

# Index-Daten
Uncert <- read.table(paste(indir, "RutishauserEtAl_2007_StatSpringPlantReconstruction_JGR_updated_2010.dat",sep=""),header=T,sep="\t")

Spring_1965_2002 <- Uncert[264:301,]
#Uncert <- export.index
y.lim<-105:135



pdf(paste(output,".pdf",sep=""),width = 11, height = 8)

plot(Uncert[,1],Uncert[,2],type="l",ylim=range(y.lim),ylab="Day of Year",xlab="Year",
main="Swiss Statistical Spring Plant 1702-2010")
# lines(Uncert[,1],filter15(Uncert[,2]),lwd=2)

lines(Uncert[,1],Uncert[,3],lwd=3,col="blue")

lines(Uncert[,1],Uncert[,4],col="blue")
lines(Uncert[,1],Uncert[,5],col="blue")

abline(h=117.6)
abline(h=117.6+10.25,lty=2)
abline(h=117.6-10.25,lty=2)

legend(1700,108,
c("Statistical Spring Plant","20-year triangular filter", "Rutishauser et al. (2007) JGR 112: G04016"),
col=c("black","blue",NA),
lty=c(1,1,NA),
lwd=c(1.5,3,NA), bty = "n")

#abline(v=1951,lty=3)

dev.off()

