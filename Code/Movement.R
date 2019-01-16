rm(list=ls())

# load libraries ----------------------------------------------------------
library(ggplot2)
library(BayesianTools)

# set workdirectory -------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location (if working in RStudio)
getwd()


# source files ------------------------------------------------------------
source('Movement.FigureScripts.R')

# load data ---------------------------------------------------------------
load('../Input/Movement.RData')

# Figure S3: repellency and expellency ---------------------------------------

transparancy = 0.5
XLIMS = rbind(c(5e-4,4e-3),c(5e-4,4e-3),c(5e-4,4e-3),
              c(0.4,0.75), c(NA,NA), c(NA,NA),c(NA,NA),
              c(0,2e-4), c(0,2e-4),c(0,7e-4),
              c(5e-4,1e-3) )

pdf('FigS3.pdf',width=7.5,height=4.75)
par(mfrow = c(1,2),mai=c(0.3,0.6,0.6,0.1) ,oma=c(3,4,1,1))  

p1 <- density( Repellency[1,], na.rm = TRUE);   p1$y <- p1$y*diff(p1$x)[1]
p2 <- density( Repellency[2,], na.rm = TRUE);  p2$y <- p2$y*diff(p2$x)[1]  
p3 <- density( Repellency[3,], na.rm = TRUE);  p3$y <- p3$y*diff(p3$x)[1]

plot(p1,col='white',main='',xlim=c(-.4,.4), ylim=c(0,7.5e-3),xaxs='i',yaxs='i',las=1)
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
mtext( bquote(bold("A")), side = 3, line = 0.1 ,at =-0.37,cex=1) 
mtext(text=expression('Repellency (' ~rho* ' )'),side=1,line=2.25)

p1 <- density( Expellency[1,], na.rm = TRUE);   p1$y <- p1$y*diff(p1$x)[1]
p2 <- density( Expellency[2,], na.rm = TRUE);  p2$y <- p2$y*diff(p2$x)[1]  
p3 <- density( Expellency[3,], na.rm = TRUE);  p3$y <- p3$y*diff(p3$x)[1]
plot(p1,col='white',main='',xlim=c(0,2),ylim=c(0,1.2e-2),xaxs='i',yaxs='i',las=1)
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
mtext( bquote(bold("B")), side = 3, line = 0.1 ,at =0.08,cex=1) 
mtext(text=expression('Expellency (' ~phi* ' )'),side=1,line=2.25)
legend('topright',c('control','low','high'), fill=c(alpha('black',transparancy),alpha('darkorange',transparancy),alpha('deeppink',transparancy)),box.lty=0)


mtext(side = 2, 'Density', outer = TRUE, line = 1)
# 
dev.off()
setwd(cd)

# outcomes described in MS--------------------------------------------------------------

# low
median(Repellency[2,],na.rm=TRUE)
HDIofMCMC(mcmc(Repellency[2,]),credMass = 0.95)
#high
median(Repellency[3,],na.rm=TRUE)
HDIofMCMC(mcmc(Repellency[3,]),credMass = 0.95)

# low
1 - median(Expellency[2,])
1 - HDIofMCMC(mcmc(Expellency[2,]),credMass = 0.95)
#high
1 - median(Expellency[3,])
1 - HDIofMCMC(mcmc(Expellency[3,]),credMass = 0.95)
