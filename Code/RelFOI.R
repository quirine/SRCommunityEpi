rm(list=ls())

# load libraries ----------------------------------------------------------
library(mgcv)
library(hBayesDM)
library(coda)
library(ggplot2)
library(Hmisc)

# Set work directory -------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location (if working in RStudio)
getwd()

# Plotting parameters -----------------------------------------------------
Coverages=seq(0,1,.01)
Dosages = seq(0,2,0.001)
colorlist= colorRampPalette(c("lightcyan3", "darkorange","deeppink"))(length(Dosages))
colorlist = c(rgb(0,0,0), colorlist)
num.draws = 1000
Names = c( 'biting rate', 'gonotrophic cycle', 'acute mortality', 'delayed mortality', 'repellency', 'expellency')

# Load models and source files-------------------------------------------------------------
load("../Input/BloodFeeding.RData")
load("../Input/Lethality.Instantaneous.RData")
load("../Input/Lethality.Residual.Exp.RData")
load("../Input/Movement.RData")
source('RelFOI.Functions.R')
source('RelFOI.Parameters.R')

# Outcomes for paper ------------------------------------------------------
median(Rho)
HDIofMCMC(mcmc(Rho),credMass = 0.95)
median(Rho.h)
HDIofMCMC(mcmc(Rho.h),credMass = 0.95)
median(Q.t)
HDIofMCMC(mcmc(Q.t),credMass = 0.95)
median(Q.t.h)
HDIofMCMC(mcmc(Q.t.h),credMass = 0.95)

# reduction of 50% FoI at a coverage of:
temp = numeric()
for( ii in 1:num.draws){     
  temp1 = 1 - rel_foi(Coverages,Mu[ii],Rho[ii],Alpha[ii],tau,Q.t[ii]*q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii])
  temp = c( temp, Coverages[ which.min(abs(temp1-0.5)) ])
}
foi.reduction = 1 - rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc)
Coverages[ which.min(abs(foi.reduction-0.5)) ]
HDIofMCMC(mcmc(temp), credMass = 0.95)
# and for high dosage
temp = numeric()
for( ii in 1:num.draws){     
  temp1 = 1 - rel_foi(Coverages,Mu.h[ii],Rho.h[ii],Alpha.h[ii],tau,Q.t.h[ii]*q.u,q.u,G.T.h[ii],g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii])
  temp = c( temp, Coverages[ which.min(abs(temp1-0.5)) ])
}
foi.reduction = 1 - rel_foi(Coverages,mu.h,rho.h,alpha.h,tau,q.u*q.t.scaler.h,q.u,g.t.h,g.u,g.u*g.tau.scaler,gc,alpha.gc.h)
Coverages[ which.min(abs(foi.reduction-0.5)) ]
HDIofMCMC(mcmc(temp), credMass = 0.95)

# max effect of product
foi.reduction.rep.bit = 1 - rel_foi(Coverages,mu=mu,rho,alpha,tau,q.t.scaler*q.u,q.u,g.t=g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=Mu[ii],Rho[ii],Alpha[ii],tau,Q.t[ii]*q.u,q.u,g.t=G.T[ii],g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit[length(Coverages)]

# max effect of product without lethality
foi.reduction.rep.bit = 1 - rel_foi(Coverages,mu=0,rho,alpha,tau,q.t.scaler*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit[length(Coverages)]

# and for high dosage
foi.reduction.rep.bit = 1 - rel_foi(Coverages,mu=0,rho.h,alpha.h,tau,q.t.scaler.h*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc.h)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho.h[ii],Alpha.h[ii],tau,Q.t.h[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit[length(Coverages)]

# max effect of product without lethality and expellency
foi.reduction.rep.bit = 1 - rel_foi(Coverages,mu=0,rho,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit[length(Coverages)]

# and for high dosage
foi.reduction.rep.bit.h = 1 - rel_foi(Coverages,mu=0,rho.h,alpha.h,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc.h)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho.h[ii],Alpha.h[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit.h[length(Coverages)]

# without probing effect
foi.reduction.rep.bit = 1 - rel_foi(Coverages,mu=0,rho,alpha=alpha.gc,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha.gc[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit[length(Coverages)]

# a high dosage
foi.reduction.rep.bit.h = 1 - rel_foi(Coverages,mu=0,rho.h,alpha=alpha.gc.h,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc.h)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho.h[ii],Alpha.gc.h[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit.h[length(Coverages)]

# max effect of product without lethality with positive expellency
q.t.scaler.positive = 1.3
foi.reduction.rep.bit.exp.pos = 1 - rel_foi(Coverages,mu=0,rho,alpha,tau,q.t=q.u*q.t.scaler.positive,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,q.t.scaler.positive*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit.exp.pos[length(Coverages)]

# max effect of product without lethality with expellency with higher outdoor mortality
g.tau.scaler.positive = 3
foi.reduction.rep.bit.exp.pos.outdoor = 1 - rel_foi(Coverages,mu=0,rho,alpha,tau,q.t=q.u*q.t.scaler.positive,q.u,g.t=g.u,g.u,g.u*g.tau.scaler.positive,gc,alpha.gc)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,q.t.scaler.positive*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler.positive,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,length(Coverages)]),]
HDIofMCMC(mcmc(1 - temp[,length(Coverages)]), credMass = 0.95 ); 
foi.reduction.rep.bit.exp.pos.outdoor[length(Coverages)]

# Fig 3 ------------------------------------------------
transparency = 0.2
YLIMS = c(0,1.2)
line.weight = 2
index.l = round(0.025 * num.draws); 
index.h = round(0.975 * num.draws); 
cov.index = which(Coverages == .5)

pdf('Fig3.pdf',width=7.5,height=3)
old.par <- par(mfrow=c(1,4),oma = c(5,4,4,2),
               mar = c(0,1.4,1,0))

plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu[ii],Rho[ii],Alpha[ii],tau,Q.t[ii]*q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')

temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu.h[ii],Rho.h[ii],Alpha.h[ii],tau,Q.t.h[ii]*q.u,q.u,G.T.h[ii],g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')

lines(Coverages,rel_foi(Coverages,mu.h,rho.h,alpha.h,tau,q.u*q.t.scaler.h,q.u,g.t.h,g.u,g.u*g.tau.scaler,gc,alpha.gc.h),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')

# mtext('all effects',3,line=1.4,cex=0.8)
mtext( bquote(bold("A")), side = 3, line = 0.1 ,at = 0.02,cex=1) 
mtext( bquote(bold("All effects")), side = 3, line = 1.1 ,cex=.9) 

plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')

temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho.h[ii],Alpha.h[ii],tau,Q.t.h[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha,tau,q.t.scaler*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
lines(Coverages,rel_foi(Coverages,mu=0,rho.h,alpha.h,tau,q.t.scaler.h*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc.h),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')

# mtext('all but mortality effects',3,line=1.4,cex=0.8)
mtext( bquote(bold("B")), side = 3, line = 0.1 ,at = 0.02,cex=1) 
mtext( bquote(bold("As A minus lethality")), side = 3, line = 1.1 ,cex=.9) 

plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')

temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho.h[ii],Alpha.h[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
lines(Coverages,rel_foi(Coverages,mu=0,rho.h,alpha.h,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc.h),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')

# mtext(paste('all but mortality \n and expellency effects'),3,line=1.4,cex=0.8)
mtext( bquote(bold("C")), side = 3, line = 0.1 ,at = 0.02,cex=1) 
mtext( bquote(bold("As B minus expellency")), side = 3, line = 1.1 ,cex=.9) 

plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho[ii],Alpha.gc[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')

temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,Rho.h[ii],Alpha.gc.h[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha.gc,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
lines(Coverages,rel_foi(Coverages,mu=0,rho.h,alpha.gc.h,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc.h),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
# mtext(paste('all but mortality, expellency, \n and probing effects'),3,line=1.4,cex=0.8)
mtext( bquote(bold("D")), side = 3, line = 0.1 ,at = 0.02,cex=1) 
mtext( bquote(bold("As C minus probing")), side = 3, line = 1.1 ,cex=.9) 

legend('topright',c('low dosage', 'high dosage'), 
       fill = c(alpha('darkorange',transparency),alpha('deeppink',transparency)), bty = 'n')

mtext('Coverage',1,line=2.25, outer = TRUE)
mtext('Mean FoI relative to baseline',2,line=1, outer = TRUE)

dev.off()


# Fig 4----------------------------------
transparency = 0.5
# Construct vectors that capture uncertainty by effect
YLIMS = c(0,1.5)
title.size = 0.7
letter.size = 0.8
line.weight = 2
index.l = round(0.025 * num.draws); 
index.h = round(0.975 * num.draws);
cov.index = which(Coverages == 0.5)
jj = 0 # counter for row letters

pdf('Fig4.pdf',width=7.5,height=7.5)
old.par <- par(mfrow=c(6,6),
               oma = c(6,5,4,2) + 0.1,
               mar = c(0,1.3,1.3,0) + 0.1) 
# biting rate
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',pty='s',xaxs='i',yaxs='i',las=1, bty = 'n')
box(lwd=3) 
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')

temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,Alpha.h[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha.h,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')

mtext(Names[1], side = 3, line = 1, cex = title.size, font = 4)
mtext(Names[1], side = 2, line = 2.25, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1

# biting and oviposition high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='',ylim=YLIMS)
mtext(Names[2], side = 3, line = 1, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# biting and instantaneous death high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu[ii],rho=0,Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext(Names[3], side = 3, line = 1, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# biting and residual death high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,Alpha[ii],tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext(Names[4], side = 3, line = 1, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# biting and repellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext(Names[5], side = 3, line = 1, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# biting and expellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,Alpha[ii],tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext(Names[6], side = 3, line = 1, cex = title.size, font = 4)  
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1

# biting and oviposition
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='',ylim=YLIMS)
mtext(Names[2], side = 2, line = 2.25, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1, bty = 'n')
box(lwd=3) 
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')

temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc.h[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')

lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc.h),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1

# oviposition and instantaneous death high dosage 
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu[ii],rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition and residual death high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition and repellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition and expellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1

# biting and instantaneous death
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu[ii],rho=0,Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext(Names[3], side = 2, line = 2.25, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition and instantaneous death 
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu[ii],rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# instantaneous death 
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1, bty = 'n')
box(lwd=3) 
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu[ii],rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,Mu.h[ii],rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
lines(Coverages,rel_foi(Coverages,mu.h,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1

# instantaneous death and residual death high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,Mu[ii],rho=0,alpha=1,tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# instantaneous death and repellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,Mu[ii],Rho[ii],alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu,rho,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# instantaneous death and expellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,Mu[ii],rho=0,alpha=1,tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha=1,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# biting and residual death
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,Alpha[ii],tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext(Names[4], side = 2, line = 2.25, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition and residual death
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# instantaneous death and residual death
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,Mu[ii],rho=0,alpha=1,tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# residual death
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1, bty = 'n')
box(lwd=3) 
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,G.T.h[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.t=q.u,q.u,g.t.h,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# residual death and repellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],alpha=1,tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha=1,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# residual death and expellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,Q.t[ii]*q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# biting and repellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],Alpha[ii],tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext(Names[5], side = 2, line = 2.25, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition and repellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# instantaneous death and repellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,Mu[ii],Rho[ii],alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu,rho,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# residual death and repellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],alpha=1,tau,q.t=q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha=1,tau,q.t=q.u,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# repellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS, xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1, bty = 'n')
box(lwd=3) 
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho.h[ii],alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')
lines(Coverages,rel_foi(Coverages,mu=0,rho.h,alpha=1,tau,q.t=q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# repellency and expellency high dosage
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,xaxt='n',yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],alpha=1,tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 0)
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha=1,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1

# biting and expellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,Alpha[ii],tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext(Names[6], side = 2, line = 2.25, cex = title.size, font = 4)
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# oviposition and expellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,Alpha.gc[ii]))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# instantaneous death and expellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,Mu[ii],rho=0,alpha=1,tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu,rho=0,alpha=1,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# residual death and expellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,Q.t[ii]*q.u,q.u,G.T[ii],g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
     type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# repellency and expellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1)
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,Rho[ii],alpha=1,tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
lines(Coverages,rel_foi(Coverages,mu=0,rho,alpha=1,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
# expellency
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=YLIMS,yaxt ='n',pty='s',xaxs='i',yaxs='i',las=1, bty = 'n')
box(lwd=3) 
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,Q.t[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('darkorange',transparency), border = 'darkorange')
temp = numeric()
for( ii in 1:num.draws){     
  temp = rbind(temp, rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,Q.t.h[ii]*q.u,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1))
}
temp = temp[order(temp[,cov.index]),]
polygon(c(Coverages, rev(Coverages)),c(temp[index.l,],rev(temp[index.h,]) ), col=alpha('deeppink',transparency), border = 'deeppink')

lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.u*q.t.scaler,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='darkorange',lwd=line.weight,xlab='',ylab='')
lines(Coverages,rel_foi(Coverages,mu=0,rho=0,alpha=1,tau,q.u*q.t.scaler.h,q.u,g.t=g.u,g.u,g.u*g.tau.scaler,gc,alpha.gc=1),
      type = 'l',col='deeppink',lwd=line.weight,xlab='',ylab='')
mtext((paste0(LETTERS[floor(jj/6)+1],letters[(jj%%6)+1])), side = 3, line = 0.1 ,at = 0.1,cex=letter.size, font = 4)
jj = jj + 1
legend('bottomright',c('low dosage', 'high dosage'), 
       fill = c(alpha('darkorange',transparency),alpha('deeppink',transparency)), bty = 'n')

mtext('Coverage',1,line=3, outer = TRUE)
mtext('Mean FoI relative to baseline',2,line=2.25, outer = TRUE)

dev.off()

# FigS4----------------------------------------------------
num.sa = 21
text.size = 0.8

colorlist.2= colorRampPalette(c("yellow","red","darkred"))(num.sa)

# Control parameters (based on Survival Analysis of Bloodfeeding and Longevity data)
Ns = seq(from = 3, to = 33, length.out = num.sa)   # eip (Range as per Chan 2012)
Bs = seq(from = 0, to = 1, length.out = num.sa)   
Cs = seq(from = 0, to = 1, length.out = num.sa)   
Xs = seq(from = 0.01, to = 0.99, length.out = num.sa) # human prevalence
Gcs = seq(from = 1, to = 14, length.out = num.sa)   # gonotrophic cycle length
Aus = seq(from = 1/5, to = 10, length.out = num.sa  )             # biting rate
G.scalers = seq(from = 0.1, to = 10, length.out = num.sa) 
G.tau.scalers = seq(from = g.tau.scaler*0.1, to = g.tau.scaler*10, length.out = num.sa) # mortality indoor vs outdoor
Qus = seq(from = q.u*10, to = q.u/10, length.out = num.sa )                # residence time
Tau.scalers = seq(from = 0.01, to = 10, length.out = num.sa)          # proportion of time spent outdoors

n.temp = n; b.temp = b; c.temp = c; X.temp = X; a.temp = a.u;

cd <- getwd()
setwd(Figure.path)

pdf('FigS4.pdf',width=10,height=4.5)
old.par <- par(mfrow=c(2,5),oma = c(5,4,2,1) + 0.1,
               mar = c(0,1,3.3,0) + 0.1)

n = Ns[1]
plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',xaxt='n', las = 1)
for( ii in 2:num.sa){
  n = Ns[ii]
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
n = n.temp
mtext( bquote(bold("A")), side = 3, line = 0.1 ,at = 0.05,cex=1) 
mtext( expression(bold('EIP')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

b = Bs[1]
plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',xaxt='n',yaxt='n')
for( ii in 2:num.sa){
  b = Bs[ii]
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
b = b.temp
mtext( bquote(bold("B")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Transmission \n prob m to h')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

c = Cs[1]
plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',xaxt='n',yaxt='n')
for( ii in 2:num.sa){
  c = Cs[ii]
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
c = c.temp
mtext( bquote(bold("C")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Transmission \n prob h to m')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

X = Xs[1]
plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',xaxt='n',yaxt='n')
for( ii in 2:num.sa){
  X = Xs[ii]
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
X = X.temp
mtext( bquote(bold("D")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Human infection \n prevalence')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,Gcs[1],alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',xaxt='n',yaxt='n')
for( ii in 2:num.sa){
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,Gcs[ii],alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
mtext( bquote(bold("E")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Gonotrophic cycle \n length')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

a.u = Aus[1]
plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='', las = 1)
for( ii in 2:num.sa){
  a.u = Aus[1]
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
a.u = a.temp
mtext( bquote(bold("F")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Biting rate')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.t=q.u,q.u,G.scalers[1]*(g.u/aft.g),G.scalers[1]*g.u,G.scalers[1]*g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',yaxt='n')
for( ii in 2:num.sa){
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.t=q.u,q.u,G.scalers[ii]*(g.u/aft.g),G.scalers[ii]*g.u,G.scalers[ii]*g.u*g.tau.scaler,gc,alpha.gc=1),col=colorlist.2[ii],type='l',ylim=c(0,3),lwd=2)
}
mtext( bquote(bold("G")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Mortality rate')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.u/aft.g,g.u,g.u*G.tau.scalers[1],gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',yaxt='n')
for( ii in 2:num.sa){
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,q.u*q.t.scaler,q.u,g.u/aft.g,g.u,g.u*G.tau.scalers[ii],gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
mtext( bquote(bold("H")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Transit vs \n indoor mortality')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,Qus[1]*q.t.scaler,Qus[1],g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',yaxt='n')
for( ii in 2:num.sa){
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,tau,Qus[ii]*q.t.scaler,Qus[ii],g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
mtext( bquote(bold("I")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Residence \n time')), side = 3, line = 0.1 ,cex=text.size, font = 4) 

plot(Coverages,rel_foi(Coverages,mu,rho,alpha,tau=(Tau.scalers[1]*(1/q.u)),q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[1],type='l',ylim=c(0,1),lwd=2, ylab ='', xlab='',yaxt='n')
for( ii in 2:num.sa){
  lines(Coverages,rel_foi(Coverages,mu,rho,alpha,Tau.scalers[ii]*(1/q.u),q.u*q.t.scaler,q.u,g.t,g.u,g.u*g.tau.scaler,gc,alpha.gc),col=colorlist.2[ii],type='l',ylim=c(0,1),lwd=2)
}
mtext( bquote(bold("J")), side = 3, line = 0.1 ,at = 0.05,cex=1)
mtext( expression(bold('Transit vs \n indoor time')), side = 3, line = 0.1 ,cex=text.size, font = 4) 
 
mtext('Coverage',1,line=2.25, outer = TRUE)
mtext('Relative FoI',2,line=2, outer = TRUE)

dev.off()
setwd(cd)
