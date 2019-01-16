# Code to run parametric survival analysis on Aedes Longevity data for different dosages
# Performs model selection, Proportional Hazards Testing and creates some summary figures

rm(list=ls())

# Load Packages -----------------------------------------------------------
library(Hmisc)
library(survival) # Note: requires 2.41-3! Figure scripts don't work with later versions
library(parfm)
library(flexsurv)
library(mgcv)
library(hBayesDM)

# Set work directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location (if working in RStudio)
getwd()

# Load source files --------------------------------------------------------
source('LongevityFunctions.R')

# load and prepare data ---------------------------------------------------
source('Longevity_PrepareData.R')

# Define Parameters -------------------------------------------------------
ng = 6 # number of groups
condDays = 15 # max number of days to examine in conditional survival analysis
con.time = 1 # days to condition on in conditional survival analysis
SR = c(0, rep(1,5))
num.draws = 5000 # to derive credible intervals
Dosages = c(0, 0.5, 0.75, 1, 1.25, 1.5)
colorlist= colorRampPalette(c("lightcyan3", "darkorange","deeppink"))(length(Dosages))
colorlist = c(rgb(0,0,0), colorlist)

# Conditional Survival Analysis -------------------------------------------
I <- which(SurvSR$TTE > con.time)
km.mod.con = survfit(Surv(TTE-con.time, Event == 'died')~ Dose, data=SurvSR[I,])
exp.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ Dose, data = SurvSR[I,], dist = "exponential")
weib.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ Dose, data = SurvSR[I,], dist = "weibull")
weib2.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ Dose + strata(Dose), data = SurvSR[I,], dist = "weibull")
lognorm.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ Dose, data = SurvSR[I,], dist = "lognormal")
gamma.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ Dose, data = SurvSR[I,], dist = "gamma")
gengamma.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ Dose, data = SurvSR[I,], dist = "gengamma")
gengamma2.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ Dose +sigma(Dose), data = SurvSR[I,], dist = "gengamma")

gengamma.orig.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ Dose, data = SurvSR[I,], dist = "gengamma.orig")

# save('km.mod.con','exp.mod.con', 'SurvSR', file = 'Lethality.Residual.Exp.RData')

# Model comparison --------------------------------------------------------
AICs.con=sapply(list(exp.mod.con, weib.mod.con, weib2.mod.con,lognorm.mod.con), extractAIC)
AICs.con = AICs.con[2,]
colnames(AICs.con, do.NULL =TRUE)
AICs.con = c(AICs.con, gamma.mod.con$AIC, gengamma.mod.con$AIC, gengamma2.mod.con$AIC)
names(AICs.con) = c("exp", "weib", "weib2", "lognorm", "gamma", "gengamma", "gengamma2")
AICs.con = rbind(AICs.con, AICs.con - AICs.con[6])

# Outcomes for paper ------------------------------------------------------
# fitted parameters for table:
temp = summary(exp.mod.con); Exp.pars = temp$table[,1]; 
Exp.pars = rbind( Exp.pars, temp$table[,1]-temp$table[,2]*1.96) 
Exp.pars = rbind( Exp.pars,temp$table[,1]+temp$table[,2]*1.96) 
summary(exp.mod.con)

temp = summary(weib.mod.con); Weib.pars = temp$table[,1]; 
Weib.pars = rbind( Weib.pars, temp$table[,1]-temp$table[,2]*1.96) 
Weib.pars = rbind( Weib.pars,temp$table[,1]+temp$table[,2]*1.96) 
Weib.pars[,3] =  1/(exp(Weib.pars[,3])); colnames(Weib.pars)[3] = "Shape"  

temp = summary(lognorm.mod.con); Logn.pars = temp$table[,1]; 
Logn.pars = rbind( Logn.pars, temp$table[,1]-temp$table[,2]*1.96) 
Logn.pars = rbind( Logn.pars,temp$table[,1]+temp$table[,2]*1.96) 
Logn.pars[,3] =  (exp(Logn.pars[,3])); colnames(Logn.pars)[3] = "Scale"     # 

temp$table = gamma.mod.con$res; Gamma.pars = temp$table[,1:3]; 

temp$table = gengamma.mod.con$res; Gengamma.pars = temp$table[,1:3]; 
temp$table = gengamma.orig.mod.con$res; Gengamma.orig.pars = temp$table[,1:3]; # same but with different parameterization

# AICs and delta AICs
AICs.con

# fractions of waiting times by dosage
Exp.coefs <- rmvn(n= num.draws, mu = exp.mod.con$coefficients, V = exp.mod.con$var)
AFT.g.exp <- exp(Exp.coefs[,2]) 
median(AFT.g.exp)
HDIofMCMC(AFT.g.exp, credMass = 0.95)
# high dosage
AFT.g.exp.h <- exp(Exp.coefs[,2]*1.5) 
median(AFT.g.exp.h)
HDIofMCMC(AFT.g.exp.h, credMass = 0.95)


Phis = rbind(exp(Exp.pars[,2]),
      exp(Weib.pars[,2]),
      exp(Logn.pars[,2]),
      1/exp(Gamma.pars[3,]),
      exp(Gengamma.pars[4,]))
rownames(Phis) = c('exp', 'weib', 'logn', 'gamma', 'gengamma')
colnames(Phis) = c('mean', '95_low', '95_high')


# Test PH and Weibull -----------------------------------------------------
Treat = numeric()
for (ii in 1:length(Dosages)){
  Treat = c(Treat,rep(Dosages[ii],km.mod.con$strata[[ii]]))
}

setwd(figure.path)

pdf('FigS2.pdf',width=4.75,height=4.75)
par(mfrow = c(1,1),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1)) 

# create log-log plot to explore PH and Weibull. 
plot(log(km.mod.con$time),log(-log(km.mod.con$surv)), type='n', xlab='', ylab='', las = 1)
for (i in 1:length(Dosages)){
  points(log(km.mod.con$time)[Treat== Dosages[i]], log(-log(km.mod.con$surv[Treat== Dosages[i]])), col=colorlist[i], lty=6, pch = 16)
  lines(lowess(log(-log(km.mod.con$surv[Treat== Dosages[i] ]))~ log(km.mod.con$time)[Treat== Dosages[i] ]), lty=1, col=colorlist[i])
}
legend('bottomright', c("control","SR Dose 0.5" ,"SR Dose 0.75", "SR Dose 1.00", "SR Dose 1.25", "SR Dose 1.5"), col=colorlist, lty=1, cex=.9, bty='n')  
  
mtext(side = 1, text = expression(ln(t)), line = 2.25, cex = 1 )
mtext(side = 2, text = expression(ln(-ln(S(t)))), line = 2.25, cex = 1  )

dev.off()
setwd(cd)

# Plot survival curves ----------------------------------------------------
setwd(figure.path)

# pdf('SurvivalCurves.pdf',width=7.50,height=2.70)
par(mfrow = c(1,5),mai=c(0.3,0,0.1,0.1) ,oma=c(3,4,3,1)) 

plot(km.mod.con, xlab='', ylab='',lty=2,col=colorlist)
for( i in 1:length(Dosages)){
  curve(pexp(x, 1/(exp(coef(exp.mod.con)[1] + coef(exp.mod.con)[2] * Dosages[i])), lower.tail=FALSE), 
        from=0, to=max(SurvSR$TTE), col=colorlist[i], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
}
mtext( bquote(bold("A")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Exponential', line = 1.1, cex = .9, font = 4)

par(yaxt='n')
plot(km.mod.con, xlab='', ylab='',lty=2,col=colorlist)
for( i in 1:length(Dosages)){
curve(pweibull(x, scale=exp(coef(weib.mod.con)[1]+ coef(weib.mod.con)[2]*Dosages[i]), shape=1/weib.mod.con$scale, lower.tail=FALSE), 
      from=0, to=max(SurvSR$TTE), col=colorlist[i], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
}
mtext( bquote(bold("B")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Weibull', line = 1.1, cex = .9, font = 4)

plot(km.mod.con, xlab='', ylab='',lty=2,col=colorlist)
for( i in 1:length(Dosages)){
  curve(plnorm(x, meanlog = (coef(lognorm.mod.con)[1]+ coef(lognorm.mod.con)[2]*Dosages[i]) , sdlog = (Logn.pars[1,3]), lower.tail=FALSE), 
        from=0, to=max(SurvSR$TTE), col=colorlist[i], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
}
mtext( bquote(bold("C")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Lognormal', line = 1.1, cex = .9, font = 4)

Dose<-Dosages

plot(km.mod.con, lty=2, col=colorlist) # xlab='time (days)', ylab=' Survival',
lines.flexsurvreg(gamma.mod.con,type="survival",newdata=data.frame(Dose),col = colorlist, col.ci=NULL, lwd=2)
mtext( bquote(bold("D")), side = 3, line = 0.1 ,at =4,cex=1)
mtext(side = 3, text = 'Gamma', line = 1.1, cex = .9, font = 4)

plot(km.mod.con, lty=2, col=colorlist) # xlab='time (days)', ylab=' Survival',
lines.flexsurvreg(gengamma.mod.con,type="survival",newdata=data.frame(Dose),col = colorlist, col.ci=NULL, lwd=2)
mtext( bquote(bold("E")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Generalized gamma', line = 1.1, cex = .9, font = 4)

legend('topright', c("0","0.5" ,"0.75", "1.00", "1.25", "1.5", "K-M"), 
       col=colorlist[c(seq(6),1)], lty=c(rep(1,6),2), cex=.9, bty='n')  

mtext('Time (days)',1,line=1, outer = TRUE)
mtext('Survival',2,line=2.25, outer = TRUE)

# dev.off()
setwd(cd)

# Dose Response Curve -----------------------------------------------------
# mark that this result converts the AFT coef to the PH coef
HR=exp(-1*weib.mod.con$coefficients[2]*Dosages * 1/weib.mod.con$scale)
# Hazard = exp(coef(weib.mod.con)[1]  + coef(weib.mod.con)[2]*Dosages)
barplot(HR,xlab = 'Dosages', ylab = 'Hazard Ratio',names.arg = Dosages)

# Hazard Ratios -----------------------------------------------------------

HR_exp=exp(-1*exp.mod.con$coefficients[2] * 1/exp.mod.con$scale)
HR_weib=exp(-1*weib.mod.con$coefficients[2] * 1/weib.mod.con$scale)
HR_loglog=exp(-1*loglog.mod.con$coefficients[2])
DelayFactor_GenGamma <- exp(gengamma.mod.con$coefficients[4]) 


# Dose Response Curve Delay Factor ----------------------------------------
Dosage <- c(seq(0,2,by=0.25))

cd <- getwd()
setwd("Figures")

pdf('SALongevityAedes.DoseDelay.gengam.pdf',width=4.5,height=4.5)
par(mar=rep(0,4),oma=c(4,4,.5,.5))

plot(Dosages,exp(gengamma.mod.con$coefficients[4]*Dosages), ylim=c(0.5,1.5), type="l", ylab = " ", lwd =4) #ylab="Biting Rate", 

mtext('Dosage',1,line=2.5)
mtext('Relative Survival Time',2,line=2.75)

dev.off()
setwd(cd)

# Save Models -------------------------------------------------------------
cd <- getwd()
setwd('Output')

save('weib.mod.con', 'weib2.mod.con', 'gengamma.mod.con','gengamma2.mod.con','gamma.mod.con','loglog.mod.con', 'exp.mod.con' , file='SALongevityAedes.Models.Con.RData')

setwd(cd)


# Plot Gen Gamma for Epidemics --------------------------------------------
Dosage <- c(seq(0,2,by=0.25))
colorlist= colorRampPalette(c("lightcyan3", "darkorange","deeppink"))(length(Dosage))
colorlist = c(rgb(0,0,0), colorlist)

# GenGamma
cd <- getwd()
setwd("Figures")

pdf('SALongevityAedes_GenGamma_Epidemics.pdf',width=4.5,height=4.5)
par(mar=rep(0,4),oma=c(4,4,.5,.5))
Dose<-Dosages
plot(km.mod.con, lty=2, col='white') # xlab='time (days)', ylab=' Survival',
lines.flexsurvreg(gengamma.mod.con,type="survival",newdata=data.frame(Dose),col = colorlist, col.ci=NULL, lwd=4)
#legend('topright', c("control","SR Dose 0.5" ,"SR Dose 0.75", "SR Dose 1.00", "SR Dose 1.25", "SR Dose 1.5", "Kaplan Meier"), col=colorlist[c(seq(6),1)], lty=c(rep(1,6),2), cex=.9, bty='n')  

mtext('Time (days)',1,line=2.1)
mtext('Survival',2,line=2.1)

dev.off()
setwd(cd)

