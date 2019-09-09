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

# Define Parameters -------------------------------------------------------
ng = 6 # number of groups
condDays = 15 # max number of days to examine in conditional survival analysis
con.time = 1 # days to condition on in conditional survival analysis
SR = c(0, rep(1,5))
num.draws = 5000 # to derive credible intervals
Dosages = c(0, 0.5, 0.75, 1, 1.25, 1.5)
colorlist= colorRampPalette(c("lightcyan3", "darkorange","deeppink"))(length(Dosages))
colorlist = c(rgb(0,0,0), colorlist)

# load and prepare data ---------------------------------------------------
source('Longevity_PrepareData.R')

# Conditional Survival Analysis (factor)-------------------------------------------
km.mod.con = survfit(Surv(TTE-con.time, Event == 'died')~ factor(Dose), data=SurvSR.con)
exp.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose), data = SurvSR.con, dist = "exponential")
weib.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose), data = SurvSR.con, dist = "weibull")
weib2.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose) + strata(Dose), data = SurvSR.con, dist = "weibull")
lognorm.mod.con <- survreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose), data = SurvSR.con, dist = "lognormal")
gamma.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose), data = SurvSR.con, dist = "gamma")
gengamma.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose), data = SurvSR.con, dist = "gengamma")
gengamma2.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose) +sigma(Dose), data = SurvSR.con, dist = "gengamma")

gengamma.orig.mod.con <- flexsurvreg(Surv(TTE-con.time, Event == 'died') ~ factor(Dose), data = SurvSR.con, dist = "gengamma.orig")

save('km.mod.con','exp.mod.con', 'SurvSR.con', file = "../Input/Lethality.Residual.Exp.RData")


# Model comparison --------------------------------------------------------
AICs.con=sapply(list(exp.mod.con, weib.mod.con, weib2.mod.con,lognorm.mod.con), extractAIC)
AICs.con = AICs.con[2,]
colnames(AICs.con, do.NULL =TRUE)
AICs.con = c(AICs.con, gamma.mod.con$AIC, gengamma.mod.con$AIC) #, gengamma2.mod.con$AIC
names(AICs.con) = c("exp", "weib", "weib2", "lognorm", "gamma", "gengamma") #, "gengamma2"
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
Weib.pars[,7] =  1/(exp(Weib.pars[,7])); colnames(Weib.pars)[7] = "Shape"  

temp = summary(lognorm.mod.con); Logn.pars = temp$table[,1]; 
Logn.pars = rbind( Logn.pars, temp$table[,1]-temp$table[,2]*1.96) 
Logn.pars = rbind( Logn.pars,temp$table[,1]+temp$table[,2]*1.96) 
Logn.pars[,7] =  (exp(Logn.pars[,7])); colnames(Logn.pars)[7] = "Scale"     # 

temp$table = gamma.mod.con$res; Gamma.pars = temp$table[,1:3]; 

temp$table = gengamma.mod.con$res; Gengamma.pars = temp$table[,1:3]; 
temp$table = gengamma.orig.mod.con$res; Gengamma.orig.pars = temp$table[,1:3]; # same but with different parameterization

# AICs and delta AICs
AICs.con

# fractions of waiting times by dosage
Exp.coefs <- rmvn(n= num.draws, mu = exp.mod.con$coefficients, V = exp.mod.con$var)
AFT.g.exp <- exp(Exp.coefs[,4]) 
median(AFT.g.exp)
HDIofMCMC(AFT.g.exp, credMass = 0.95)
# high dosage
AFT.g.exp.h <- exp(Exp.coefs[,6]) 
median(AFT.g.exp.h)
HDIofMCMC(AFT.g.exp.h, credMass = 0.95)

Phis = rbind(exp(Exp.pars[,4]),
             exp(Weib.pars[,4]),
             exp(Logn.pars[,4]),
             1/exp(Gamma.pars[5,]),
             exp(Gengamma.pars[6,]))
rownames(Phis) = c('exp', 'weib', 'logn', 'gamma', 'gengamma')
colnames(Phis) = c('mean', '95_low', '95_high')


# Test PH and Weibull -----------------------------------------------------
Treat = numeric()
for (ii in 1:length(Dosages)){
  Treat = c(Treat,rep(Dosages[ii],km.mod.con$strata[[ii]]))
}

pdf('FigS2.pdf',width=4.75,height=4.75)
par(mfrow = c(1,1),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1)) 

# create log-log plot to explore PH and Weibull. 
plot(log(km.mod.con$time),log(-log(km.mod.con$surv)), type='n', xlab='', ylab='', las = 1)
for (i in c(1,4,6)){
  points(log(km.mod.con$time)[Treat== Dosages[i]], log(-log(km.mod.con$surv[Treat== Dosages[i]])), col=colorlist[i], lty=6, pch = 16)
  lines(lowess(log(-log((km.mod.con$surv+0.00000001)[Treat== Dosages[i] ]))~ log(km.mod.con$time)[Treat== Dosages[i] ]), lty=1, col=colorlist[i])
}
legend('bottomright', c("control","low", "high"), col=colorlist[c(1,4,6)], lty=1, cex=.9, bty='n')  
  
mtext(side = 1, text = expression(ln(t)), line = 2.25, cex = 1 )
mtext(side = 2, text = expression(ln(-ln(S(t)))), line = 2.25, cex = 1  )

dev.off()

# Plot survival curves (factor)----------------------------------------------------
setwd(figure.path)

pdf('SurvivalCurves.pdf',width=7.50,height=2.70)
par(mfrow = c(1,5),mai=c(0.3,0,0.1,0.1) ,oma=c(3,4,3,1)) 

plot(km.mod.con[c(1,4,6)], xlab='', ylab='',lty=2,col=colorlist[c(1,4,6)])
curve(pexp(x, 1/(exp(coef(exp.mod.con)[1] )), lower.tail=FALSE), 
      from=0, to=max(SurvSR.con$TTE), col=colorlist[1], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
for( i in c(4,6)){
  curve(pexp(x, 1/(exp(coef(exp.mod.con)[1] + coef(exp.mod.con)[i] )), lower.tail=FALSE), 
        from=0, to=max(SurvSR.con$TTE), col=colorlist[i], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
}
mtext( bquote(bold("A")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Exponential', line = 1.1, cex = .9, font = 4)

par(yaxt='n')
plot(km.mod.con[c(1,4,6)], xlab='', ylab='',lty=2,col=colorlist[c(1,4,6)])
curve(pweibull(x, scale=exp(coef(weib.mod.con)[1]), shape=1/weib.mod.con$scale, lower.tail=FALSE), 
      from=0, to=max(SurvSR.con$TTE), col=colorlist[1], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
for( i in c(4,6)){
  curve(pweibull(x, scale=exp(coef(weib.mod.con)[1]+ coef(weib.mod.con)[i]), shape=1/weib.mod.con$scale, lower.tail=FALSE), 
        from=0, to=max(SurvSR.con$TTE), col=colorlist[i], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
}
mtext( bquote(bold("B")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Weibull', line = 1.1, cex = .9, font = 4)

plot(km.mod.con[c(1,4,6)], xlab='', ylab='',lty=2,col=colorlist[c(1,4,6)])
curve(plnorm(x, meanlog = (coef(lognorm.mod.con)[1]) , sdlog = (Logn.pars[1,3]), lower.tail=FALSE), 
      from=0, to=max(SurvSR.con$TTE), col=colorlist[1], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
for( i in c(4,6)){
  curve(plnorm(x, meanlog = (coef(lognorm.mod.con)[1]+ coef(lognorm.mod.con)[i]) , sdlog = (Logn.pars[1,3]), lower.tail=FALSE), 
        from=0, to=max(SurvSR.con$TTE), col=colorlist[i], ylab =expression(hat(S)(t)), xlab='t',add=T, lwd=2)
}
mtext( bquote(bold("C")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Lognormal', line = 1.1, cex = .9, font = 4)

Dose<-Dosages

plot(km.mod.con[c(1,4,6)], xlab='', ylab='',lty=2,col=colorlist[c(1,4,6)])
lines(gamma.mod.con, col = c(colorlist[1],rep('white',2),colorlist[4],'white',colorlist[6]), col.ci=NULL, lwd=2)
lines(km.mod.con[c(1,4,6)], xlab='', ylab='',lty=2,col=colorlist[c(1,4,6)])
mtext( bquote(bold("D")), side = 3, line = 0.1 ,at =4,cex=1)
mtext(side = 3, text = 'Gamma', line = 1.1, cex = .9, font = 4)

plot(km.mod.con[c(1,4,6)], lty=2, col=colorlist[c(1,4,6)]) # xlab='time (days)', ylab=' Survival',
lines(gengamma.mod.con, col = c(colorlist[1],rep('white',2),colorlist[4],'white',colorlist[6]), col.ci=NULL, lwd=2)
lines(km.mod.con[c(1,4,6)], xlab='', ylab='',lty=2,col=colorlist[c(1,4,6)])
mtext( bquote(bold("E")), side = 3, line = 0.1 ,at =4,cex=1) 
mtext(side = 3, text = 'Generalized gamma', line = 1.1, cex = .9, font = 4)

legend('topright',c("control","low", "high"), col=colorlist[c(1,4,6)], lty=c(rep(1,6),2), cex=.9, bty='n')  

mtext('Time (days)',1,line=1, outer = TRUE)
mtext('Survival',2,line=2.25, outer = TRUE)

dev.off()
setwd(cd)

# Hazard Ratios -----------------------------------------------------------

HR_exp=exp(-1*exp.mod.con$coefficients[2] * 1/exp.mod.con$scale)
HR_weib=exp(-1*weib.mod.con$coefficients[2] * 1/weib.mod.con$scale)
HR_loglog=exp(-1*loglog.mod.con$coefficients[2])
DelayFactor_GenGamma <- exp(gengamma.mod.con$coefficients[4]) 

# Save Models -------------------------------------------------------------
cd <- getwd()
setwd('Output')

save('weib.mod.con', 'weib2.mod.con', 'gengamma.mod.con','gengamma2.mod.con','gamma.mod.con','loglog.mod.con', 'exp.mod.con' , file='SALongevityAedes.Models.Con.RData')

setwd(cd)

