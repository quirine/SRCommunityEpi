# Code to derive instantaneous death rate on Aedes Longevity data for different dosages
# Performs model selection, Proportional Hazards Testing and creates some summary figures

rm(list=ls())
# Load Packages -----------------------------------------------------------
library(emdbook)
library(plotrix) 
library(bbmle)
library(lattice)
library(MASS)
library(Hmisc,warn=FALSE)
library(car)
library(mgcv)

# Set work directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location (if working in RStudio)
getwd()

# Load source files --------------------------------------------------------
source('LongevityFunctions.R')

# Define Parameters -------------------------------------------------------
ng = 6 # number of groups
SR = c(0, rep(1,5))
con.time = 1
Initial = 125
Dosages = c(0, 1, 1.5)
num.samples = 1e4
colorlist= colorRampPalette(c("lightcyan3", "darkorange","deeppink"))(length(Dosages))
colorlist = c(rgb(0,0,0), colorlist)

# load and prepare data ---------------------------------------------------
source('Longevity_PrepareData.R')
parnames(binom.instant.kd) <- c("beta0", "beta1")

K <- numeric()
for (ii in Dosages){
  Ind = which(SurvSR.init$Dose==ii)  
  K = c(K,length(Ind))
}
 
# Fit rates ---------------------------------------------------------------

kd.dose <- mle2(minuslogl = binom.instant.kd, start = c( beta0 = -7.5, beta1 = -1), data = list(N = Initial, k = K, dose = Dosages))

save('kd.dose', file = '../Input/Lethality.Instantaneous.RData')

# Get confidence intervals of outcomes ------------------------------------

MV.Coefs = rmvn(n=1e4, mu = kd.dose@coef, V = kd.dose@vcov)

# Estimates for MS -----------------------------------------------------
Probs.by.dose <- exp(kd.dose@coef[1] + kd.dose@coef[2]*Dosages) 

Probs.by.dose[which(Dosages==1)]
Probs.by.dose.uncertainty <- exp(MV.Coefs[,1] + MV.Coefs[,2]) 
median(Probs.by.dose.uncertainty)
HDIofMCMC(mcmc(Probs.by.dose.uncertainty),credMass = 0.95)

Probs.by.dose[which(Dosages==1.5)]
Probs.by.dose.uncertainty.15 <- exp(MV.Coefs[,1] + MV.Coefs[,2] * 1.5) 
median(Probs.by.dose.uncertainty.15)
HDIofMCMC(mcmc(Probs.by.dose.uncertainty.15),credMass = 0.95)

temp = summary(kd.dose); Exp.pars = temp@coef[,1]; 
Exp.pars = rbind( Exp.pars, temp@coef[,1]-temp@coef[,2]*1.96) 
Exp.pars = rbind( Exp.pars,temp@coef[,1]+temp@coef[,2]*1.96) 
Probs.by.dose[which(Dosages==1.5)]
# for table
summary(kd.dose)
confint(kd.dose)


