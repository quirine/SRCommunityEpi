# Multinomial Bloodfeeding: estimating the probability of partial and full feeds concurrently
rm(list=ls())

# Loading Packages  -------------------------------------------------------
library(emdbook)
library(plotrix) 
library(bbmle)
library(lattice)
library(MASS)
library(Hmisc)
library(car)
library(lmtest)
library(binom)
library(mgcv)

# Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location (if working in RStudio)
getwd()


# Load source files -------------------------------------------------------
source('BloodfeedingFunctions.R')
source('LongevityFunctions.R')

# Data loading and processing ---------------------------------------------
BfData <- read.csv('../Data/SR_Bloodfeeding.csv')
Time.min <- BfData$Time * 60 + 10 # adding exposure time to time points
BfData <- data.frame(BfData, Time.min)
attach(BfData)

# Plotting parameters -----------------------------------------------------
num.draws = 10000
Dosage <- c(seq(0,2,by=0.25))
Dosage.small <- c(0,1,1.5)
colorlist= colorRampPalette(c("lightcyan3", "darkorange","deeppink"))(length(Dosage))
colorlist = c(rgb(0,0,0), colorlist)

# Multinomial Model -------------------------------------------------------

multinomNLL1 = function(params,k_partial,k_full,k_non,N,dose,t){ 
  beta0_part <- params[1]
  beta1_part <- params[2]
  beta0_full <- params[3]
  beta1_full <- params[4]
  
  rate_part = 1/(exp(beta0_part + beta1_part * dose ))  # mean number of partial bites per time step
  rate_full = 1/(exp(beta0_full + beta1_full * dose )) 
  
  prob_part = 1 - dpois(0,rate_part * t)  # probability of at least one partial bite over t 
  prob_full =  pexp(t, rate_full) # probability of a full blood meal after time = t
  
  prob_part_nofull = prob_part * (1 - prob_full) 
  prob_non = (1 - prob_full) * (1 - prob_part)  
  
  -sum(dmultinom(c(k_partial,k_full,k_non),prob=c(prob_part_nofull,
                                                  prob_full,
                                                  prob_non),
                 size=NULL,log=TRUE))
}

parnames(multinomNLL1) <- c("beta0_part", "beta1_part", "beta0_full", "beta1_full")

# Fit Multinomial Models ------------------------------------------------------
# without dose effect
bf.nodose <- mle2(minuslogl = multinomNLL1, start = c( beta0_part = 7, beta1_part = 0.15, beta0_full = 7, beta1_full = 0.15),
                  data = list(N = Initial, k_partial = Partial_Bloodfed, k_full = Full_Bloodfed, k_non = Non_Bloodfed, dose = rep(0,nrow(BfData)), t = Time.min))

# with dose effect
bf <- mle2(minuslogl = multinomNLL1, start = c( beta0_part = 7, beta1_part = 0.15, beta0_full = 7, beta1_full = 0.15),
           data = list(N = Initial, k_partial = Partial_Bloodfed, k_full = Full_Bloodfed, k_non = Non_Bloodfed, dose = Dose, t = Time.min))

# for each dose separately
ind = which(BfData$Dose == 0)
bf.control <- mle2(minuslogl = multinomNLL1, start = c( beta0_part = 7, beta1_part = 0.15, beta0_full = 7, beta1_full = 0.15),
                   data = list(N = Initial[ind], k_partial = Partial_Bloodfed[ind], k_full = Full_Bloodfed[ind], k_non = Non_Bloodfed[ind], dose = rep(0,length(ind)), t = Time.min[ind]))

ind = which(BfData$Dose == 1.0)
bf.low <- mle2(minuslogl = multinomNLL1, start = c( beta0_part = 7, beta1_part = 0.15, beta0_full = 7, beta1_full = 0.15),
               data = list(N = Initial[ind], k_partial = Partial_Bloodfed[ind], k_full = Full_Bloodfed[ind], k_non = Non_Bloodfed[ind], dose = rep(0,length(ind)), t = Time.min[ind]))

ind = which(BfData$Dose == 1.5)
bf.high <- mle2(minuslogl = multinomNLL1, start = c( beta0_part = 7, beta1_part = 0.15, beta0_full = 7, beta1_full = 0.15),
                data = list(N = Initial[ind], k_partial = Partial_Bloodfed[ind], k_full = Full_Bloodfed[ind], k_non = Non_Bloodfed[ind], dose = rep(0,length(ind)), t = Time.min[ind]))


# Likelihood-ratio tests on dose effects -------------------
df <- 4 - 2
L1 = -bf.nodose@details$value
L2 = -bf@details$value
teststat <- -2 * ( L1 - L2 )
pchisq(teststat, df = df, lower.tail = FALSE)

df <- (3*2) - 2
L1 = -bf.nodose@details$value
L2 = sum(-bf.control@details$value,
         -bf.low@details$value,
         -bf.high@details$value)
teststat <- -2 * (L1 - L2)
pchisq(teststat, df = df, lower.tail = FALSE)

# Estimate Exponentials of Biting Time (with dose effect)------------------------------------
# p-value dosage effect:
summary(bf)@coef[c(2,4),4]

Rate.full.bf=1/(exp(bf@coef[3] + bf@coef[4] * Dosage)) # full biting rate
Rate.part.bf=1/(exp(bf@coef[1] + bf@coef[2] * Dosage)) # partial biting rate
Rate.bf=Rate.full.bf + Rate.part.bf  # overall biting rate

Prob.full.bf <- sapply(Dosage,function(x)pexp(seq(0,max(Time.min)),rate=1/(exp(bf@coef[3] + bf@coef[4] * x) ))) 
Prob.part.bf <- sapply(Dosage,function(x)1 - (dpois(0,(1/(exp(bf@coef[1] + bf@coef[2] * x)) * seq(0,max(Time.min)))))) 
Prob.part.nofull.bf <- Prob.part.bf * (1 - Prob.full.bf)   
Prob.nobites.bf <- 1 - (Prob.full.bf + Prob.part.nofull.bf) 

# Estimate Exponentials of Biting Time (by dosage)------------------------------------
Rates.control <- get.feeding.rates(bf.control)
Rates.low <- get.feeding.rates(bf.low)
Rates.high <- get.feeding.rates(bf.high)

RATES.control <- get.feeding.rates.uncertainty(bf.control, num.draws = num.draws )
RATES.low <- get.feeding.rates.uncertainty(bf.low, num.draws = num.draws)
RATES.high <- get.feeding.rates.uncertainty(bf.high, num.draws = num.draws)

Probs.control <- get.feeding.probabilities.over.time(bf.control,time.vector = seq(0,max(Time.min)))
Probs.low <- get.feeding.probabilities.over.time(bf.low,time.vector = seq(0,max(Time.min)))
Probs.high <- get.feeding.probabilities.over.time(bf.high,time.vector = seq(0,max(Time.min)))

# similarly with uncertainty
PROBS.control <- get.feeding.probabilities.over.time.uncertainty(bf.control, time.vector = seq(0,max(Time.min)), num.draws = num.draws)
PROBS.low <- get.feeding.probabilities.over.time.uncertainty(bf.low, time.vector = seq(0,max(Time.min)), num.draws = num.draws)
PROBS.high <- get.feeding.probabilities.over.time.uncertainty(bf.high, time.vector = seq(0,max(Time.min)), num.draws = num.draws)

# Outcomes for paper ------------------------------------------------------
# time till blood feeding increased by:

(1/Rates.low['Rate.overall']) / (1/Rates.control['Rate.overall']) - 1  # FAR 1
HDIofMCMC( (1/RATES.low$Rate.overall) / (1/RATES.control$Rate.overall) - 1 ) # FAR 1

(1/Rates.high['Rate.overall']) / (1/Rates.control['Rate.overall']) - 1  # FAR 1.5
HDIofMCMC( (1/RATES.high$Rate.overall) / (1/RATES.control$Rate.overall) - 1 ) # FAR 1.5


# mostly result from full blood feeds:
(1/Rates.low['Rate.full']) / (1/Rates.control['Rate.full']) - 1  # FAR 1
HDIofMCMC( (1/RATES.low$Rate.full) / (1/RATES.control$Rate.full) - 1 ) # FAR 1

(1/Rates.high['Rate.full']) / (1/Rates.control['Rate.full']) - 1  # FAR 1.5
HDIofMCMC( (1/RATES.high$Rate.full) / (1/RATES.control$Rate.full) - 1 ) # FAR 1.5

# the partial biting rate increased by:
1 - (1/Rates.low['Rate.part']) / (1/Rates.control['Rate.part'])  # FAR 1
HDIofMCMC( 1 - (1/RATES.low$Rate.part) / (1/RATES.control$Rate.part) ) # FAR 1

1 - (1/Rates.high['Rate.part']) / (1/Rates.control['Rate.part'])   # FAR 1.5
HDIofMCMC( 1 - (1/RATES.high$Rate.part) / (1/RATES.control$Rate.part)) # FAR 1.5

(Rates.low['Rate.overall']) / (Rates.control['Rate.overall'])  # alpha FAR 1
HDIofMCMC ( RATES.low$Rate.overall / RATES.control$Rate.overall )
(Rates.high['Rate.overall']) / (Rates.control['Rate.overall']) # alpha FAR 1.5
HDIofMCMC ( RATES.high$Rate.overall / RATES.control$Rate.overall )

(Rates.low['Rate.full']) / (Rates.control['Rate.full'])  # alpha.gc FAR 1
HDIofMCMC ( RATES.low$Rate.full / RATES.control$Rate.full )
(Rates.high['Rate.full']) / (Rates.control['Rate.full']) # alpha.gc FAR 1.5
HDIofMCMC ( RATES.high$Rate.full / RATES.control$Rate.full )

# Table 2 -----------------------------------------------------------------

MV.Coefs = rmvn(n=1e4, mu = bf.control@coef, V = bf.control@vcov)
bf.control@coef['beta0_part']; HDIofMCMC(MV.Coefs[,1])
bf.control@coef['beta0_full']; HDIofMCMC(MV.Coefs[,3])

MV.Coefs = rmvn(n=1e4, mu = bf.low@coef, V = bf.low@vcov)
bf.low@coef['beta0_part']; HDIofMCMC(MV.Coefs[,1])
bf.low@coef['beta0_full']; HDIofMCMC(MV.Coefs[,3])

MV.Coefs = rmvn(n=1e4, mu = bf.high@coef, V = bf.high@vcov)
bf.high@coef['beta0_part']; HDIofMCMC(MV.Coefs[,1])
bf.high@coef['beta0_full']; HDIofMCMC(MV.Coefs[,3])


# Save --------------------------------------------------------------------

objects.to.save = apropos('bf', ignore.case = F)
save( list = objects.to.save, file = '../Input/BloodFeeding.RData')
save(list=ls(), file = '../Input/BloodFeedingWorkSpace.RData')

# Data Biting Probabilities with Data (Fig 2)---------------------------------------------------------------
transparency = 0.2
colors = Cs(black, darkorange, deeppink)
pdf('Fig2.pdf',width=6,height=3)

old.par <- par(mfrow=c(1,3),oma = c(5,4,4,2) + 0.1,
               mar = c(0,0,1,0) + 0.1)
plot(seq(0,nrow(Probs.control)-1),Probs.control[,'Prob.full'],type = "l" ,
     xlab="time in minutes", ylab="Probability ", col=colorlist[1], ylim=c(0,1), las = 1)
temp = apply(PROBS.control$PROB.full, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('black',transparency), border = alpha('black',transparency)) 

lines(seq(0,nrow(Probs.control)-1),Probs.low[,'Prob.full'],type = "l" ,
      col=colorlist[which(Dosage==1)])
temp = apply(PROBS.low$PROB.full, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('darkorange',transparency), border = alpha('darkorange',transparency)) 

lines(seq(0,nrow(Probs.control)-1),Probs.high[,'Prob.full'],type = "l" ,
      col=colorlist[which(Dosage==1.5)])
temp = apply(PROBS.high$PROB.full, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('deeppink',transparency), alpha('deeppink',transparency)) 

# add data
for (jj in c(0,1,1.5)){
  Ind = which(Dose == jj)
  Probs.temp = Full_Bloodfed[Ind]/Initial[Ind]
  CIs = binom.confint(Full_Bloodfed[Ind], Initial[Ind], methods = 'exact')[,c('lower', 'upper')] 
  points(Time.min[Ind],Probs.temp,pch=15, col = colors[which(Dosage.small==jj)])
  for( ii in 1:length(Time.min[Ind])){
    lines(rep(Time.min[Ind][ii],2),CIs[ii,],lwd = 1, col = colors[which(Dosage.small==jj)]) 
  }
}

mtext( bquote(bold("A")), side = 3, line = 0.1 ,at = 100,cex=1)
mtext(side = 3, text = 'Fully blood-fed', line = 1.1, cex = .9, font = 4)

plot(seq(0,nrow(Probs.control)-1),Probs.control[,'Prob.part.notfull'],type = "l" ,
     xlab="time in minutes", ylab="Probability ", col=colorlist[1], ylim=c(0,1), las = 1)
temp = apply(PROBS.control$PROB.part.notfull, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('black',transparency), border = alpha('black',transparency)) 

lines(seq(0,nrow(Probs.control)-1),Probs.low[,'Prob.part.notfull'],type = "l" ,
      col=colors[2])
temp = apply(PROBS.low$PROB.part.notfull, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('darkorange',transparency), border = alpha('darkorange',transparency)) 

lines(seq(0,nrow(Probs.control)-1),Probs.high[,'Prob.part.notfull'],type = "l" ,
      col=colors[3])
temp = apply(PROBS.high$PROB.part.notfull, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('deeppink',transparency), alpha('deeppink',transparency)) 


for (jj in c(0,1,1.5)){
  Ind = which(Dose == jj)
  Probs.temp = Partial_Bloodfed[Ind]/Initial[Ind]
  CIs = binom.confint(Partial_Bloodfed[Ind], Initial[Ind], methods = 'exact')[,c('lower', 'upper')] 
  points(Time.min[Ind],Probs.temp,pch=15, col = colors[which(Dosage.small==jj)])
  for( ii in 1:length(Time.min[Ind])){
    lines(rep(Time.min[Ind][ii],2),CIs[ii,],lwd = 1, col = colors[which(Dosage.small==jj)]) 
  }
}

mtext( bquote(bold("B")), side = 3, line = 0.1 ,at = 100,cex=1)
mtext(side = 3, text = 'Partially blood-fed', line = 1.1, cex = .9, font = 4)
legend('topright', c("control","low" ,"high"),
       col=colors, lty=c(rep(1,3)), cex=.9, bty='n')  

plot(seq(0,nrow(Probs.control)-1),Probs.control[,'Prob.nobites'],type = "l" ,
     xlab="time in minutes", ylab="Probability ", col=colorlist[1], ylim=c(0,1), las = 1)
temp = apply(PROBS.control$PROB.nobites, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('black',transparency), border = alpha('black',transparency)) 

lines(seq(0,nrow(Probs.control)-1),Probs.low[,'Prob.nobites'],type = "l" ,
      col=colors[2])
temp = apply(PROBS.low$PROB.nobites, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('darkorange',transparency), border = alpha('darkorange',transparency)) 

lines(seq(0,nrow(Probs.control)-1),Probs.high[,'Prob.nobites'],type = "l" ,
      col=colors[3])
temp = apply(PROBS.high$PROB.nobites, 1, HDIofMCMC)
polygon(c(seq(0,nrow(Probs.control)-1),rev(seq(0,nrow(Probs.control)-1)) ),c(temp[2,],rev(temp[1,])),col=alpha('deeppink',transparency), alpha('deeppink',transparency)) 


for (jj in c(0,1,1.5)){
  Ind = which(Dose == jj)
  Probs.temp = Non_Bloodfed[Ind]/Initial[Ind]
  CIs = binom.confint(Non_Bloodfed[Ind], Initial[Ind], methods = 'exact')[,c('lower', 'upper')] 
  points(Time.min[Ind],Probs.temp,pch=15, col = colors[which(Dosage.small==jj)])
  for( ii in 1:length(Time.min[Ind])){
    lines(rep(Time.min[Ind][ii],2),CIs[ii,],lwd = 1, col = colors[which(Dosage.small==jj)]) 
  }
}


mtext( bquote(bold("C")), side = 3, line = 0.1 ,at = 100,cex=1)
mtext(side = 3, text = 'Unfed', line = 1.1, cex = .9, font = 4)
mtext(side = 1, text = 'Time (min)', line = 2.25, cex = .9, outer = TRUE)
mtext(side = 2, text = 'Probability', line = 2.25, cex = .9, outer = TRUE)

dev.off()
