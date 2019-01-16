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

# Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location (if working in RStudio)
getwd()

# Data loading and processing ---------------------------------------------
BfData <- read.csv('../Data/SR_Bloodfeeding.csv')
Time.min <- BfData$Time * 60 + 10 # adding exposure time to time points
BfData <- data.frame(BfData, Time.min)
attach(BfData)

# Plotting parameters -----------------------------------------------------
Dosage <- c(seq(0,2,by=0.25))
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
#   test = prob_full + prob_part_nofull + prob_non
#   print(test)

    -sum(dmultinom(c(k_partial,k_full,k_non),prob=c(prob_part_nofull,
                                                  prob_full,
                                                  prob_non),
                 size=NULL,log=TRUE))
}

parnames(multinomNLL1) <- c("beta0_part", "beta1_part", "beta0_full", "beta1_full")

# Fit Multinomial Models ------------------------------------------------------

bf <- mle2(minuslogl = multinomNLL1, start = c( beta0_part = 7, beta1_part = 0.15, beta0_full = 7, beta1_full = 0.15),
           data = list(N = Initial, k_partial = Partial_Bloodfed, k_full = Full_Bloodfed, k_non = Non_Bloodfed, dose = Dose, t = Time.min))


# Estimate Exponentials of Biting Time ------------------------------------
# p-value dosage effect:
summary(bf)@coef[c(2,4),4]

Rate.full.bf=1/(exp(bf@coef[3] + bf@coef[4] * Dosage)) # full biting rate
Rate.part.bf=1/(exp(bf@coef[1] + bf@coef[2] * Dosage)) # partial biting rate
Rate.bf=Rate.full.bf + Rate.part.bf  # overall biting rate

Prob.full.bf <- sapply(Dosage,function(x)pexp(seq(0,max(Time.min)),rate=1/(exp(bf@coef[3] + bf@coef[4] * x) ))) 
Prob.part.bf <- sapply(Dosage,function(x)1 - (dpois(0,(1/(exp(bf@coef[1] + bf@coef[2] * x)) * seq(0,max(Time.min)))))) 
Prob.part.nofull.bf <- Prob.part.bf * (1 - Prob.full.bf)   
Prob.nobites.bf <- 1 - (Prob.full.bf + Prob.part.nofull.bf)  

# Outcomes for paper ------------------------------------------------------
# time till blood feeding increased by:
(1/Rate.bf[which(Dosage==1)]) / (1/Rate.bf[which(Dosage==0)]) - 1  # FAR 1
(1/Rate.bf[which(Dosage==1.5)]) / (1/Rate.bf[which(Dosage==0)]) - 1  # FAR 1.5
# mostly result from full blood feeds:
(1/Rate.full.bf[which(Dosage==1)]) / (1/Rate.full.bf[which(Dosage==0)]) - 1  # FAR 1
(1/Rate.full.bf[which(Dosage==1.5)]) / (1/Rate.full.bf[which(Dosage==0)]) - 1  # FAR 1.5
# the partial biting rate increased by:
1 - (1/Rate.part.bf[which(Dosage==1)]) / (1/Rate.part.bf[which(Dosage==0)])   # FAR 1
1 - (1/Rate.part.bf[which(Dosage==1.5)]) / (1/Rate.part.bf[which(Dosage==0)])   # FAR 1.5

(Rate.bf[which(Dosage==1)]) / (Rate.bf[which(Dosage==0)])   # alpha FAR 1
(Rate.bf[which(Dosage==1.5)]) / (Rate.bf[which(Dosage==0)])   # alpha FAR 1

(Rate.full.bf[which(Dosage==1)]) / (Rate.full.bf[which(Dosage==0)])  # alpha.gc FAR 1
(Rate.full.bf[which(Dosage==1.5)]) / (Rate.full.bf[which(Dosage==0)])  # alpha.gc FAR 1

save( 'bf', file = 'BloodFeeding.RData')
save(list=ls(), file = 'BloodFeedingWorkSpace.RData')

# Data Biting Probabilities with Data (Fig 2)---------------------------------------------------------------

pdf('Fig2.pdf',width=6,height=3)

old.par <- par(mfrow=c(1,3),oma = c(5,4,4,2) + 0.1,
               mar = c(0,0,1,0) + 0.1)
plot(seq(0,dim(Prob.full.bf)[1]-1),Prob.full.bf[,1],type = "l" ,
     xlab="time in minutes", ylab="Probability ", col=colorlist[1], ylim=c(0,1), las = 1)
for (i in 2:length(Dosage)){
  lines(seq(0,dim(Prob.full.bf)[1]-1),Prob.full.bf[,i], col=colorlist[i])
}
Ind = which(Dose == 0)
points(Time.min[Ind],Full_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==0)])
Ind = which(Dose == 1)
points(Time.min[Ind],Full_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==1)])
Ind = which(Dose == 1.5)
points(Time.min[Ind],Full_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==1.5)])

mtext( bquote(bold("A")), side = 3, line = 0.1 ,at = 100,cex=1)
mtext(side = 3, text = 'Fully blood-fed', line = 1.1, cex = .9, font = 4)

plot(seq(0,dim(Prob.part.nofull.bf)[1]-1),Prob.part.nofull.bf[,1],type = "l" ,
     xlab="time in minutes", ylab = " ",col=colorlist[1], ylim=c(0,1),yaxt='n', las = 1)
for (i in 2:length(Dosage)){
  lines(seq(0,dim(Prob.part.nofull.bf)[1]-1),Prob.part.nofull.bf[,i], col=colorlist[i])
}
Ind = which(Dose == 0)
points(Time.min[Ind],Partial_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==0)])
Ind = which(Dose == 1)
points(Time.min[Ind],Partial_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==1)])
Ind = which(Dose == 1.5)
points(Time.min[Ind],Partial_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==1.5)])

mtext( bquote(bold("B")), side = 3, line = 0.1 ,at = 100,cex=1)
mtext(side = 3, text = 'Partially blood-fed', line = 1.1, cex = .9, font = 4)
legend('topright', c("control","SR dose 0.5" ,"SR dose 0.75", "SR dose 1.00", "SR dose 1.25", "SR dose 1.5"),
       col=colorlist[c(seq(6))], lty=c(rep(1,6)), cex=.9, bty='n')  

plot(seq(0,dim(Prob.nobites.bf)[1]-1),Prob.nobites.bf[,1],type = "l" ,
     xlab="time in minutes", ylab = " ",col=colorlist[1], ylim=c(0,1),yaxt='n', las = 1)
for (i in 2:length(Dosage)){
  lines(seq(0,dim(Prob.nobites.bf)[1]-1),Prob.nobites.bf[,i], col=colorlist[i])
}
Ind = which(Dose == 0)
points(Time.min[Ind],Non_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==0)])
Ind = which(Dose == 1)
points(Time.min[Ind],Non_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==1)])
Ind = which(Dose == 1.5)
points(Time.min[Ind],Non_Bloodfed[Ind]/Initial[Ind],pch=15, col = colorlist[which(Dosage==1.5)])

mtext( bquote(bold("C")), side = 3, line = 0.1 ,at = 100,cex=1)
mtext(side = 3, text = 'Unfed', line = 1.1, cex = .9, font = 4)
mtext(side = 1, text = 'Time (min)', line = 2.25, cex = .9, outer = TRUE)
mtext(side = 2, text = 'Probability', line = 2.25, cex = .9, outer = TRUE)

dev.off()
