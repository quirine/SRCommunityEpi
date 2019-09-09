rm(list=ls())

# load libraries ----------------------------------------------------------
library(deSolve)

# set wd ------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location (if working in RStudio)

# load source files -------------------------------------------------------
source('SR.Markov.Models.R')

# set parameters ----------------------------------------------------------
delta = 0.01
Init <- rep(1/3, 3)
Time <- seq(0,50,by=delta)
q.u <- 0.3
q.t <- 0.1
tau <- 1
rho <- .2
C <- .5

params = c(q.u, q.t, tau, rho, C)


# Run model ---------------------------------------------------------------

out = ode(Init,Time,sr.markov.model,parms=params,method = 'ode45')


# Calculate equilibrium states --------------------------------------------

Pi.tau = -( (tau*q.t*q.u)  / (q.t*(C-tau*q.u-1) + C*(rho-1)*q.u ))
Pi.u =  ((1-C)*Pi.tau  ) / (tau * q.u) 
Pi.t =  ((1-rho)*C*Pi.tau  ) / (tau * q.t)

# Plot --------------------------------------------------------------------
plot(out[,1], out[,2], col = 'black', type = 'l',ylim=c(0,1), xlab ='', ylab='', lwd = 2)
lines(out[,1], out[,3], col = 'green', lwd = 2)
lines(out[,1], out[,4], col = 'red', lwd = 2)
abline(h = Pi.tau, lty = 3, col = 'black', lwd = 2)
abline(h = Pi.u, lty = 3, col = 'green', lwd = 2)
abline(h = Pi.t, lty = 3, col = 'red', lwd = 2)
mtext(side =1, text = "Time", line = 2.25)
mtext(side =2, text = "Probability", line = 2.25)



