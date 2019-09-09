

# Defaults / Baseline  ----------------------------------------------------------------
# all quantities in days

n = 14                 # eip 
b = .5                 # probability of infection mosquito to human
c = .5                 # probability of infection human to mosquito
X = .125               # human prevalence
gc = 4                 # gonotrophic cycle
g.u = 1-0.82           # mortality
g.tau.scaler = 1       # outdoor hazard fraction proportional to indoor mortality
a.u = 0.76             # biting rate
q.u = 1/2              # 1 over residence time
tau.scaler = 0.3       # proportion of time spent outdoors
high.dose = 1.5

# Scalers from experiments ---------------------------------------------------------------
mu = exp(kd.dose@coef[1] + kd.dose@coef[2])  -    exp(kd.dose@coef[1])         # instantaneous death rate due to product
mu.h = exp(kd.dose@coef[1] + kd.dose@coef[2] * high.dose)  -    exp(kd.dose@coef[1]) 
aft.g <- exp(exp.mod.con$coefficients[4])                   # fraction of waiting time till death upon treatment
aft.g.h <- exp(exp.mod.con$coefficients[6] ) 

Rate.full.bf= c(1/(exp(bf.control@coef[3] )), 1/(exp(bf.low@coef[3] )), 1/(exp(bf.high@coef[3] )) ); names(Rate.full.bf) = Cs(control, low, high)
Rate.part.bf=c(1/(exp(bf.control@coef[1] )), 1/(exp(bf.low@coef[1] )), 1/(exp(bf.high@coef[1] )) ); names(Rate.part.bf) = Cs(control, low, high)
Rate.bf=Rate.full.bf + Rate.part.bf  # overall biting rate

alpha <- Rate.bf['low'] / Rate.bf['control']   # FAR 1
alpha.h <- Rate.bf['high'] / Rate.bf['control']   # FAR 1.5
alpha.gc <- Rate.full.bf['low'] / Rate.full.bf['control']
alpha.gc.h <- Rate.full.bf['high'] / Rate.full.bf['control']

rho = median(Repellency[2,], na.rm = TRUE) # repellency effect
rho.h = median(Repellency[3,], na.rm = TRUE)
q.t.scaler = median(Expellency[2,], na.rm = TRUE)
q.t.scaler.h = median(Expellency[3,], na.rm = TRUE)

# Derived parameters -----------------------------------------------------------------
# DONE IN THE FUNCTIONS
g.t = g.u / aft.g 
g.t.h = g.u / aft.g.h
g.tau = g.u * g.tau.scaler                              # death rate during transit
q.t = q.u  * q.t.scaler
q.t.h = q.u  * q.t.scaler.h
tau = tau.scaler*(1/q.u)  # (1-0.77)/0.77 for Thailand, idem with 0.9 for Puerto Rico

# Derive uncertainties ----------------------------------------------------
Mu.coefficients = rmvn(n= num.draws, mu = kd.dose@coef, V = kd.dose@vcov)
Mu = exp(Mu.coefficients[,1] + Mu.coefficients[,2]) -  exp(Mu.coefficients[,1])
Mu.h = exp(Mu.coefficients[,1] + Mu.coefficients[,2] * high.dose) -  exp(Mu.coefficients[,1])

Aft.g.coefficients <- rmvn(n= num.draws, mu = exp.mod.con$coefficients, V = exp.mod.con$var)
AFT.g.exp <- exp(Aft.g.coefficients[,4]) 
AFT.g.exp.h <- exp(Aft.g.coefficients[,6]) 
G.T = g.u / AFT.g.exp
G.T.h = g.u / AFT.g.exp.h

Bf.full.coefficients = rmvn(n= num.draws, mu = bf.control@coef[3:4], V = bf.control@vcov[3:4,3:4]); 
RATE.full.bf=1/(exp(Bf.full.coefficients[,1]))
Bf.part.coefficients = rmvn(n= num.draws, mu = bf.control@coef[1:2], V = bf.control@vcov[1:2,1:2])
RATE.part.bf=1/(exp(Bf.part.coefficients[,1]))
Bf.full.coefficients = rmvn(n= num.draws, mu = bf.low@coef[3:4], V = bf.low@vcov[3:4,3:4]); 
RATE.full.bf=cbind(RATE.full.bf, 1/(exp(Bf.full.coefficients[,1])))
Bf.part.coefficients = rmvn(n= num.draws, mu = bf.low@coef[1:2], V = bf.low@vcov[1:2,1:2])
RATE.part.bf=cbind(RATE.part.bf, 1/(exp(Bf.part.coefficients[,1])))
Bf.full.coefficients = rmvn(n= num.draws, mu = bf.high@coef[3:4], V = bf.high@vcov[3:4,3:4]); 
RATE.full.bf=cbind(RATE.full.bf, 1/(exp(Bf.full.coefficients[,1]))); colnames(RATE.full.bf) = Cs(control, low, high)
Bf.part.coefficients = rmvn(n= num.draws, mu = bf.high@coef[1:2], V = bf.high@vcov[1:2,1:2])
RATE.part.bf=cbind(RATE.part.bf, 1/(exp(Bf.part.coefficients[,1]))); colnames(RATE.part.bf) = Cs(control, low, high)

RATE.bf=RATE.full.bf + RATE.part.bf; colnames(RATE.bf) = Cs(control, low, high)  # overall biting rates with uncertainty

Alpha <- (RATE.bf[,'low']) / (RATE.bf[,'control'])  # FAR 1
Alpha.gc <- (RATE.full.bf[,'low']) / (RATE.full.bf[,'control'])  
Alpha.h <- (RATE.bf[,'high']) / (RATE.bf[,'control'])  # FAR 1.5
Alpha.gc.h <- (RATE.full.bf[,'high']) / (RATE.full.bf[,'control'])  

Rho = Repellency[2,complete.cases(Repellency[2,])]
Rho.h = Repellency[3,complete.cases(Repellency[3,])]

Q.t = Expellency[2,complete.cases(Expellency[2,])]  # q.t.scaler uncertainty
Q.t.h = Expellency[3,complete.cases(Expellency[3,])]


