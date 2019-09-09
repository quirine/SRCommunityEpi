
get.feeding.rates <- function(bbmle.model){
  Rate.full=1/(exp(bbmle.model@coef[3])) # full biting rate
  Rate.part=1/(exp(bbmle.model@coef[1])) # partial biting rate
  Rate.overall=Rate.full + Rate.part  # overall biting rate
  
  out = c( Rate.full, Rate.part, Rate.overall )
  names(out) = Cs( Rate.full, Rate.part, Rate.overall )
  return(out)
  
}

get.feeding.rates.uncertainty <- function(bbmle.model, num.draws){
  COEFS = rmvn(n= num.draws, mu = bbmle.model@coef[c(1,3)], V = bbmle.model@vcov[c(1,3),c(1,3)])
  RATE.full=1/(exp(COEFS[,2])) # full biting rate
  RATE.part=1/(exp(COEFS[,1])) # partial biting rate
  RATE.overall=RATE.full + RATE.part  # overall biting rate
  
  out = list( RATE.full, RATE.part, RATE.overall )
  names(out) = Cs( Rate.full, Rate.part, Rate.overall )
  return(out)
  
}

get.feeding.probabilities.over.time <- function(bbmle.model, time.vector){
  Prob.full <- feed.probs.full(bbmle.model@coef[3], time.vector)
  Prob.part <- feed.probs.part(bbmle.model@coef[1], time.vector)
  Prob.part.notfull <- Prob.part * (1 - Prob.full)   
  Prob.nobites <- 1 - (Prob.full + Prob.part.notfull)
  
  out = cbind( Prob.full, Prob.part, Prob.part.notfull, Prob.nobites)
  return(out)
}

get.feeding.probabilities.over.time.uncertainty <- function(bbmle.model, time.vector, num.draws){
  COEFS.full = rmvn(n= num.draws, mu = bbmle.model@coef[3:4], V = bbmle.model@vcov[3:4,3:4])
  COEFS.part = rmvn(n= num.draws, mu = bbmle.model@coef[1:2], V = bbmle.model@vcov[1:2,1:2])
  PROB.full <- sapply(COEFS.full[,1],feed.probs.full, time.vector)
  PROB.part <- sapply(COEFS.part[,1],feed.probs.part, time.vector)
  PROB.part.notfull <- PROB.part * (1 - PROB.full)   
  PROB.nobites <- 1 - (PROB.full + PROB.part.notfull)
  out <- list(COEFS.full, COEFS.part, PROB.full, PROB.part, PROB.part.notfull, PROB.nobites)
  names(out) <- Cs(COEFS.full, COEFS.part, PROB.full, PROB.part, PROB.part.notfull, PROB.nobites)
  return(out)
}

feed.probs.full <- function(coefs,time.vector){
  Prob.full <- pexp(time.vector,rate=1/(exp(coefs) ))
}

feed.probs.part <- function(coefs, time.vector){
  Prob.part <- 1 - (dpois(0,(1/(exp(coefs )) * time.vector)))  
}


