sr.markov.model <- function (t, x, parms) {
  
  D <- x[1];  U <- x[2];  T <- x[3];   
  q.u <- parms[1];  q.t <- parms[2];  tau <- parms[3];    rho <- parms[4];     C <- parms[5]
  
  # differential equations
  dD <- ((1/tau)*(rho*C-1)) * D   + q.u * U     + q.t * T 
  dU <- ( (1/tau)*(1-C)) * D      - q.u * U             
  dT <- ((1/tau)*((1-rho)*C)) * D               - q.t * T
    
  list(c(dD,dU,dT))
}
