

# with 
# mu: probability of death if SR is there
# rho: probability of being repelled by SR 
# alpha: delay in biting by SR expressed as proportion of feeding cycle 

rel_foi = function(C,mu,rho,alpha,tau,q.t,q.u,g.t,g.u,g.tau,gc,alpha.gc){ 
  g = Death.Rate(0,0,tau,q.t,q.u,g.tau,g.u,g.t,mu)
  g.c = Death.Rate(C,rho,tau,q.t,q.u,g.tau,g.u,g.t,mu)
  a = Biting.rate(0,0,a.u,alpha,q.t,q.u)  
  a.c = Biting.rate(C,rho,a.u,alpha,q.t,q.u)
  gc.rate = Biting.rate(0,0,1/gc,alpha=alpha.gc,q.t,q.u)  
  gc.rate.c = Biting.rate(C,rho,1/gc,alpha=alpha.gc,q.t,q.u)
  m = gc.rate / g^2;   m.c = gc.rate.c / g.c^2
  foi = ((b * a^2 * m * c * X) / (g + a * c * X)) * exp(-g * n)
  foi.c = ((b * a.c^2 * m.c * c * X) / (g.c + a.c * c * X)) * exp(-g.c * n)
  if( any(m.c<0) ) { print('Warning: mosquito density is negative') }
  if( any(g.c<0) ) { print('Warning: mosquito death rate is negative') }
  if( any(a.c<0) ) { print('Warning: mosquito biting rate is negative') }
  if( any(gc.rate.c<0) ) { print('Warning: mosquito gonotrophic cycle is negative') }
  return(foi.c/foi)
}

# Nested functions

Death.Rate = function(C,rho,tau,q.t,q.u,g.tau,g.u,g.t,mu){
  Pi = Pies(C,rho,tau,q.t,q.u)
  g = Pi$pi.t * (mu * q.t + g.t) + 
    Pi$pi.u * g.u +
    Pi$pi.tau * g.tau
  return(g)
}

Pies = function(C,rho,tau,q.t,q.u){
  pi.tau = -( (tau*q.t*q.u)  / (q.t*(C-tau*q.u-1) + C*(rho-1)*q.u ))
  pi.u =  ((1-C)*pi.tau  ) / (tau * q.u) 
  pi.t =  ((1-rho)*C*pi.tau  ) / (tau * q.t)
  return(data.frame(pi.tau, pi.u, pi.t))
}

Biting.rate = function(C,rho,a.u,alpha,q.t,q.u){ 
  prob.fail.t = (q.t) / ((a.u * alpha) + q.t )  
  prob.fail.u = q.u / (a.u + q.u )
  D = C * rho + C * (1 - rho) * prob.fail.t + (1 - C) * prob.fail.u 
  geom.mean = D / ( 1 - D )           
  delta = ( tau * C * rho + 
              (1/q.t  + tau) * (C * (1-rho) * prob.fail.t) +
              (1/q.u + tau) * ((1-C) * prob.fail.u) ) / D
  a = 1 / ( geom.mean * delta + 1/a.u) 
  return(a)
  
}


