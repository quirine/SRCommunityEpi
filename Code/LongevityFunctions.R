

get.TTE <- function(Data,dosage){
  
  TTE=numeric(); ID=numeric(); Event=numeric(); Treat =numeric(); Dose =numeric(); Group =numeric(); PhysState =character()
  
  for( j in 1:ng){
    TempDiff <- diff( Data[,1+j],lag = 1)
    if(any(TempDiff>0)) message('Warning: There are zombie mosquitoes in your data!')
    TempDiff <- abs( TempDiff)
    TempN <- sum(TempDiff)
    Nums <- array(data = 0, dim= TempN)
    m=1
    for( i in 1:length(TempDiff)){
      if(TempDiff[i] > 0){
        Nums[m:(m+TempDiff[i]-1)] <- rep.int(i,TempDiff[i]) 
        m <- m+TempDiff[i]
      }
    }
    message(c(TempN,"",i))
    TTE <-c(TTE,Nums); ID <- c(ID,seq(TempN)); Event <- c(Event, rep(1,TempN)); Treat <- c(Treat, rep(0,TempN))
    Dose <- c(Dose, rep(0.00,TempN)); Group <- c(Group, rep(j,TempN)); PhysState <- c(PhysState, rep('sugar fed',TempN))
  }
  
  for( j in 1:ng){
    TempDiff <- diff( Data[,9+j],lag = 1)
    if(any(TempDiff>0)) message('Warning: There are zombie mosquitoes in your data!')
    TempDiff <- abs( TempDiff)
    TempN <- sum(TempDiff)
    Nums <- array(data = 0, dim= TempN)
    m=1
    for( i in 1:length(TempDiff)){
      if(TempDiff[i] > 0){
        Nums[m:(m+TempDiff[i]-1)] <- rep.int(i,TempDiff[i]) 
        m <- m+TempDiff[i]
      }
    }
    message(c(TempN,"",i))
    TTE <-c(TTE,Nums); ID <- c(ID,seq(TempN)); Event <- c(Event, rep(1,TempN)); Treat <- c(Treat, rep(1,TempN))
    Dose <- c(Dose, rep(dosage,TempN)); Group <- c(Group, rep(6+j,TempN)); PhysState <- c(PhysState, rep('sugar fed',TempN))
  }
  DataFrame <- data.frame(ID,TTE,Event,Treat,Dose,Group,PhysState)
  return(DataFrame)
} 

binom.instant.kd = function(params, k,N,dose){ 
  beta0 <- params[1]
  beta1 <- params[2]
  rate = exp(beta0 + beta1 * dose )
  probability = pexp(1,rate)
  -sum(dbinom(k,prob=probability,size=N,log=TRUE))
}

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}