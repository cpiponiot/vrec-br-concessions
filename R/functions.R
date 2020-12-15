# volume prediction depending on stand maturity (t) and other model parameters
volume = function(t, ag , am, bg, bm, th, pdef = 0) {
  (ag / th * (1 - (th * exp(-bg * t) - bg * exp(-th * t)) / (th - bg)) - am /
     th * (1 - (th * exp(-bm * t) - bm * exp(-th * t)) / (th - bm))) * (1 - pdef)
}

# prediction of volume loss depending on extracted volume (m3/ha)
deltaVPrediction <- function(V0, Vext, om0, psi, e) {
 
  mu = om0 ^ (1 - psi)
  Var = (1 - mu) * mu ^ 2 * e
  
  alpha =   (mu / Var - 1 / mu) * mu ^ 2
  beta = alpha * (1 / mu - 1)
  
  # ome = proportion of commercial volume in volume loss due to logging operations
  ome = apply(cbind(alpha, beta, om0), 1, function(X) {
    if (X[3] %in% c(0,1)) 
      return(X[3]) 
    else 
      return(rtrunc(1, "beta", shape1 = X[1], shape2 = X[2],  a = X[3]))
  })  
  
  deltaV = apply(cbind(Vext / ome, V0), 1, min)
  deltaV[om0 == 0] <- 0
  return(deltaV)

}

# prediction of post-logging maturity t1 
# solve equation: volume(t1) = volume(t0) - deltaV
t0Prediction <- function(t0, deltaV, aG, aM, bG, bM, theta, pdef = 0) {
  equ = function(x)
    volume(x, aG, aM, bG, bM, theta, pdef) - (volume(t0, aG, aM, bG, bM, theta, pdef) - deltaV)
  return(uniroot(equ, c(0, t0))$root)
}
# same function applied to a parameter vector X
t0Prediction2 <- function(X) {
  t0     <- X[1]
  deltaV <- X[2]
  aG     <- X[3]
  aM     <- X[4]
  bG     <- X[5]
  bM     <- X[6]
  theta  <- X[7]
  pdef   <- X[8]
  equ = function(x)
    volume(x, aG, aM, bG, bM, theta, pdef) - (volume(t0, aG, aM, bG, bM, theta, pdef) - deltaV)
  return(uniroot(equ, c(0, t0))$root)
}

# evolution of timber volume proportion in time
#### Arguments ####
## tvec: maturity vector
## t1: first maturity
## om1: first (post-logging) proportion of timber species in the total volume
## omR: proportion of timber species in recruits' volume
## ag, am, bg, bm, th: volume model parameters
## int, slope: logging precision model parameters (intercept, slope)

omega_t = function(tvec, t1, om1, omR, ag, am, bg, bm, th, int, slope){
  
  t = seq.int(t1,max(tvec)); ord_t = floor(tvec - t1 + 1)
  
  V = ag/th*(1-(th*exp(-bg*t)-bg*exp(-th*t))/(th-bg))-am/th*(1-(th*exp(-bm*t)-bm*exp(-th*t))/(th-bm))
  
  dVG = ag*bg/(th-bg)*(exp(-bg*t)-exp(-th*t))+am*(1-(th*exp(-bm*t)-bm*exp(-th*t))/(th-bm))
  
  dVM = am*(1-exp(-bm*t))
  
  pR = 1/(1+exp(-(int+slope*log(V))))
  
  om = om1
  
  for (i in 1:(length(t)-1)) {
    om2 = min((om[i]*V[i] + dVG[i]*(pR[i]*omR + (1-pR[i])*om[i]) - dVM[i]*om[i])/V[i+1], 1)
    om = c(om,om2)
  }
  
  return(om[ord_t])
}