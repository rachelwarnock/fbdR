
#' @export
fbdskyline.probability<-function(frs,b,d,s,rho,k,int.min){
  frs<-frs
  numFossils<-k

  lambda<<-b
  mu<<-d
  psi<<-s
  rho<<-rho
  intervals.min<<-int.min

  num.lineages = length(frs$sp)
  num.extinct = num.lineages - sum(frs$extant)
  ot = max(frs$bi)
  oi = frs$bint[which(frs$bi == ot)]

  lk = num.lineages * log(rho[1]) - log(1 - fbdSkylineP(oi, ot))

  # fossils
  for(i in 1:length(intervals.min)){
    if(k[i] > 0)
      lk = lk + k[i] * log(lambda[i])
  }

  # extinction events
  for(i in 1:length(frs$sp)){
    d = frs$di[i]
    if(d != 0){
      j = frs$dint[i]
      lk = lk + log(mu[j])
    }
  }

  # speciation events
  for(i in 1:length(frs$sp)){
    b = frs$bi[i]
    if(b != ot){
      j = frs$bint[i]
      lk = lk + log(lambda[j])
    }
  }

  # range probabilities
  for(i in 1:length(frs$sp)){
    g = frs$gamma[i]
    b = frs$bi[i]
    d = frs$di[i]
    o = frs$oi[i]
    bint = frs$bint[i]
    dint = frs$dint[i]
    oint = frs$oint[i]

    lk = lk + log(g) + fbdSkylineQtildaLog(oint, o) - fbdSkylineQtildaLog(dint, d)
    lk = lk + fbdSkylineQLog(bint, b) - fbdSkylineQLog(oint, o)

  }

  return(lk)
}

#' @export
recount.k<-function(fossils, intervals){
  k = c()
  for(i in 1:length(intervals)){
    k = c(k, 0)
  }
  for(i in 1:length(f$h)){
    if(f$h[i] != 0){
      j = assign.interval(intervals, f$h[i])
      k[j] = k[j] + 1
    }
  }
  return(k)
}

#' @export
reassign.intervals<-function(frs,intervals){
  frs = cbind(frs, bint=-1, dint=-1, oint=-1)
  # birth
  for(i in 1:length(frs$bi)){
    frs$bint[i] = assign.interval(intervals, frs$bi[i])
  }
  # death
  for(i in 1:length(frs$di)){
    frs$dint[i] = assign.interval(intervals, frs$di[i])
  }
  # FA
  for(i in 1:length(frs$oi)){
    frs$oint[i] = assign.interval(intervals, frs$oi[i])
  }
  return(frs)
}

assign.interval<-function(intervals, t){
  i = -1
  for(j in 1:length(intervals)){
    if(t >= intervals[j])
      i = j
  }
  return(i)
}

fbdSkylineA<-function(i){

  Ai = sqrt ( (lambda[i] - mu[i] - psi[i]) * (lambda[i] - mu[i] - psi[i]) + 4 * lambda[i] * psi[i] )
  return(Ai)

}

fbdSkylineB<-function(i){

  Ai = fbdSkylineA(i)
  ti = intervals.min[i]

  Bi = ( (1 - 2 * (1 - rho[i]) * fbdSkylineP(i-1, ti)) * lambda[i] + mu[i] + psi[i] ) / Ai

  return(Bi)

}

fbdSkylineBLog<-function(i){

  Ai = fbdSkylineA(i)
  ti = intervals.min[i]

  a = (1 - (2 * (1 - rho[i]) * fbdSkylineP(i-1, ti))) * lambda[1]
  b = mu[i]
  c = psi[i]
  x = max(a, b, c)
  Bi = log(x) + log(a/x + b/x + c/x) - log(Ai)

  return(Bi)

}

fbdSkylineP<-function(i, t){

  if(i == 0 || t == 0)
    return(1)

  Ai = fbdSkylineA(i)
  Bi = fbdSkylineB(i)
  ti = intervals.min[i]

  p1 = lambda[i] + mu[i] + psi[i]
  p2 = Ai * ( exp(Ai * (t - ti)) * (1 + Bi) - (1 - Bi) ) / ( exp(Ai * (t - ti)) * (1 + Bi) + (1 - Bi) )
  p = (p1 - p2) / 2 * lambda[i]

  return(p)
}

fbdSkylineQ<-function(i, t){

  Ai = fbdSkylineA(i)
  Bi = fbdSkylineB(i)
  ti = intervals.min[i]

  q = (4 * exp(Ai * (t - ti))) / (( exp(-Ai * (t - ti)) * (1 - Bi) + (1 + Bi) ) * ( exp(-Ai * (t - ti)) * (1 - Bi) + (1 + Bi) ) )

  return(q)
}

# this does work but shouldn't
fbdSkylineQLog<-function(i, t){

  Ai = fbdSkylineA(i)
  Bi = fbdSkylineB(i)
  ti = intervals.min[i]

  q = (log(4) + (Ai * (t - ti))) - log(( exp(-Ai * (t - ti)) * (1 - Bi) + (1 + Bi) ) * ( exp(-Ai * (t - ti)) * (1 - Bi) + (1 + Bi) ) )

  return(q)
}

fbdSkylineQtilda<-function(i, t){

  Ai = fbdSkylineA(i)
  Bi = fbdSkylineB(i)
  #ti = intervals.min[i]

  q1a = 4 * exp(-t * (lambda[i] + mu[i] + psi[i]) ) * exp(-t * Ai)
  q1b = 4 * exp(-t * Ai) + (1 - Bi * Bi) * (1 - exp(-t * Ai)) *  (1 - exp(-t * Ai))
  q2a = (1 + Bi) * exp(-t * Ai) + (1 - Bi)
  q2b = (1 - Bi) * exp(-t * Ai) + (1 + Bi)

  qt = sqrt( (q1a/q1b) * (q2a/q2b) )
  return(qt)
}

fbdSkylineQtildaLog<-function(i, t){

  Ai = fbdSkylineA(i)
  Bi = fbdSkylineB(i)

  q1aLog = log(4) + (-t * (lambda[i] + mu[i] + psi[i]) ) + (-t * Ai)
  q1b = 4 * exp(-t * Ai) + (1 - Bi * Bi) * (1 - exp(-t * Ai)) *  (1 - exp(-t * Ai))
  q2a = (1 + Bi) * exp(-t * Ai) + (1 - Bi)
  q2b = (1 - Bi) * exp(-t * Ai) + (1 + Bi)

  qt = 0.5 * (q1aLog - log(q1b) + log(q2a) - log(q2b))
  return(qt)
}

revbayes_example<-function(i, t){

  if(t == 0)
    return(1.0)

  b = lambda[i]
  d = mu[i]
  f = psi[i]
  r = rho[i]
  ti = intervals.min[i]

  diff = b - d - f
  bp = b*f
  dt = t - ti

  A = sqrt( diff*diff + 4*bp )
  B = ( (1.0 - 2.0*(1.0-r)*revbayes_example(i-1, ti))*b + d + f ) / A

  e = exp(A*dt)
  tmp = b + d + f - A*(e*(1.0+B)-(1.0-B)) / (e*(1.0+B)+(1.0-B))

  return(tmp / 2*b )
}


