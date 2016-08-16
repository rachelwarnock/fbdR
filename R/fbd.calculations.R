
#' @export
recount.gamma<-function(frs){
  #frs<-fossilRanges
  frs<-frs

  ot=max(frs$bi)

  for(f in 1:length(frs$bi)){
    bf=frs$bi[f]
    df=frs$di[f]
    g=0
    for(j in 1:length(frs$bi)){
      bj=frs$bi[j]
      dj=frs$di[j]

      if(bj > bf & dj < bf)
        g=g+1
    }
    if(bf==ot)
      g=g+1

    frs$gamma[f]=g
  }
  #fossilRanges <<- frs
  return(frs)
  #eof
}

#' @export
recount.extant<-function(frs){
  #frs<-fossilRanges
  frs<-frs

  for(f in 1:length(frs$bi)){
    df=frs$di[f]

    if(df==0)
      frs$extant[f]=1
    else
      frs$extant[f]=0
  }
  #fossilRanges <<- frs
  return(frs)
  #eof
}

#### probability functions

#' @export
fbd.probability<-function(frs,b,d,s,k,rho=1,complete=F){
  frs<-frs
  lambda<<-b
  mu<<-d
  psi<<-s
  rho<<-rho
  numFossils<-k
  complete<-complete

  ot = max(frs$bi) # define the origin
  extinctLineages = length(frs$extant) - sum(frs$extant)

  pr = numFossils*log(psi)
  pr = pr + extinctLineages*log(mu)
  pr = pr - log(lambda * (1 - fbdPfxn(ot) ) )

  # for complete sampling (b_i = o_i)
  if(complete)
    rp=sum(unlist(lapply(1:length(frs$bi),function(x){ rangePrComplete(frs$gamma[x],frs$bi[x],frs$di[x]) })))
  else
    rp=sum(unlist(lapply(1:length(frs$bi),function(x){ rangePr(frs$gamma[x],frs$bi[x],frs$di[x],frs$oi[x]) })))

  pr = pr + rp
  return(pr)

}

fbdC1fxn<-function(){

  c1 = abs ( sqrt( (lambda - mu - psi)^2 + (4*lambda*psi) ) )

  return(c1)
}

fbdC2fxn<-function(){

  c1 = fbdC1fxn()
  c2 = -(lambda - mu - (2*lambda*rho) - psi ) / c1

  return(c2)

}

fbdC3fxn<-function(){

  c3 = lambda * (-psi+rho*(mu+lambda*(-1+rho)+psi))

  return(c3)

}

fbdC4fxn<-function(){

  c1 = fbdC1fxn()
  c3 = fbdC3fxn()

  c4 = c3/(c1^2)

  return(c4)
}

fbdPfxn<-function(t){
  t<-t

  c1 = fbdC1fxn()
  c2 = fbdC2fxn()

  # p = 1 + ( (-(lambda-mu-psi) +  (c1 * ( ( exp(-c1*t)*(1-c2)-(1+c2) ) / ( exp(-c1*t)*(1-c2)+(1+c2) ) ) ) ) / (2*lambda) )
  # p = ( ((lambda+mu+psi) +  (c1 * ( ( exp(-c1*t)*(1-c2)-(1+c2) ) / ( exp(-c1*t)*(1-c2)+(1+c2) ) ) ) ) / (2*lambda) )
  p = ( ((lambda+mu+psi) +  (c1 * ( ( exp(-c1*t)*(1-c2)/(1+c2) - 1 ) / ( exp(-c1*t)*(1-c2)/(1+c2) + 1 ) ) ) ) / (2*lambda) )

  return(p)
}

fbdQTildaFxnLog<-function(t){
  t<-t

  c1 = fbdC1fxn()
  c2 = fbdC2fxn()
  c4 = fbdC4fxn()

  f1aLog = -t*(lambda+mu+psi+c1) #log f1a
  f1b = c4 * (1-exp(-t*c1))^2 - exp(-t*c1)
  f2a = ((1-c2) * (2*c4*(exp(-t*c1)-1))+exp(-t*c1)*(1-c2^2))
  f2b = ((1+c2) * (2*c4*(exp(-t*c1)-1))+exp(-t*c1)*(1-c2^2))

  f = f1aLog + log(-1/f1b * (f2a/f2b) )

  q = 0.5*f + log(rho)

  return(q)

}

fbdQfxnLog<-function(t){
  t<-t

  c1 = fbdC1fxn()
  c2 = fbdC2fxn()

  #log(4 * exp (-c1*t) / (((exp(-c1*t) * (1-c2)) + (1+c2))^2))

  f1 = log(4) + (-c1 * t)
  f2 = 2 * (log( (exp(-c1*t) * (1-c2)) + (1+c2) ))

  v = f1 - f2;

  return(v)
}

rangePrComplete<-function(gamma,bi,di){
  bi<-bi
  di<-di
  gamma<-gamma

  rp = log(lambda*gamma) + fbdQTildaFxnLog(bi) - fbdQTildaFxnLog(di)

  return(rp)

}

rangePr<-function(gamma,bi,di,oi){
  bi<-bi
  di<-di
  oi<-oi
  gamma<-gamma

  rp = log(lambda*gamma) + fbdQTildaFxnLog(oi) - fbdQTildaFxnLog(di) + fbdQfxnLog(bi) - fbdQfxnLog(oi)

  return(rp)

}

qt_heath<-function(t){
  t<-t

  c1=fbdC1fxn()
  c2=fbdC2fxn()

  v1 = 2 * ( 1 - c2^2 )
  v2 = exp(-c1 * t) * ((1 - c2)^2)
  v3 = exp(c1 * t) * ((1 + c2)^2)

  v = v1 + v2 + v3

  return(v)
}

# q=exp(fbdQfxnLog(told))/exp(fbdQfxnLog(tyoung)) - exp(fbdQTildaFxnLog(told))/exp(fbdQTildaFxnLog(tyoung))
# exp(fbdQTildaFxnLog(told))/exp(fbdQTildaFxnLog(tyoung)) - exp(-(told-tyoung) * (lambda+mu+psi) )
