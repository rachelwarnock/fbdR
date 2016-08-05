
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
fbd.probability<-function(frs,b,d,s,k,rho=1){
  frs<-frs
  lambda<<-b
  mu<<-d
  psi<<-s
  rho<<-rho
  numFossils<-k

  ot = max(frs$bi) # define the origin
  extinctLineages = length(frs$extant) - sum(frs$extant)

  pr = numFossils*log(psi)
  pr = pr + extinctLineages*log(mu)
  pr = pr - log(lambda * (1 -fbdPfxn(ot) ) )

  #   for (fr in 1:length(frs$bi)) {
  #     gamma=frs$gamma[fr]
  #     bi=frs$bi[fr]
  #     di=frs$di[fr]
  #
  #     rangePr = log(lambda*gamma) + fbdQTildaFxnLog(bi) - fbdQTildaFxnLog(di)
  #
  #     pr = pr + rangePr
  #   }

  # this makes it a bit faster
  rp=sum(unlist(lapply(1:length(frs$bi),function(x){ rangePr(frs$gamma[x],frs$bi[x],frs$di[x]) })))

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

  p = 1 + ( (-(lambda-mu-psi) +  (c1 * ( ( exp(-c1*t)*(1-c2)-(1+c2) ) / ( exp(-c1*t)*(1-c2)+(1+c2) ) ) ) ) / (2*lambda) )

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

rangePr<-function(gamma,bi,di){
  bi<-bi
  di<-di
  gamma<-gamma

  rp = log(lambda*gamma) + fbdQTildaFxnLog(bi) - fbdQTildaFxnLog(di)

  return(rp)

}


