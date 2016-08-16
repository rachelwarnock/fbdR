
#' @export
fbd.likelihood.est<-function(b,d,s,k,frs,est.b=T,est.d=T,est.s=F,lower.b=0.001,lower.d=0.001,lower.s=0.001,upper.b=10,upper.d=10,upper.s=100,complete=F){
  b<-b # fixed or starting value
  d<-d # fixed or starting value
  s<-s # fixed or starting value
  k<-k
  frs<-frs
  est.b<-est.b
  est.d<-est.d
  est.s<-est.s
  lower.b<-lower.b
  lower.d<-lower.d
  lower.s<-lower.s
  upper.b<-upper.b
  upper.d<-upper.d
  upper.s<-upper.s
  complete<-complete

  # estimate b, d and s
  if(est.b & est.d & est.s) {

    est<-est.bds(b,d,s,k,frs,lower.b,upper.b,lower.d,upper.d,lower.s,upper.s,complete)
    return(list(lambda=est$par[1],mu=est$par[2],psi=est$par[3],likelihood=est$value))

  }
  # estimate b and d only, fix s
  else if (est.b & est.d){

    est<-est.bd(b,d,s,k,frs,lower.b,upper.b,lower.d,upper.d,complete)
    return(list(lambda=est$par[1],mu=est$par[2],likelihood=est$value))

  }
  # estimate b, d or s only
  else if (est.b){

    est<-est.d(b,d,s,k,frs,lower.b,upper.b)
    return(list(lambda=est$maximum,likeilhood=est$objective))

  }
  else if (est.d){

    est<-est.d(b,d,s,k,frs,lower.d,upper.d)
    return(list(mu=est$maximum,likeilhood=est$objective))

  }
  else if (est.s){

    est<-est.s(b,d,s,k,frs,lower.d,upper.d)
    return(list(psi=est$maximum,likeilhood=est$objective))

  }
  else {
    return("You didn't ask me to estimate anything!")
  }

  #eof
}

# function to co-estimate lambda, mu, psi

est.bds<-function(b,d,s,k,frs,lower.b,upper.b,lower.d,upper.d,lower.s,upper.s,complete){
  b<-b # starting value
  d<-d # starting value
  s<-s # starting value
  k<-k
  frs<-frs
  lower.b<-lower.b
  upper.b<-upper.b
  lower.d<-lower.d
  upper.d<-upper.d
  lower.s<-lower.s
  upper.s<-upper.s
  complete<-complete

  est<-optim(c(b,d,s),fbd.likelihood.est.bds,k=k,frs=frs,complete=complete,control=list(fnscale=-1,maxit=500),
             method="L-BFGS-B", # this method allows bounds on parameters
             lower=c(lower.b,lower.d,lower.s),
             upper=c(upper.b,upper.d,upper.s)
  )
  return(est)
  #eof
}

fbd.likelihood.est.bds<-function(p,k,frs,complete){
  b=p[1]
  d=p[2]
  s=p[3]
  k<-k
  frs<-frs
  complete<-complete

  lk=fbd.probability(frs,b,d,s,k,rho=1,complete=complete)
  return(lk)
  #eof
}

# functions to co-estimate lambda & mu

est.bd<-function(b,d,s,k,frs,lower.b,upper.b,lower.d,upper.d){
  b<-b # starting value
  d<-d # starting value
  s<-s # fixed
  k<-k
  frs<-frs
  lower.b<-lower.b
  lower.d<-lower.d
  upper.b<-upper.b
  upper.d<-upper.d
  complete<-complete

  est<-optim(c(b,d),fbd.likelihood.est.bd,s=s,k=k,frs=frs,complete=complete,control=list(fnscale=-1,maxit=500), # fnscale=-1 tells the function to maximise
         method="L-BFGS-B", # this method allows bounds on parameters
         lower=c(lower.b,lower.d),
         upper=c(upper.b,upper.d)
        )
  return(est)
  #eof
}

fbd.likelihood.est.bd<-function(p,s,k,frs,complete){
  b=p[1]
  d=p[2]
  s<-s
  k<-k
  frs<-frs
  complete<-complete

  lk=fbd.probability(frs,b,d,s,k,rho=1,complete=complete)
  return(lk)
  #eof
}

# function to estimate lambda

est.b<-function(b,d,s,k,frs,lower.b,upper.b){
  b<-b
  d<-d
  s<-s
  k<-k
  frs<-frs
  lower.b<-lower.b
  upper.b<-upper.b

  est<-optimize(function(b) {fbd.probability(frs=frs,b,d=d,s=s,k=k)},
                interval = c(lower.b, upper.b),
                maximum = TRUE)

  return(est)
  #eof
}

# function to estimate mu

est.d<-function(b,d,s,k,frs,lower.d,upper.d){
  b<-b
  d<-d
  s<-s
  k<-k
  frs<-frs
  lower.d<-lower.d
  upper.d<-upper.d

  est<-optimize(function(d) {fbd.probability(frs=frs,b=b,d,s=s,k=k)},
             interval = c(lower.d, upper.d),
             maximum = TRUE)

  return(est)
  #eof
}

# function to estimate psi

est.s<-function(b,d,s,k,frs,lower.s,upper.s){
  b<-b
  d<-d
  s<-s
  k<-k
  frs<-frs
  lower.s<-lower.s
  upper.s<-upper.s

  est<-optimize(function(s) {fbd.probability(frs=frs,b=b,d=d,s,k=k)},
                interval = c(lower.s,upper.s),
                maximum = TRUE)

  return(est)
  #eof
}

# zero indicates successful covergence
# est$convergence == 0


