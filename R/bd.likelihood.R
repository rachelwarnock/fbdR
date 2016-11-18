
######## fossil range estimates

#' @export
est.bd.ranges<-function(b,d,frs,lower.b=0.001,upper.b=10,lower.d=0.001,upper.d=10){
  b<-b # starting value
  d<-d # starting value
  frs<-frs
  lower.b<-lower.b
  lower.d<-lower.d
  upper.b<-upper.b
  upper.d<-upper.d
  method =  1 # 1: "Nelder-Mead" 2: "L-BFGS-B"

  if(method==2){
    est<-optim(c(b,d),bd.likelihood.est.ranges,frs=frs,control=list(fnscale=-1,maxit=1000), # fnscale=-1 tells the function to maximise
               method="L-BFGS-B", # this method allows bounds on parameters
               lower=c(lower.b,lower.d),
               upper=c(upper.b,upper.d)
    )
  }
  else{
    est<-optim(c(b,d),bd.likelihood.est.ranges,frs=frs,control=list(fnscale=-1,maxit=1000))
  }

  est$par = abs(est$par)

  return(est)
  #eof
}

bd.likelihood.est.ranges<-function(p,frs){
  b=abs(p[1])
  d=abs(p[2])
  frs<-frs

  lk=bd.probability.range(frs,b,d)
  return(lk)
  #eof
}

######## phylogenetic estimates

#' @export
est.bd.extant<-function(b,d,tree,lower.b=0.001,upper.b=10,lower.d=0.001,upper.d=10){
  b<-b # starting value
  d<-d # starting value
  tree<-tree
  lower.b<-lower.b
  lower.d<-lower.d
  upper.b<-upper.b
  upper.d<-upper.d
  method =  1 # 1: "Nelder-Mead" 2: "L-BFGS-B"

  if(method==2){
    est<-optim(c(b,d),bd.likelihood.est.extant,tree=tree,control=list(fnscale=-1,maxit=1000), # fnscale=-1 tells the function to maximise
             method="L-BFGS-B", # this method allows bounds on parameters
             lower=c(lower.b,lower.d),
             upper=c(upper.b,upper.d)
             )
  }
  else{
    est<-optim(c(b,d),bd.likelihood.est.extant,tree=tree,control=list(fnscale=-1,maxit=1000))
  }

  est$par = abs(est$par)

  return(est)
  #eof
}

bd.likelihood.est.extant<-function(p,tree){
  b=abs(p[1])
  d=abs(p[2])
  tree<-tree

  if((b < d) || (b == d))
   return(-10^100)

  lk=bd.probability.extant(tree,b,d)
  return(lk)
  #eof
}

######## combined estimates

#' @export
est.bd.combined<-function(tree,frs,b=0.3,d=0.1,b.star=1,d.star=0.1,nd=0.1,constrained=FALSE){
  b<-b # starting value
  d<-d # starting value
  b.star<-b.star # starting value
  d.star<-d.star # starting value
  nd<-nd # starting value

  if(!constrained)
    est<-optim(c(b,d,b.star,d.star),bd.likelihood.est.combined,tree=tree,frs=frs,control=list(fnscale=-1,maxit=1000)) # fnscale=-1 tells the function to maximise
  else
    est<-optim(c(b,b.star,nd),bd.likelihood.est.combined.const,tree=tree,frs=frs,control=list(fnscale=-1,maxit=1000)) # fnscale=-1 tells the function to maximise

  est$par = abs(est$par)

  return(est)
  # eof
}

bd.likelihood.est.combined<-function(p,tree,frs){
  b=abs(p[1])
  d=abs(p[2])
  b.star=abs(p[3])
  d.star=abs(p[4])
  tree<-tree
  frs<-frs

  if((b < d) || (b == d))
    return(-10^100)

  lk = bd.probability.range(frs,b.star,d.star) + bd.probability.extant(tree,b,d)

  return(lk)
}

bd.likelihood.est.combined.const<-function(p,tree,frs){
  b=abs(p[1])
  b.star=abs(p[2])
  nd=abs(p[3])
  tree<-tree
  frs<-frs

  if(b > b.star)
    return(-10^100)

  if(nd >= b || nd >= b.star)
    return(-10^100)

  d = b - nd
  d.star = b.star - nd

  lk = bd.probability.range(frs,b.star,d.star) + bd.probability.extant(tree,b,d)

  return(lk)
}


