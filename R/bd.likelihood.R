
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

  est<-optim(c(b,d),bd.likelihood.est.ranges,frs=frs,control=list(fnscale=-1,maxit=500), # fnscale=-1 tells the function to maximise
             method="L-BFGS-B", # this method allows bounds on parameters
             lower=c(lower.b,lower.d),
             upper=c(upper.b,upper.d)
  )
  return(est)
  #eof
}

bd.likelihood.est.ranges<-function(p,frs){
  b=p[1]
  d=p[2]
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

  est<-optim(c(b,d),bd.likelihood.est.extant,tree=tree,control=list(fnscale=-1,maxit=10000), # fnscale=-1 tells the function to maximise
             method="L-BFGS-B", # this method allows bounds on parameters
             lower=c(lower.b,lower.d),
             upper=c(upper.b,upper.d)
  )
  return(est)
  #eof
}

bd.likelihood.est.extant<-function(p,tree){
  b=p[1]
  d=p[2]
  tree<-tree

  if((b < d) || (b == d))
    return(-10^100)

  lk=bd.probability.extant(tree,b,d)
  return(lk)
  #eof
}

######## combined estimates
#' @export
est.bd.combined<-function(b,d,b.star,d.star,tree,frs,
                         lower.b=0.001,upper.b=10,lower.d=0.001,upper.d=10,
                         lower.b.star=0.001,upper.b.star=10,lower.d.star=0.001,upper.d.star=10,
                         constrained=FALSE){
  b<-b # starting value
  d<-d # starting value
  b.star<-b.star # starting value
  d.star<-d.star # starting value

  est<-optim(c(b,d,b.star,d.star),bd.likelihood.est.combined,tree=tree,frs=frs,control=list(fnscale=-1,maxit=20000), # fnscale=-1 tells the function to maximise
             method="L-BFGS-B", # this method allows bounds on parameters
             lower=c(lower.b,lower.d,lower.b.star,lower.d.star),
             upper=c(upper.b,upper.d,upper.b.star,upper.d.star)
  )

  return(est)
  # eof
}

bd.likelihood.est.combined<-function(p,tree,frs){
  b=p[1]
  d=p[2]
  b.star=p[3]
  d.star=p[4]
  tree<-tree
  frs<-frs

  if((b < d) || (b == d))
    return(-10^100)

  lk = bd.probability.range(frs,b.star,d.star) + bd.probability.extant(tree,b,d)

  return(lk)
}
