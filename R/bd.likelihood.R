
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

########

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
    return(-10000000)

  lk=bd.probability.extant(tree,b,d)
  return(lk)
  #eof
}

# constrOptim(c(1,0.1),bd.likelihood.est.extant,NULL,ui=rbind(c(-1,0),c(0,-1),c(1,0),c(0,1),c(1,-1)),
#             ci=c(-10,-10,0.001,0.001,0.001),
#             tree=tmp[[1]],control=list(fnscale=-1,maxit=10000)
# )
