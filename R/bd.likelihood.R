
######## fossil range estimates

#' Maximum likelihood estimation of birth and death rates for a given set of stratigraphic ranges
#'
#' This function maximises the likelihood of the function bd.probability.range using optim (Nelder-Mead option)
#'
#' @param b Initial value for birth rate
#' @param d Initial value for death rate
#' @param frs Dataframe of species ranges
#' @param lower.b Lower birth rate
#' @param upper.b Upper birth rate
#' @param lower.d Lower death rate
#' @param upper.d Upper death rate
#' @param crown If TRUE assume the process begins at the first speciation event and not the origin (default = FALSE)
#' @return Named list including maximum likelihood estimates of birth and death rate
#' @examples
#' # simulate tree & assume complete sampling
#' t = TreeSim::sim.bd.taxa(100,1,1,0.1)[[1]]
#' # add symmetric speciation events & generate fossil range dataframe
#' f = 0.5
#' ages <- FossilSim::mixed.ages(t, f, root.edge = TRUE)
#' # add anagenic speciation events
#' lambda.a = 0.1
#' frs <- FossilSim::anagenic.species(ages, lambda.a)
#' # estimate birth & death rates
#' birth = runif(1)
#' death = runif(1)
#' out = est.bd.ranges(birth, death, frs)
#' # ML birth rate
#' out$par[1]
#' # ML death rate
#' out$par[2]
#' @export
est.bd.ranges<-function(b,d,frs,lower.b=0.001,upper.b=10,lower.d=0.001,upper.d=10,crown=FALSE){
  b<-b # starting value
  d<-d # starting value
  frs<-frs
  lower.b<-lower.b
  lower.d<-lower.d
  upper.b<-upper.b
  upper.d<-upper.d
  crown<-crown
  method =  1 # 1: "Nelder-Mead" 2: "L-BFGS-B"

  if(method==2){
    est<-optim(c(b,d),bd.likelihood.est.ranges,frs=frs,crown=crown,control=list(fnscale=-1,maxit=1000), # fnscale=-1 tells the function to maximise
               method="L-BFGS-B", # this method allows bounds on parameters
               lower=c(lower.b,lower.d),
               upper=c(upper.b,upper.d)
    )
  }
  else{
    est<-optim(c(b,d),bd.likelihood.est.ranges,frs=frs,crown=crown,control=list(fnscale=-1,maxit=1000))
  }

  est$par = abs(est$par)

  return(est)
  #eof
}

bd.likelihood.est.ranges<-function(p,frs,crown=FALSE){
  b=abs(p[1])
  d=abs(p[2])
  frs<-frs
  crown<-crown

  lk=bd.probability.range(frs,b,d,crown=crown)
  return(lk)
  #eof
}

######## phylogenetic estimates

#' Maximum likelihood estimation of birth and death rates for a given phylogeny of extant taxa
#'
#' This function maximises the likelihood of the function bd.probability.extant using optim (Nelder-Mead option)
#'
#' @param b Initial value for birth rate
#' @param d Initial value for death rate. Must be < b
#' @param tree Phylo object of extant taxa (the function will remove any extinct taxa prior to calculating the likelihood)
#' @param lower.b Lower birth rate
#' @param upper.b Upper birth rate
#' @param lower.d Lower death rate
#' @param upper.d Upper death rate
#' @param crown If TRUE the process is conditioned on the crown (default = F)
#' @param rho Extant species sampling probability (default = 1)
#' @return Named list including maximum likelihood estimates of birth and death rate
#' @examples
#' # simulate tree
#' t = TreeSim::sim.bd.taxa(100,1,1,0.1)[[1]]
#' # estimate birth & death rates
#' birth = runif(1)
#' death = runif(1, 0, birth)
#' out = est.bd.extant(birth, death, t)
#' # ML birth rate
#' out$par[1]
#' # ML death rate
#' out$par[2]
#' @export
est.bd.extant<-function(b,d,tree,lower.b=0.001,upper.b=10,lower.d=0.001,upper.d=10,crown=FALSE,rho=1){
  b<-b # starting value
  d<-d # starting value
  tree<-tree
  lower.b<-lower.b
  lower.d<-lower.d
  upper.b<-upper.b
  upper.d<-upper.d
  crown<-crown
  rho<-rho
  method =  1 # 1: "Nelder-Mead" 2: "L-BFGS-B"

  if(method==2){
    est<-optim(c(b,d),bd.likelihood.est.extant,tree=tree,crown=crown,rho=rho,control=list(fnscale=-1,maxit=1000), # fnscale=-1 tells the function to maximise
               method="L-BFGS-B", # this method allows bounds on parameters
               lower=c(lower.b,lower.d),
               upper=c(upper.b,upper.d)
    )
  }
  else{
    est<-optim(c(b,d),bd.likelihood.est.extant,tree=tree,crown=crown,rho=rho,control=list(fnscale=-1,maxit=1000))
  }

  est$par = abs(est$par)

  return(est)
  #eof
}

bd.likelihood.est.extant<-function(p,tree,crown=FALSE,rho=1){
  b=abs(p[1])
  d=abs(p[2])
  tree<-tree
  crown<-crown
  rho<-rho

  if((b < d) || (b == d))
   return(-10^100)

  lk=bd.probability.extant(tree,b,d,crown=crown,rho=rho)
  return(lk)
  #eof
}

######## phylogenetic estimates constrained

est.bd.extant.constr<-function(nd,tree,crown=FALSE,rho=1){
  nd<-nd # fixed value
  tree<-tree
  crown<-crown
  rho<-rho

  # use optimize for one-dimensional optimization
  est<-optimize(bd.likelihood.est.extant.constr,tree=tree,nd=nd,crown=crown,rho=rho,c(0,1000),maximum = TRUE)

  #est$par = abs(est$par)

  return(est)
  # eof
}

bd.likelihood.est.extant.constr<-function(b,tree,nd,crown=FALSE,rho=1){
  b<-abs(b)
  tree<-tree
  nd<-nd
  crown<-crown
  rho<-rho

  if(nd >= b)
    return(-10^100)

  d = b - nd

  lk=bd.probability.extant(tree,b,d,crown=crown,rho=rho)
  return(lk)
  #eof
}

######## combined estimates

#' Maximum likelihood estimation of birth and death rates for a given phylogeny of extant taxa and set of stratigraphic ranges
#'
#' This function maximises the likelihood of the function bd.probability.extant and bd.probability.range using optim, either treating stratigraphic range and phylogenetic birth and death rates
#' independently (\code{constrained = FALSE}) or constraining the parameters under the birth-death chronospecies model (\code{constrained = TRUE}).
#'
#' @param tree Phylo object of extant taxa (the function will remove any extinct taxa prior to calculating the likelihood)
#' @param frs Dataframe of species ranges
#' @param b Initial value for the phylogenetic birth rate
#' @param d Initial value for the phylogenetic death rate. Must be < b
#' @param b.star Initial value for the stratigraphic range birth rate
#' @param d.star Initial value for the stratigraphic range death rate
#' @param nd Net diversification (b - d). Must be > 0
#' @param constrained If FALSE fossil and phylogenetic rates are treated as independent. If TRUE rates are constrained under the birth-death chronospecies model
#' @param constrained.p Number of parameters to be estimated under the constrained model. If p = 2 rates are equal for fossils and phylogenies (this is the equal rates model, a special case of the BD chronospecies model).
#' If p = 3 rates are not equal but b - d = b.star - d.star (this the compatible rates BD chronospecies model)
#' @param crown If TRUE the process is conditioned on the crown instead of the origin (default = F)
#' @param rho Extant species sampling probability (default = 1)
#' @return Named list including maximum likelihood estimates of diversification rate parameters.
#' Under the independent rates model (\code{constrained = FALSE}) the function returns estimates of birth, death, birth.star and death.star.
#' Under the compatible rates model (\code{constrained = TRUE, constrained.p = 3}) the function returns birth and birth.star estimates for the phylogeny and stratigraphic ranges, respectively, and a single diverification rate estimate.
#' Under the equal rates model (\code{constrained = TRUE, constrained.p = 2}) the function returns a single set of birth and death rates applicable to both the phylogeny and stratigraphic ranges.
#' @examples
#' # simulate tree & assume complete sampling
#' t = TreeSim::sim.bd.taxa(100,1,1,0.1)[[1]]
#' # add symmetric speciation events & generate fossil range dataframe
#' beta = 0.5
#' ages <- FossilSim::mixed.ages(t, beta, root.edge = TRUE)
#' # add anagenic speciation events
#' lambda.a = 0.1
#' frs <- FossilSim::anagenic.species(ages, lambda.a)
#'
#' # estimate birth & death rates
#'
#' # The equal rates birth-death chronospecies model
#' out = est.bd.combined(t, frs, constrained = TRUE, constrained.p = 2)
#' # ML birth rate
#' out$par[1]
#' # ML death rate
#' out$par[2]
#'
#' # The compatible rates birth-death chronospecies model
#' out = est.bd.combined(t, frs, constrained = TRUE, constrained.p = 3)
#' # ML birth rate
#' out$par[1]
#' # ML birth.star rate
#' out$par[2]
#' # ML death rate
#' out$par[1] - out$par[3]
#' # ML death.star rate
#' out$par[2] - out$par[3]
#'
#' # The indpendent rates birth-death model
#' out = est.bd.combined(t, frs, constrained = FALSE)
#' # ML birth rate
#' out$par[1]
#' # ML birth.star rate
#' out$par[2]
#' # ML death rate
#' out$par[3]
#' # ML death.star rate
#' out$par[4]
#' @export
est.bd.combined<-function(tree,frs,b=0.3,d=0.1,b.star=1,d.star=0.1,nd=0.1,constrained=FALSE,constrained.p=3,crown=FALSE,rho=1){
  b<-b # starting value
  d<-d # starting value
  b.star<-b.star # starting value
  d.star<-d.star # starting value
  nd<-nd # starting value
  crown<-crown
  rho<-rho

  if(!constrained)
    est<-optim(c(b,d,b.star,d.star),bd.likelihood.est.combined,tree=tree,frs=frs,crown=crown,rho=rho,control=list(fnscale=-1,maxit=1000)) # fnscale=-1 tells the function to maximise
  else if (constrained.p == 3)
    est<-optim(c(b,b.star,nd),bd.likelihood.est.combined.const,tree=tree,frs=frs,crown=crown,rho=rho,control=list(fnscale=-1,maxit=1000)) # fnscale=-1 tells the function to maximise
  else if (constrained.p == 2)
    est<-optim(c(b,d),bd.likelihood.est.combined.const.simple,tree=tree,frs=frs,crown=crown,rho=rho,control=list(fnscale=-1,maxit=1000)) # fnscale=-1 tells the function to maximise

  est$par = abs(est$par)

  return(est)
  # eof
}

# l, m, l*, m*
bd.likelihood.est.combined<-function(p,tree,frs,crown=FALSE,rho=1){
  b=abs(p[1])
  d=abs(p[2])
  b.star=abs(p[3])
  d.star=abs(p[4])
  tree<-tree
  frs<-frs
  crown<-crown
  rho<-rho

  if((b < d) || (b == d))
    return(-10^100)

  lk = bd.probability.range(frs,b.star,d.star,crown=crown) + bd.probability.extant(tree,b,d,crown=crown,rho=rho)

  return(lk)
}

# l, l*, d
bd.likelihood.est.combined.const<-function(p,tree,frs,crown=FALSE,rho=1){
  b=abs(p[1])
  b.star=abs(p[2])
  nd=abs(p[3])
  tree<-tree
  frs<-frs
  crown<-crown
  rho<-rho

  if(b > b.star)
    return(-10^100)

  if(nd >= b || nd >= b.star)
    return(-10^100)

  d = b - nd
  d.star = b.star - nd

  lk = bd.probability.range(frs,b.star,d.star,crown=crown) + bd.probability.extant(tree,b,d,crown=crown,rho=rho)

  return(lk)
}

# l, m
bd.likelihood.est.combined.const.simple<-function(p,tree,frs,crown=FALSE,rho=1){
  b=abs(p[1])
  d=abs(p[2])
  tree<-tree
  frs<-frs
  crown<-crown
  rho<-rho

  if(b < d)
    return(-10^100)

  b.star = b
  d.star = d

  lk = bd.probability.range(frs,b.star,d.star,crown=crown) + bd.probability.extant(tree,b,d,crown=crown,rho=rho)

  return(lk)
}

