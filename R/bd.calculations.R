

#' Birth-death probability (Keiding, 1975)
#'
#' This is the underlying birth-death process used in Silvestro et al. 2014, eq. 9
#'
#' @param frs Dataframe of species ranges
#' @param b Rate of speciation
#' @param d Rate of extinction
#' @return Log likelihood
#' @export
bd.probability.range<-function(frs,b,d,crown=FALSE){
  frs<-frs
  b<-b
  d<-d
  crown<-crown # condition on the crown, instead of the origin

  B = 0 # total number of speciation events # int
  D = 0 # total number of extinction events # int
  S = 0 # total lineage duration

  if(crown)
    B = length(frs$sp)
  else
    B = length(frs$sp)-1 # note the origin is not a speciation event

  D = length(which(frs$end!=0))
  S = sum(frs$start-frs$end)

  #likelihood = (lambda^B) * (mu^D) * exp(-(lambda+mu)*S)
  ll = (B * log(b)) + (D * log(d)) + (-(b + d) * S)

  return(ll)

}

#' Birth-death probability (Stadler, 2010)
#'
#' Probability of extant species phylogeny conditioned on the origin (Stadler, 2010, eq. 2)
#' or the crown (crown = T) (Stadler, 2010, eq. 5)
#'
#' @param tree Phylo object of extant taxa only
#' @param b Rate of speciation (branching)
#' @param d Rate of extinction (branch termination)
#' @param rho Extant species sampling
#' @param crown Boolean stating whether the process is conditioned on the crown (default = F)
#' @return Log likelihood
#' @export
bd.probability.extant<-function(tree,b,d,rho=1,mpfr=FALSE,crown=FALSE){
  tree<-tree
  bits = 128
  if(mpfr){
    lambda<<-mpfr(b, bits)
    mu<<-mpfr(d, bits)
    rho<<-mpfr(rho, bits)
  }
  else{
    lambda<<-b
    mu<<-d
    rho<<-rho
  }
  rho<-rho
  crown<-crown

  tree<-geiger::drop.extinct(tree)

  # calculate node ages
  node.ages = TreeSim::getx(tree)

  # equation 5
  if(crown){

    # take care of the crown
    c = max(node.ages)
    ll = 2 * ( bdP1Log(c) - log(1 - bdP0(c)) )

    # eliminate the root node
    node.ages = node.ages[-which(node.ages==c)]

  }
  else { # equation 2

    # take care of the origin
    origin = max(node.ages)+tree$root.edge
    ll = bdP1Log(origin) - log(1 - bdP0(origin))

  }

  for(i in node.ages){
    ll = ll + log(lambda) + bdP1Log(i)
  }

  return(ll)
}

# bdP1 function Stadler 2010, 322
bdP1<-function(t){
  t<-t

  p1 = rho * (lambda - mu)^2 * exp(-(lambda-mu)*t)
  p2 = p1 / (((rho * lambda) + (((lambda * (1 - rho)) - mu) * exp(-(lambda-mu)*t)))^2)

  return(p2)
}

bdP1Log<-function(t){
  t<-t

  p1 = log(rho) + (2 * log(lambda - mu) ) + (-(lambda-mu)*t)

  p2 = p1 - (2 * log((rho * lambda) + (((lambda * (1 - rho)) - mu) * exp(-(lambda-mu)*t))) )

  return(p2)
}

# bdP0 function Stadler 2010, 322 or Phat in Heath et al. 2014
bdP0<-function(t){
  t<-t

  p = 1 - ( (rho * (lambda-mu)) / ((rho * lambda) + ((lambda*(1-rho)-mu)*exp(-(lambda-mu)*t)) ) )

  return(p)
}

# this doesn't work if 1-(b-d) < 0
bdP0Log<-function(t){
  t<-t

  p = log(1 - (rho * (lambda-mu)) ) - log( ((rho * lambda) + ((lambda*(1-rho)-mu)*exp(-(lambda-mu)*t)) ) )

  return(p)
}

