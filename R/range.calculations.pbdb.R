#' Generate presence/absence matrix from simulated data
#'
#' @param fossils Fossil dataframe.
#' @param intervals Intervals dataframe.
#'
#' @return Dataframe recording presences & absences during different time intervals.
#'
#' @examples
#' set.seed(123)
#' # simulate tree
#' t = TreeSim::sim.bd.taxa(100,1,0.01,0.005)[[1]]
#'
#' # simulate fossils
#' f = FossilSim::sim.fossils.poisson(t, 1)
#'
#' # asymmetric species assignment
#' f = FossilSim::asymmetric.fossil.mapping(t, f)
#'
#' # generate presence.absence matrix
#' presence.absence.matrix(f, stages)
#'
#' @export
presence.absence.matrix<-function(fossils, intervals, add.extant = TRUE){

  df = data.frame(matrix(ncol = length(intervals$name) + 1, nrow = 0))
  colnames(df) = c("sp", as.character(intervals$name))

  for(i in unique(fossils$sp)){

    h = fossils$h[which(fossils$sp == i)]

    s = c(i)

    for(j in 1:length(intervals$name)){
      if(add.extant){
        if(any(h >= intervals$end[j] & h < intervals$start[j]))
          s = c(s, 1)
        else s = c(s, 0)
      } else{
        if(any(h > intervals$end[j] & h < intervals$start[j]))
          s = c(s, 1)
        else s = c(s, 0)
      }
    }
    df[nrow(df) + 1,] = s
  }
  return(df)
  #eof
}

##

#' Work out first and last appearance intervals
#'
#' @param presence.absence Dataframe of presence/absence data. Each row represents a taxon and each column represents a time interval from oldest to youngest.
#'
#' @return Dataframe containing first and last appearance intervals. Intervals are numbered oldest (= 1) to youngest (> 1).
#'
#' @examples
#' set.seed(123)
#' # simulate tree
#' t = TreeSim::sim.bd.taxa(100,1,0.01,0.005)[[1]]
#'
#' # simulate fossils
#' f = FossilSim::sim.fossils.poisson(t, 1)
#'
#' # asymmetric species assignment
#' f = FossilSim::asymmetric.fossil.mapping(t, f)
#'
#' # generate presence.absence matrix
#' pa = presence.absence.matrix(f, stages)
#'
#' # work out first and last appearance intervals
#' fa = first.last.appearances.pa(pa)
#'
#' @export
first.last.appearances.pa<-function(presence.absence){

  n = length(presence.absence[,1])

  fa.la<-data.frame(sp = numeric(), fa = numeric(), la = numeric()) # dataframe for FAs & LAs

  for(i in 1:n){

    sp = presence.absence[1][i,]

    if(any(presence.absence[-1][i,] == 1)){
      bins = which(presence.absence[-1][i,] == 1)
      FA = bins[1]
      LA = bins[length(bins)]
      fa.la <- rbind(fa.la, data.frame(sp = sp,FA = FA, LA = LA))
    }
  }
  return(fa.la)
  #eof
}

## Categorise interval types

#' Categorise per interval taxon types using boundary crosser categories based on presence/absence data
#'
#' @details
#' Taxa types are detailed in Foote (2000). \cr
#' Nbt - bottom and top boundary crossers \cr
#' Nbl - bottom boundary crossers only (e.g. species goes extinct) \cr
#' NFt - top boundary crossers only (e.g. species originates) \cr
#' NFl - species originates and becomes extinct during the same interval \cr
#'
#' @param presence.absence Dataframe of presence/absence data. Each row represents a taxon and each column represents a time interval from oldest to youngest.
#' @param intervals Dataframe of geological intervals.
#' @return dataframe of per interval taxon types.
#'
#' @references
#' Foote, M. 2000. Origination and extinction components of taxonomic diversity: General problems. Paleobiology 26: 74-102.
#'
#' @examples
#' set.seed(123)
#' # simulate tree
#' t = TreeSim::sim.bd.taxa(100,1,0.01,0.005)[[1]]
#'
#' # simulate fossils
#' f = FossilSim::sim.fossils.poisson(t, 1)
#'
#' # asymmetric species assignment
#' f = FossilSim::asymmetric.fossil.mapping(t, f)
#'
#' # generate presence.absence matrix
#' pa = presence.absence.matrix(f, stages)
#'
#' # categorise interval types
#' interval.types.bc.pa(pa, stages)
#'
#' @export
interval.types.bc.pa<-function(presence.absence, intervals){

  FAs<-first.last.appearances.pa(presence.absence)

  taxa.types<-data.frame(int = numeric(),
                         NFl = numeric(),
                         NFt = numeric(),
                         Nbl = numeric(),
                         Nbt = numeric()) # dataframe for taxon types

  # during each interval
  for(h in 1:length(intervals$name)){

    NFls = c()
    NFts = c()
    Nbls = c()
    Nbts = c()

    # for each lineage
    for(i in 1:length(FAs$sp)) {

      fa = FAs[,2][i] #FA
      la = FAs[,3][i] #LA
      id = FAs[,1][i] #sp

      # define NFls (species originate & go extinct)
      # or NFts (top boundary crossers - species orginates)
      # if the first appearance is equivalent to the current horizon:
      if(fa == h) {
        # and the first appearance is equivalent to the last appearance
        # species is categorised as NFL
        if(fa == la) {
          NFls = c(NFls, FAs[,1][i])
          #useful.un = c(useful.un, id)
        }
        # else species is categorised as NFt
        else {
          NFts = c(NFts, FAs[,1][i])
          #useful.bc = c(useful.bc, id)
        }
      }
      # else define Nbl (bottom boundary crossers - species goes extinct)
      else if (la == h) {
        Nbls = c(Nbls, FAs[,1][i])
        #useful.bc = c(useful.bc, id)
      }
      # else define Nbt (bottom and top boundary crossers)
      else if ( (fa < h) & (la > h) ) {
        Nbts = c(Nbts, FAs[,1][i])
        #useful.bc = c(useful.bc, id)
      }
    }
    taxa.types<-rbind(taxa.types, data.frame(int = h,
                                             NFl = length(NFls),
                                             NFt = length(NFts),
                                             Nbl = length(Nbls),
                                             Nbt = length(Nbts)))
  }
  return(taxa.types)
  #eof
}

## Calculate speciation and extinction rates

#' Calculate speciation and extinction rates using the boundary crosser approach
#'
#' @param presence.absence Dataframe of presence/absence data. Each row represents a taxon and each column represents a time interval from oldest to youngest.
#' @param intervals Dataframe of geological intervals.
#' @param continuous If TRUE calculate continuous rates (i.e. account for interval length)
#'
#' @return Dataframe containing the following columns:
#' \itemize{
#' \item \code{int} interval number. Intervals are numbered from oldest (=1) to youngest (>1)
#' \item \code{int.name} Interval name
#' \item \code{NFl} singletons
#' \item \code{NFt} top boundary crossers
#' \item \code{Nbl} bottom boundary crossers
#' \item \code{Nbt} range through taxa
#' \item \code{p} speciation rate
#' \item \code{q} extinction rate
#' }
#' Note this approach does not return rates for the first interval.
#'
#' @references
#' Foote, M. 2000. Origination and extinction components of taxonomic diversity: General problems. Paleobiology 26: 74-102.
#'
#' @seealso \code{\link{interval.types.bc.pa}}
#'
#' @examples
#'
#' # option 1: generate simulated data
#'
#' set.seed(111)
#'
#' # simulate tree
#' birth = 0.02
#' death = 0.01
#' tips = 400
#' t = TreeSim::sim.bd.taxa(tips, 1, birth, death)[[1]]
#'
#' # simulate fossils
#' f = FossilSim::sim.fossils.poisson(t, 1)
#'
#' # add extant occurrences
#' f = FossilSim::add.extant.occ(t, f)
#'
#' # asymmetric species assignment (this may take a while esp. for large trees)
#' f = FossilSim::asymmetric.fossil.mapping(t, f)
#'
#' # generate presence.absence matrix
#' pa = presence.absence.matrix(f, stages)
#'
#' # option 2: load precooked presence absence matrix
#' pa = presence.absence.precooked
#'
#' # calculate speciation & extinction rates
#' out = boundary.crosser.rates.pa(pa, stages)
#'
#' # rates
#' mean(out$p, na.rm = TRUE)
#' mean(out$q, na.rm = TRUE)
#'
#' # plot the output
#' snum = length(stages$name)
#' plot(1:snum, out$p, pch = 16, bty = "n", xlab = "stage #", ylab = "speciation")
#' lines(1:snum, rep(birth, snum), lty = 2, lwd = "2", col = 2)
#'
#' plot(1:snum, out$q, pch = 16, bty = "n", xlab = "stage #", ylab = "extinction")
#' lines(1:snum, rep(death, snum), lty = 2, lwd = "2", col = "deepskyblue")
#'
#' @export
boundary.crosser.rates.pa<-function(presence.absence, intervals, continuous = TRUE) {
  return.global = FALSE

  taxa.types = interval.types.bc.pa(presence.absence, intervals)

  # calculate delta t
  s1 = intervals$start - intervals$end

  o.e.rates<-data.frame(int = numeric(),
                        int.name = character(),
                        NFl = numeric(),
                        NFt = numeric(),
                        Nbl = numeric(),
                        Nbt = numeric(),
                        p = numeric(),
                        q = numeric()) # dataframe for speciation and extinction rates

  k.total = 0 # range through taxa
  n.total = 0 # range through taxa + specation events

  e.total = 0 # range through taxa
  Ntot.total = 0 # range through taxa + extinction events

  for(h in 1:length(taxa.types[,1])){

    NFl = taxa.types[,2][h]
    NFt = taxa.types[,3][h]
    Nbl = taxa.types[,4][h]
    Nbt = taxa.types[,5][h]

    k.total = k.total + Nbt
    n.total = n.total + Nbt + NFt

    e.total = e.total + Nbt
    Ntot.total = Ntot.total + Nbt + Nbl

    if(Nbt > 0){
      if(!continuous){
        # calculate orgination rates
        p.hat = -log(Nbt/(Nbt+NFt))
        p.hat = round(p.hat,3)
        # calculate extinction rates
        q.hat = -log(Nbt/(Nbt+Nbl))
        q.hat = round(q.hat,3)
      } else {
        # calculate orgination rates
        p.hat = -log(Nbt/(Nbt+NFt)) / s1[h]
        p.hat = round(p.hat,3)
        # calculate extinction rates
        q.hat = -log(Nbt/(Nbt+Nbl)) / s1[h]
        q.hat = round(q.hat,3)
      }
    } else{
      p.hat = NA
      q.hat = NA
    }
    o.e.rates<-rbind(o.e.rates, data.frame(int = h,
                                           int.name = intervals$name[h],
                                           NFl = NFl,
                                           NFt = NFt,
                                           Nbl = Nbl,
                                           Nbt = Nbt,
                                           p = p.hat, q = q.hat))
  }

  if(k.total > 0){
    if(!continuous){
      p.total = (-log(k.total/n.total))
      q.total = (-log(e.total/Ntot.total))
    }
    else{
      p.total = (-log(k.total/n.total)) / (intervals$start[1] - intervals$end[length(intervals$name)])
      q.total = (-log(e.total/Ntot.total)) / (intervals$start[1] - intervals$end[length(intervals$name)])
    }
  }
  else{
    p.total = NA
    q.total = NA
  }

  if(return.global)
    return(list(speciation = p.total, extinction = q.total))
  else
    return(per.interval.rates = o.e.rates)
  #eof
}

#TODO do a test with rate variation across intervals
#TODO set up a massive tree rep on euler
#TODO message Alex

#END

#TODO add output = proportion of useful taxa
#TODO make interval order more flexible
