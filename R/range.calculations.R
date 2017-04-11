#' Calculate out first and last appearances
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @return Function returns the first and last appearances for each lineage, represented by the max secure age of the sampling horizon.
#' @export
first.last.appearances<-function(fossils) {
  fossils<-fossils

  fa.la<-data.frame(node=numeric(),fa=numeric(),la=numeric()) # dataframe for FAs

  for (i in unique(fossils$sp)) {
    a = i

    # identify fossils
    f.row=which(fossils[,2]==a)
    occs=c()
    if(length(f.row) > 0){
      occs=c(occs,fossils$h[f.row])
    }

    # find the oldest and youngest occurrences
    if(length(occs) > 0) {
      max.fa=max(occs)
      max.la=min(occs)
      fa.la<-rbind(fa.la,data.frame(node=a,fa=max.fa,la=max.la))
    }
  }
  return(fa.la)
  #eof
}

## Categorize interval types

#' Categorize per interval taxon types using the boundary crosser approach
#'
#' @details
#' Taxa types are detailed in Foote (2000). \cr
#' Nbt - bottom and top boundary crossers \cr
#' Nbl - bottom boundary crossers only (e.g. species goes extinct) \cr
#' NFt - top boundary crossers only (e.g. species originates) \cr
#' NFl - species originates and becomes extinct during the same interval \cr
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @param basin.age Maximum age of the oldest stratigraphic interval
#' @param strata Number of stratigraphic horizons
#' @return dataframe of per interval taxon types.
#' @export
interval.types.bc<-function(fossils,basin.age,strata) {
  fossils<-fossils
  basin.age<-basin.age
  strata<-strata

  FAs<-first.last.appearances(fossils)

  s1=basin.age/strata # max age of youngest horizon/horizon length
  horizons<-seq(s1, basin.age, length=strata)

  taxa.types<-data.frame(horizons=numeric(),NFl=numeric(),NFt=numeric(),Nbl=numeric(),Nbt=numeric()) # dataframe for taxon types

  # during each horizon
  for(h in horizons){ # loop 1

    NFls=c()
    NFts=c()
    Nbls=c()
    Nbts=c()

    # for each lineage
    for(i in 1:length(FAs[,1])) {
      fa=FAs[,2][i]
      la=FAs[,3][i]

      # define NFls (species originate & go extinct)
      # or NFts (top boundary crossers - species orginates)
      # if the first appearance is equivalent to the current horizon:
      if(fa==h) {
        # and the first appearance is equivalent to the last appearance
        # species is categorized as NFL
        if(fa==la) {
          NFls=c(NFls,FAs[,1][i])
        }
        # else species is categorized as NFt
        else {
          NFts=c(NFts,FAs[,1][i])
        }
      }
      # else define Nbl (bottom boundary crossers - species goes extinct)
      else if (la==h) {
        Nbls=c(Nbls,FAs[,1][i])
      }
      # else define Nbt (bottom and top boundary crossers)
      else if ( (fa > h) & (la < h) ) {
        Nbts=c(Nbts,FAs[,1][i])
      }
    }
    taxa.types<-rbind(taxa.types,data.frame(horizons=h,NFl=length(NFls),NFt=length(NFts),Nbl=length(Nbls),Nbt=length(Nbts)))
  }
  return(taxa.types)
  #eof
}

#' Categorize per interval taxon types using the three-timer approach (Alroy, 2008)
#'
#' @details
#' Taxa types are detailed in Alroy (2000) and Alroy (2014) \cr
#' 1. taxa sampled at all in a focal bin (Ns)
#' 2. taxa sampled in a bin but not immediately before or after (one-timers, or 1t) \cr
#' 3a. taxa sampled immediately before and within the ith bin (two-timers, or 2ti) * referred to as two_t_a \cr
#' 3b. or within and immediately after the ith bin (2ti+1) * referred to as two_t_b \cr
#' 4. taxa sampled in three consecutive bins (three-timers, or 3t) \cr
#' 5. taxa sampled before and after but not within a bin (part-timers, or Pt) \cr
#' n.b. all three timers are also two timers \cr
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @param basin.age Maximum age of the oldest stratigraphic interval
#' @param strata Number of stratigraphic horizons
#' @return dataframe of per interval taxon types.
#' @export
interval.types.3t<-function(fossils,basin.age,strata) {
  fossils<-fossils
  basin.age<-basin.age
  strata<-strata

  lineages<-unique(fossils$sp)

  s1=basin.age/strata # max age of youngest horizon/horizon length
  horizons<-seq(s1, basin.age, length=strata)

  taxa.types<-data.frame(horizons=numeric(),Ns=numeric(),one_t=numeric(),two_t_a=numeric(),two_t_b=numeric(),three_t=numeric(),Pt=numeric()) # data fram for taxon types

  # during each horizon
  for(h in horizons){ # 1
    ha=h+s1 # horizon before
    hb=h-s1 # horizon after

    Ns=c()  # total number of taxa sampled in the bin
    one_t=c() # one timers
    two_t_a=c() # two timers = extinction event
    two_t_b=c() # two timers = speciation event
    three_t=c() # three timers
    Pt=c() # part timers

    # for each lineage
    for(l in lineages){ # 2

      # bins in which lineage l occurrs
      i=fossils[which(fossils$sp==l),]$h

      # do I exist in the bin?
      if(any(i==as.character(h))){

        Ns=c(Ns,l)

        # do I exist in the previous bin?
        if(any(i==as.character(ha))){

          # do I also exist in the next bin?
          if(any(i==as.character(hb))){
            three_t=c(three_t,l)
          }
          else {
            two_t_a=c(two_t_a,l)
          }
        }
        # do I exist in the next bin but not the one before?
        else if (any(i==as.character(hb))){
          two_t_b=c(two_t_b,l)
        }

        # otherwise record taxon as one timer
        else {
          one_t=c(one_t,l)
        }

      }

      #else, do I exist in the bin before and after?
      else if( (any(i==as.character(ha))) & (any(i==as.character(hb))) ){
        Pt=c(Pt,l)
      }

    } # 2

    two_t_a=c(two_t_a,three_t)
    two_t_b=c(two_t_b,three_t)

    taxa.types<-rbind(taxa.types,data.frame(horizons=h,Ns=length(Ns),one_t=length(one_t),two_t_a=length(two_t_a),two_t_b=length(two_t_b),three_t=length(three_t),Pt=length(Pt)))

  } # 1
  return(taxa.types)
  #eof
}

#' Categorize per interval taxon types using the gap-filler approach (Alroy, 2014)
#'
#' @details
#' Taxa types are detailed in Alroy (2014) \cr
#' 1. taxa sampled at all in the focal bin (Ns) \cr
#' 2. taxa sampled in a bin but not immediately before or after (one-timers, or 1t) \cr
#' 3a. taxa sampled immediately before and within the ith bin (two-timers, or 2ti) * referred to as two_t_a \cr
#' 3b. or within and immediately after the ith bin (2ti+1) * referred to as two_t_b \cr
#' 4. taxa sampled in three consecutive bins (three-timers, or 3t) \cr
#' 5. taxa sampled before and after but not within a bin (part-timers, or Pt) \cr
#' 6. taxa sampled in i − 1 and i + 2 but not i + 1 (gap-fillers, a) * referred to as gf_a # may also be sampled in i \cr
#' 7. taxa sampled in i + 1 and i - 2 but not i - 1 (gap-fillers, b) * referred to as gf_b # may also be sampled in i \cr
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @param basin.age Maximum age of the oldest stratigraphic interval
#' @param strata Number of stratigraphic horizons
#' @return dataframe of per interval taxon types.
#' @export
interval.types.gf<-function(fossils,basin.age,strata) {
  fossils<-fossils
  basin.age<-basin.age
  strata<-strata

  lineages<-unique(fossils$sp)

  s1=basin.age/strata # max age of youngest horizon/horizon length
  horizons<-seq(s1, basin.age, length=strata)

  taxa.types<-data.frame(horizons=numeric(),Ns=numeric(),one_t=numeric(),two_t_a=numeric(),two_t_b=numeric(),three_t=numeric(),Pt=numeric(),gf_a=numeric(),gf_b=numeric()) # data fram for taxon types

  # during each horizon
  for(h in horizons){ # 1
    ha=h+s1 # horizon before
    ha2=ha+s1
    hb=h-s1 # horizon after
    hb2=hb-s1

    Ns=c()  # total number of taxa sampled in the bin
    one_t=c() # one timers
    two_t_a=c() # two timers = extinction event
    two_t_b=c() # two timers = speciation event
    three_t=c() # three timers
    Pt=c() # part timers
    gf_a=c() # gap fillers = imaginary extinction event
    gf_b=c() # gap fillers = imaginary speciation event

    # for each lineage
    for(l in lineages){ # 2

      # bins in which lineage l occurrs
      i=fossils[which(fossils$sp==l),]$h

      # do I exist in the bin?
      if(any(i==as.character(h))){

        Ns=c(Ns,l)

        # do I exist in the previous bin?
        if(any(i==as.character(ha))){

          # do I also exist in the next bin?
          if(any(i==as.character(hb))){
            three_t=c(three_t,l)
          }
          else {
            two_t_a=c(two_t_a,l)
            if(any(i==as.character(hb2))){
              gf_a=c(gf_a,l)
            }
          }
        }
        # do I exist in the next bin but not the one before?
        else if (any(i==as.character(hb))){
          two_t_b=c(two_t_b,l)
          if(any(i==as.character(ha2))){
            gf_b=c(gf_b,l)
          }
        }

        # otherwise record taxon as one timer
        else {
          one_t=c(one_t,l)
        }

      }

      #else, do I exist in the bin before and after?
      else if( (any(i==as.character(ha))) & (any(i==as.character(hb))) ){
        Pt=c(Pt,l)
      }

      # define gap fillers a and b
      else if ( (any(i==as.character(ha))) & (any(i==as.character(hb2))) ){
        gf_a=c(gf_a,l)
      }
      else if ( (any(i==as.character(hb))) & (any(i==as.character(ha2))) ){
        gf_b=c(gf_b,l)
      }

    } # 2

    two_t_a=c(two_t_a,three_t)
    two_t_b=c(two_t_b,three_t)

    taxa.types<-rbind(taxa.types,data.frame(horizons=h,Ns=length(Ns),one_t=length(one_t),two_t_a=length(two_t_a),two_t_b=length(two_t_b),three_t=length(three_t),Pt=length(Pt),gf_a=length(gf_a),gf_b=length(gf_b)))

  } # 1
  return(taxa.types)
}

## Calculate speciation and extinction rates

#' Calculate speciation and extinction rates using the boundary crosser approach
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @param basin.age Maximum age of the oldest stratigraphic interval
#' @param strata Number of stratigraphic horizons
#' @param continuous If TRUE calculate continuous rates
#' @param return.intervals If TRUE return per interval estimates
#' @return named list with the overall speciation rate, overall extinction rate and a dataframe of per interval estimtes if return.intervals = TRUE.
#' Note this approach cannot estimate rates for the first interval.
#' @export
boundary.crosser.rates<-function(fossils,basin.age,strata,continuous=T,return.intervals=F) {
  fossils<-fossils
  basin.age<-basin.age
  strata<-strata
  continuous<-continuous
  return.intervals<-return.intervals

  taxa.types = interval.types.bc(fossils,basin.age,strata)

  s1=basin.age/strata # equivalent to delta t
  o.e.rates<-data.frame(horizons=numeric(),NFl=numeric(),NFt=numeric(),Nbl=numeric(),Nbt=numeric(),p=numeric(),q=numeric()) # dataframe for FAs

  k.total=0 # range through taxa
  n.total=0 # range through taxa + specation events

  e.total=0 # range through taxa
  Ntot.total=0 # range through taxa + extinction events

  for(h in 1:length(taxa.types[,1])){
    NFl=taxa.types[,2][h]
    NFt=taxa.types[,3][h]
    Nbl=taxa.types[,4][h]
    Nbt=taxa.types[,5][h]


    k.total=k.total+Nbt
    n.total=n.total+Nbt+NFt

    e.total=e.total+Nbt
    Ntot.total=Ntot.total+Nbt+Nbl

    if(Nbt > 0){
      if(!continuous){
        # calculate orgination rates
        p.hat=-log(Nbt/(Nbt+NFt))
        p.hat=round(p.hat,3)
        # calculate extinction rates
        q.hat=-log(Nbt/(Nbt+Nbl))
        q.hat=round(q.hat,3)
      }
      else {
        # calculate orgination rates
        p.hat=-log(Nbt/(Nbt+NFt))/s1
        p.hat=round(p.hat,3)
        # calculate extinction rates
        q.hat=-log(Nbt/(Nbt+Nbl))/s1
        q.hat=round(q.hat,3)
      }
    }
    else{
      p.hat=NaN
      q.hat=NaN
    }
    o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],NFl=NFl,NFt=NFt,Nbl=Nbl,Nbt=Nbt,p=p.hat,q=q.hat))
  }
  if(k.total > 0){
    if(!continuous){
      p.total=(-log(k.total/n.total))
      q.total=(-log(e.total/Ntot.total))
    }
    else{
      p.total=(-log(k.total/n.total))/s1
      q.total=(-log(e.total/Ntot.total))/s1
    }
  }
  else{
    p.total = NaN
    q.total = NaN
  }

  if(return.intervals)
    return(list(speciation=p.total,extinction=q.total,per.interval.rates=o.e.rates))
  else
    return(list(speciation=p.total,extinction=q.total))
  #eof
}

#' Calculate speciation and extinction rates using the uncorrected approach
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @param basin.age Maximum age of the oldest stratigraphic interval
#' @param strata Number of stratigraphic horizons
#' @param continuous If TRUE calculate continuous rates
#' @param return.intervals If TRUE return per interval estimates
#' @return named list with the overall speciation rate, overall extinction rate and a dataframe of per interval estimtes if return.intervals = TRUE.
#' @export
uncorrected.rates<-function(fossils,basin.age,strata,continuous=T,return.intervals=F) {
  fossils<-fossils
  basin.age<-basin.age
  strata<-strata
  continuous<-continuous
  return.intervals<-return.intervals

  taxa.types = interval.types.bc(fossils,basin.age,strata)

  s1=basin.age/strata # equivalent to delta t
  o.e.rates<-data.frame(horizons=numeric(),NFl=numeric(),NFt=numeric(),Nbl=numeric(),Nbt=numeric(),p=numeric(),q=numeric()) # dataframe for FAs

  k.total=0 # total number of speciation events
  n.total=0 # total diversity

  e.total=0 # total number of extinction events
  Ntot.total=0 # total diversity (note this denominator may be different in alternative models)

  for(h in 1:length(taxa.types[,1])){
    NFl=taxa.types[,2][h]
    NFt=taxa.types[,3][h]
    Nbl=taxa.types[,4][h]
    Nbt=taxa.types[,5][h]

    Ntot=NFl+Nbl+NFt+Nbt

    k.total=k.total+NFl+NFt
    n.total=n.total+Ntot

    e.total=e.total+NFl+Nbl
    Ntot.total=Ntot.total+Ntot

    if(!continuous){
      # calculate orgination rates
      p.hat=((NFl+NFt)/Ntot)
      p.hat=round(p.hat,3)
      # calculate extinction rates
      q.hat=((NFl+Nbl)/Ntot)
      q.hat=round(q.hat,3)
    }
    else{
      # calculate orgination rates
      p.hat=((NFl+NFt)/Ntot)/s1
      p.hat=round(p.hat,3)
      # calculate extinction rates
      q.hat=((NFl+Nbl)/Ntot)/s1
      q.hat=round(q.hat,3)
    }
    o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],NFl=NFl,NFt=NFt,Nbl=Nbl,Nbt=Nbt,p=p.hat,q=q.hat))
  }

  if(!continuous){
    p.total=(k.total/n.total)
    q.total=(e.total/Ntot.total)
  }
  else{
    p.total=(k.total/n.total)/s1
    q.total=(e.total/Ntot.total)/s1
  }
  if(return.intervals)
    return(list(speciation=p.total,extinction=q.total,per.interval.rates=o.e.rates))
  else
    return(list(speciation=p.total,extinction=q.total))
  #eof
}

#' Calculate speciation and extinction rates using the three-timer approach
#'
#'
#' @details
#' The overall sampling probability Ps = 3t / (3t + Pt), where 3t and Pt are summed across the entire dataset \cr
#' The per-interval speciation rate lamda = log(2ti+1/3t) + log(Ps) \cr
#' The per-interval extinction rate mu = log(2ti/3t) + log(Ps) \cr
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @param basin.age Maximum age of the oldest stratigraphic interval
#' @param strata Number of stratigraphic horizons
#' @param continuous If TRUE calculate continuous rates
#' @param return.intervals If TRUE return per interval estimates
#' @return named list with the overall speciation rate, overall extinction rate, overall sampling rate and a dataframe of per interval estimtes if return.intervals = TRUE.
#' Note this approach cannot estimate rates for the first interval.
#' @export
three.timer.rates<-function(fossils,basin.age,strata,continuous=T,return.intervals=F) {
  fossils<-fossils
  basin.age<-basin.age
  strata<-strata
  continuous<-continuous
  return.intervals<-return.intervals

  taxa.types = interval.types.3t(fossils,basin.age,strata)

  s1=basin.age/strata # equivalent to delta t
  o.e.rates<-data.frame(horizons=numeric(),p=numeric(),q=numeric()) # dataframe for p (speciation) and q (extinction) rates

  # calculate tree wide preservation
  three_t_total=sum(taxa.types$three_t)
  Pt_total=sum(taxa.types$Pt)
  Ps=three_t_total/(three_t_total+Pt_total)

  two_t_a.total=0
  two_t_b.total=0
  three_t.total=0

  for(h in 1:length(taxa.types[,1])){
    Ns=taxa.types$Ns[h]
    one_t=taxa.types$one_t[h]
    two_t_a=taxa.types$two_t_a[h]
    two_t_b=taxa.types$two_t_b[h]
    three_t=taxa.types$three_t[h]

    if(h==length(taxa.types[,1])){ # the first sampled interval from the beginning of the process
      o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],p=NaN,q=NaN))
    }
    else {
      if(three_t > 0){
        if(!continuous){
          # calculate orgination rates
          p.hat=( log(two_t_b/three_t)+log(Ps) )
          p.hat=round(p.hat,3)
          # calculate extinction rates
          q.hat=( log(two_t_a/three_t)+log(Ps) )
          q.hat=round(q.hat,3)
        }
        else{
          # calculate orgination rates
          p.hat=(( log(two_t_b/three_t)+log(Ps) )/s1 )
          p.hat=round(p.hat,3)
          # calculate extinction rates
          q.hat=(( log(two_t_a/three_t)+log(Ps) ) /s1 )
          q.hat=round(q.hat,3)
        }
      }
      else{
        p.hat = NaN
        q.hat = NaN
      }

      if( (!is.na(p.hat)) && (p.hat < 0) ){
        p.hat=0
      }
      if( (!is.na(q.hat)) && (q.hat < 0) ){
        q.hat=0
      }
      o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],p=p.hat,q=q.hat))
    }
    two_t_a.total = two_t_a.total+two_t_a
    two_t_b.total = two_t_b.total+two_t_b
    three_t.total = three_t.total+three_t
  }

  if(three_t.total > 0){
    if(!continuous){
      p.total=( log(two_t_b.total/three_t.total) + log(Ps) )
      q.total=( log(two_t_a.total/three_t.total) + log(Ps) )
    }
    else {
      p.total=( log(two_t_b.total/three_t.total) + log(Ps))/s1
      q.total=( log(two_t_a.total/three_t.total) + log(Ps))/s1
    }
  }
  else{
    p.total = NaN
    q.total = NaN
  }

  if( (!is.na(p.total)) && (p.total < 0) ){
    p.total=0
  }
  if( (!is.na(q.total)) && (q.total < 0) ){
    q.total=0
  }
  if(return.intervals)
    return(list(speciation=p.total,extinction=q.total,sampling=Ps,per.interval.rates=o.e.rates))
  else
    return(list(speciation=p.total,extinction=q.total,sampling=Ps))
  # eof
}

#' Calculate speciation and extinction rates using the gap-filler approach
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels. h = ages.)
#' @param basin.age Maximum age of the oldest stratigraphic interval
#' @param strata Number of stratigraphic horizons
#' @param return.intervals If TRUE return per interval estimates
#' @return named list with the overall speciation rate, overall extinction rate and a dataframe of per interval estimtes if return.intervals = TRUE.
#' Note this approach cannot estimate specation for the first or second interval and cannot estimate extinction for the first or last interval.
#' @export
gap.filler.rates<-function(fossils,basin.age,strata,return.intervals=F) {
  fossils<-fossils
  basin.age<-basin.age
  strata<-strata
  return.intervals<-return.intervals

  taxa.types = interval.types.gf(fossils,basin.age,strata)

  s1=basin.age/strata # equivalent to delta t
  o.e.rates<-data.frame(horizons=numeric(),p=numeric(),q=numeric()) # dataframe for p (speciation) and q (extinction) rates

  two_t_a.total=0
  two_t_b.total=0
  three_t.total=0
  p_t.total=0
  gap.filler_a.total=0
  gap.filler_b.total=0

  for(h in 1:length(taxa.types[,1])){
    Ns = taxa.types$Ns[h]
    one_t = taxa.types$one_t[h]
    two_t_a = taxa.types$two_t_a[h]
    two_t_b = taxa.types$two_t_b[h]
    three_t = taxa.types$three_t[h]
    p_t = taxa.types$Pt[h]
    gf_a = taxa.types$gf_a[h]
    gf_b = taxa.types$gf_b[h]

    if(sum(c(three_t,p_t,gf_a)) > 0){
      # calculate orgination rates
      p.hat=(log( (two_t_b + p_t) / (three_t + p_t + gf_a) )/s1 )
      p.hat=round(p.hat,3)
    }
    else{
      p.hat = NaN
    }
    if(sum(c(three_t,p_t,gf_b)) > 0){
      # calculate extinction rates
      q.hat=(log( (two_t_a + p_t) / (three_t + p_t + gf_b) )/s1 )
      q.hat=round(q.hat,3)
    }
    else{
      q.hat = NaN
    }

    if( (!is.na(p.hat)) && (p.hat < 0) ){
      p.hat=0
    }
    if( (!is.na(q.hat)) && (q.hat < 0) ){
      q.hat=0
    }

    if(h==length(taxa.types$horizons)){ # the first sampled interval after the beginning of the process
      o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],p=NaN,q=NaN))
    }
    else if(h==(length(taxa.types$horizons)-1)){ # the second sampled interval after the beginning of the process
      o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],p=NaN,q=q.hat))
    }
    else if(h==1){ # the last sampled interval before the present
      o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],p=p.hat,q=NaN))
    }
    else{
      o.e.rates<-rbind(o.e.rates,data.frame(horizons=taxa.types$horizons[h],p=p.hat,q=q.hat))
    }

    #if(!(h %in% c(1,19,20)))
      two_t_a.total = two_t_a.total+two_t_a
      two_t_b.total = two_t_b.total+two_t_b
      three_t.total = three_t.total+three_t
      p_t.total = p_t.total+p_t
      gap.filler_a.total = gap.filler_a.total+gf_a
      gap.filler_b.total = gap.filler_b.total+gf_b

  }

  if(sum(c(three_t.total,p_t.total,gap.filler_a.total)) > 0){
    p.total = ( log( (two_t_b.total + p_t.total) / (three_t.total + p_t.total + gap.filler_a.total) ) ) /s1
  }
  else{
    p.total = NaN
  }
  if(sum(c(three_t.total,p_t.total,gap.filler_b.total)) > 0){
    q.total=( log( (two_t_a.total + p_t.total)/ (three_t.total + p_t.total + gap.filler_b.total) ) ) /s1
  }
  else{
    q.total = NaN
  }

  if( (!is.na(p.total)) && (p.total < 0) ){
    p.total=0
  }
  if( (!is.na(q.total)) && (q.total < 0) ){
    q.total=0
  }

  if(return.intervals)
    return(list(speciation=p.total,extinction=q.total,per.interval.rates=o.e.rates))
  else
    return(list(speciation=p.total,extinction=q.total))
  #eof
}

# Note functions return NA if insufficient information is available to estimate rates (not zero).
# e.g. zero means the estimated rate(s) = zero.

