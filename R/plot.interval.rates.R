
# the following functions can be applied to output generate using the function replicates.main

# function to plot the output from replicates.main: x-axis = absolute & realtive time, multiple options for plotting the the proportion of useful replicates
#' @export
draw.p.q.est.ab<-function(in1,in2,title="Title",col="grey75",p.true=lambda,q.true=mu,psv=F,sampling,abs=F,trees,upp=1){
  h.rates<-in1 # interval estimates
  t.rates<-in2 # rate estimates
  title<-title # graph title
  col<-col # polygon color
  p.true<-p.true # true speciation
  q.true<-q.true  # true extinction
  # for variable extinction rates
  if(length(p.true) > 1) {
    p.true=p.true[1]
  }
  # plot sampling value
  if(psv){
    ps.vals<-sampling
  }
  # if ab = T, use absolute not relative time -
  # option not active
  #if(abs){
  #  trees<-trees
  #}
  upp<-upp # options = 1-4 for plotting % of useful replicates
  # -------
  if (col=="grey75"){
    col1="grey75"
    col2="red"
  }
  else if(col=="red"){
    col1=rgb(255/255,51/255,51/255, alpha=0.5)
    col2="red"
  }
  else if(col=="blue"){
    col1=rgb(0/255,128/255,255/255, alpha=0.5)
    col2="blue"
  }
  else if(col=="green"){
    col1=rgb(76/255,153/255,0, alpha=0.5)
    col2="green"
  }
  else if(col=="orange"){
    col1=rgb(255/255,200/255,0/255, alpha=0.5)
    col2=rgb(255/255,205/255,0/255)
  }
  else if(col=="purple"){
    col1=rgb(142/255,104/255,158/255, alpha=0.5)
    col2="purple"
  }
  else {
    col1=col
    col2=col
  }
  # color of true lambda, mu values
  col.tb=rgb(0,0,0, alpha=0.6)

  # -------
  h.rates<-h.rates[,c(40:1)]
  # -------
  # part 1.
  # calculate tree wide mean speciation & extinction rates
  # sepciation (=p)
  # replace Inf estimates values with NaN
  t.rates$p.total[t.rates$p.total==Inf]<-NaN
  valid=t.rates$p.total[which(t.rates$p.total!="NaN")]
  # calculate the mean & CIs
  #s=round(sd(t.rates$p.total,na.rm=T),2)
  #n=length(t.rates$p.total)
  if(length(valid) > 1){
    m=round(mean(t.rates$p.total,na.rm=T),2)
    hpds=coda::HPDinterval(coda::as.mcmc(t.rates$p.total))
    lower=round(hpds[1],2)
    upper=round(hpds[2],2)
    tree.p=paste("Tree wide estimate = ",m," (",lower,", ",upper,")",sep="")
  }
  else{
    tree.p=paste("Tree wide estimate = NA")
  }
  # extinction (=q)
  t.rates$q.total[t.rates$q.total==Inf]<-NaN
  valid=t.rates$q.total[which(t.rates$q.total!="NaN")]
  if(length(valid) > 1){
    m=round(mean(t.rates$q.total,na.rm=T),2)
    hpds=coda::HPDinterval(coda::as.mcmc(t.rates$q.total))
    lower=round(hpds[1],2)
    upper=round(hpds[2],2)
    tree.q=paste("Tree wide estimate = ",m," (",lower,", ",upper,")",sep="")
  }
  else{
    tree.q=paste("Tree wide estimate = NA")
  }
  # -------

  # part 2.
  # calculate the mean & CIs per bin
  hp=0 # time bin/horizon counter
  p.rates<-data.frame(h=numeric(),mean=numeric(),L=numeric(),U=numeric(),use=numeric())
  q.rates<-data.frame(h=numeric(),mean=numeric(),L=numeric(),U=numeric(),use=numeric())

  for (t in 1:length(h.rates)) {
    # for extinction rates:
    if(is.wholenumber(t/2)==0) {
      hp=hp+1
      # replace Inf values with NaN ## this issue requires further attention
      h.rates[[t]][h.rates[[t]]==Inf]<-NaN
      # calculate the mean and sd for each bin, excluding NaN/missing values
      #m=mean(h.rates[[t]],na.rm=T)
      #s=sd(h.rates[[t]],na.rm=T)
      hr=na.omit(h.rates[t])
      if(length(hr[,1]) > 1) {
        m=mean(hr[,1])
        hpds=coda::HPDinterval(coda::as.mcmc(hr))
        lower=hpds[1]
        upper=hpds[2]
        useful=length(hr[,1])/length(h.rates[[t]])
      }
      else {
        lower=NaN
        upper=NaN
        m=NaN
        useful=0
      }
      q.rates<-rbind(q.rates,data.frame(h=hp,mean=m,L=lower,U=upper,use=useful))
    }
    # for speciation rates:
    else {
      # replace Inf values with NaN ## this issue requires further attention
      h.rates[[t]][h.rates[[t]]==Inf]<-NaN
      # calculate the mean for each bin, excluding NaN/missing values
      #m=mean(h.rates[[t]],na.rm=T)
      # the standard deviation for each bin, excluding NaN/missing values
      #s=sd(h.rates[[t]],na.rm=T)
      hr=na.omit(h.rates[t])
      if(length(hr[,1]) > 1) {
        m=mean(hr[,1])
        hpds=coda::HPDinterval(coda::as.mcmc(hr))
        lower=hpds[1]
        upper=hpds[2]
        useful=length(hr[,1])/length(h.rates[[t]])
      }
      else {
        lower=NaN
        upper=NaN
        useful=0
        m=NaN
      }
      p.rates<-rbind(p.rates,data.frame(h=hp,mean=m,L=lower,U=upper,use=useful))
    }
  }
  # -------

  # -------
  # part 3.
  # generate subsets of polygons
  # speciation
  polygons.p<-polygon.subsets(p.rates)
  # extinction
  polygons.q<-polygon.subsets(q.rates)
  # -------

  # -------
  # part 4.
  # plot the output
  # -------
  par(mfrow = c(2, 3)) # 2x3 layout
  par(oma=c(2, 0, 2, 0)) # to add an outer margin to the top and bottom of the graph -- bottom, left, top, right
  par(xpd=NA) # allow content to protrude into outer margin (and beyond)
  par(mar=c(1, 3.5, 0, 0.5)) # to reduce the margins around each plot - bottom, left, top, right -- this is harder to manipulate
  par(mgp=c(1.5, .5, 0)) # to reduce the spacing between the figure plotting region and the axis labels -- axis label at 1.5 rows distance, tick labels at .5 row

  # speciation rates
  # polygons
  # define ylims
  df<-p.rates
  if(length(na.omit(df$L))==0){
    max.u=p.true+.5
    min.l=0
  }
  else {
    min.l=min(df$L,na.rm=T)
    max.u=max(df$U,na.rm=T)
    if(max.u < p.true) {
      max.u=p.true+.5
    }
    if(min.l > 0) {
      min.l=0
    }
  }
  if(upp==4){
    par(fig=c(0,0.33,0.625,1))
  }
  plot(df$h, df$mean, ylim = c(min.l,max.u), xlim=c(0,max(df$h)), type="n",xlab="",ylab="Speciation rate HPD",bty="n",xaxt="n",tcl=-0.2)
  lines(df$h,y=rep(p.true,length(df$h)),col=col.tb,lwd=3,lty=3) # true speciation rate
  draw.polygons(polygons.p,col=col) # add color option
  title(title,line=0.5,cex.main=1.2,adj=0,font.main=1) # title options
  # useful proportion
  if(upp==1){
    lines(df$h,y=df$use,col="darkgrey",lwd=3,lty=1)
  }
  else if (upp==2){
    par(new=T)
    plot(x=df$h,y=df$use,axes=F,xlab=NA,ylab=NA,col="darkgrey",lwd=3,lty=1,type="l",col.axis="black")
    axis(4,col="black",at=c(0,1),las=2,tcl=-0.2)
    mtext(side=4,"% of useful replicates",line=0.5,cex=.65)
  }
  else if (upp==3){
    lines(df$h,y=df$use,col="darkgrey",lwd=3,lty=1)
    axis(4,col="black",at=c(0,1),las=2,tcl=-0.2)
    mtext(side=4,"% URs",line=0.5,cex=.65,adj=.15)
  }
  else if (upp==4){
    par(fig=c(0,0.33,0.5,0.625),new=T)
    plot(df$h, df$mean, ylim = c(0,1), xlim=c(0,max(df$h)), type="n",xlab=NA,ylab=NA,bty="n",xaxt="n",yaxt="n") # tcl=-0.2
    lines(df$h,y=df$use,col="darkgrey",lwd=3,lty=1)
    axis(side=2,at=c(0,0.5,1),tcl=-0.2)
    mtext(side=2,"% URs",cex=.65,line=1.5)
  }
  # -------
  # speciation
  # line plots
  if(upp==4){
    par(fig=c(0.33,0.66,0.5,1),new=T)
  }
  tr.no=length(h.rates[,1])
  # define the maximum limit on the y-axis
  if(any(!is.na(h.rates[1:tr.no,][seq(from=2,to=40,by=2)]))){
    max.u=(max(h.rates[1:tr.no,][seq(from=2,to=40,by=2)],na.rm=T))
    if(max.u < p.true) {
      max.u=p.true+.5
    }
  }
  else {
    max.u=p.true+.5
  }
  # rates vs relative time
  if(upp==4){
    plot(df$h, df$mean, ylim = c(min.l,max.u+1), xlim=c(0,max(df$h)), type="n",xlab="",ylab="Speciation rate",bty="n",xaxt="n",tcl=-0.2)
  }
  else {
    plot(df$h, df$mean, ylim = c(min.l,max.u+1), xlim=c(0,max(df$h)), type="n",xlab="",ylab="",bty="n",xaxt="n",tcl=-0.2)
  }
  for (t in 1:length(h.rates$t20.q)) {
    y=h.rates[t,][seq(from=2,to=40,by=2)]
    lines(x=df$h,y=y,col=col1) # type="o"
  }
  lines(df$h,y=rep(p.true,length(df$h)),col=col.tb,lwd=3,lty=3) # true speciation rate
  text(10,max.u+.5,tree.p)
  # -------
  # sampling probability
  if(psv){
    ps.col=rgb(90/255,90/255,90/255)
    tree.ps=tree.reps.ps(ps.vals)
    text(10,(max.u+.5)*0.9,tree.ps,col=ps.col)
  }
  # -------
  # rates vs absolute time
  if(upp==4){
    par(fig=c(0.66,1,0.5,1),new=T)
  }
  # df = dataframe containing the mean and CIs for p or q estimates
  # df$h = horizon IDs - this has to change to plot the rates against absolute time
  # note the y-axis limits do not need to change
  # define the maximum limit on the x-axis
  basin.age.max=round(max(trees[[2]]),1)+0.1 # define the maximum basin age
  ba.mx=plyr::round_any(basin.age.max,5,ceiling)

  # define an arbitrary set of horizon ages
  horizon.ages.arb=seq(from=basin.age.max,to=0,length.out=20)
  # initiate the plot -- note df$mean here is not meaningful
  plot(horizon.ages.arb, df$mean, ylim = c(min.l,max.u+1), xlim=c(ba.mx,0), type="n",xlab="",ylab="",bty="n",xaxt="n",tcl=-0.2)
  # plot the estimates for each replicate
  for (t in 1:length(h.rates$t20.q)) {
    basin.age=round(trees[[2]][t],1)+0.1
    s1=basin.age/strata # max age of youngest horizon/horizon length
    horizons=seq(basin.age, s1, length=strata)-(s1/2)
    y=h.rates[t,][seq(from=2,to=40,by=2)]
    lines(x=horizons,y=y,col=col1)
  }
  lines(horizon.ages.arb,y=rep(p.true,length(df$h)),col=col.tb,lwd=3,lty=3) # true speciation rate

  # -------
  # extinction rates
  # polygons
  # define ylims
  df<-q.rates
  if(length(na.omit(df$L))==0){
    max.u=p.true+.5
    min.l=0
  }
  else {
    min.l=min(df$L,na.rm=T)
    max.u=max(df$U,na.rm=T)
    if(max.u < p.true) {
      max.u=p.true+.5
    }
    if(min.l > 0) {
      min.l=0
    }
  }
  if(upp==4){
    par(fig=c(0,0.33,0.125,0.5),new=T)
  }
  if(upp==4){
    plot(df$h, df$mean, ylim = c(min.l,max.u), xlim=c(0,max(df$h)), type="n",xlab="",ylab="Extinction rate HPD",bty="n",xaxt="n",tcl=-0.2)
  }
  else {
    plot(df$h, df$mean, ylim = c(min.l,max.u), xlim=c(0,max(df$h)), type="n",xlab="Time bin ID",ylab="Extinction rate HPD",bty="n",tcl=-0.2)
  }
  if(length(q.true)==1){
    lines(df$h,y=rep(q.true,length(df$h)),col=col.tb,lwd=3,lty=3) # true extinction
  }
  else {
    lines(x=c(0,1,2),y=rep(q.true[1],3),col=col.tb,lwd=3,lty=3)
    lines(x=c(18,19,20),y=rep(q.true[1],3),col=col.tb,lwd=3,lty=3)
  }
  draw.polygons(polygons.q,col=col) # add color option
  # useful proportion
  if(upp==1){
    lines(df$h,y=df$use,col="darkgrey",lwd=3,lty=1)
  }
  else if (upp==2){
    par(new=T)
    plot(x=df$h,y=df$use,axes=F,xlab=NA,ylab=NA,col="darkgrey",lwd=3,lty=1,type="l",col.axis="black")
    axis(4,col="black",at=c(0,1),las=2,tcl=-0.2)
    mtext(side=4,"% of useful replicates",line=0.5,cex=.65)
  }
  else if (upp==3){
    lines(df$h,y=df$use,col="darkgrey",lwd=3,lty=1)
    axis(4,col="black",at=c(0,1),las=2,tcl=-0.2)
    mtext(side=4,"% URs",line=0.5,cex=.65,adj=.15)
  }
  else if (upp==4){
    par(fig=c(0,0.33,0,0.125),new=T)
    #plot(df$h, df$mean, ylim = c(0,1), xlim=c(0,max(df$h)), type="n",xlab="",ylab="% URs",bty="n",xaxt="n",tcl=-0.2)
    plot(df$h, df$mean, ylim = c(0,1), xlim=c(0,max(df$h)), type="n",xlab="Time bin ID (relative tree age)",ylab=NA,bty="n",yaxt="n",tcl=-0.2)
    lines(df$h,y=df$use,col="darkgrey",lwd=3,lty=1)
    axis(side=2,at=c(0,0.5,1),tcl=-0.2)
    mtext(side=2,"% URs",cex=.65,line=1.5)
  }
  # -------
  #  extinction
  # line plots
  # rates vs relative time
  if(upp==4){
    par(fig=c(0.33,0.66,0,0.5),new=T)
  }
  if(any(!is.na(h.rates[1:tr.no,][seq(from=1,to=40,by=2)]))){
    max.u=max(h.rates[1:tr.no,][seq(from=1,to=40,by=2)],na.rm=T)
    if(max.u < p.true) {
      max.u=p.true+.5
    }
  }
  else {
    max.u=p.true+.5
  }
  if(upp==4){
    plot(df$h, df$mean, ylim = c(min.l,max.u+1), xlim=c(0,max(df$h)), type="n",xlab="Time bin ID (relative tree age)",ylab="Extinction rate",bty="n",tcl=-0.2)
  }
  else{
    plot(df$h, df$mean, ylim = c(min.l,max.u+1), xlim=c(0,max(df$h)), type="n",xlab="Time bin ID (relative tree age)",ylab="",bty="n",tcl=-0.2)
  }
  for (t in 1:length(h.rates$t20.q)) {
    y=h.rates[t,][seq(from=1,to=40,by=2)]
    lines(x=df$h,y=y,col=col1)
  }
  if(length(q.true)==1){
    lines(df$h,y=rep(q.true,length(df$h)),col=col.tb,lwd=3,lty=3) # true extinction
  }
  else {
    lines(x=c(0,1,2),y=rep(q.true[1],3),col=col.tb,lwd=3,lty=3)
    lines(x=c(18,19,20),y=rep(q.true[1],3),col=col.tb,lwd=3,lty=3)
  }
  text(10,max.u+.5,tree.q)

  # rates vs absolute time
  if(upp==4){
    par(fig=c(0.66,1,0,0.5),new=T)
  }
  # the maximum limit on the x-axis is define above
  # initiate the plot -- note df$mean here is not meaningful
  plot(horizon.ages.arb, df$mean, ylim = c(min.l,max.u+1), xlim=c(ba.mx,0), type="n",xlab="Absolute tree age",ylab="",bty="n",tcl=-0.2)
  # plot the estimates for each replicate
  for (t in 1:length(h.rates$t20.q)) {
    basin.age=round(trees[[2]][t],1)+0.1
    s1=basin.age/strata # max age of youngest horizon/horizon length
    horizons=seq(basin.age, s1, length=strata)-(s1/2)
    y=h.rates[t,][seq(from=1,to=40,by=2)]
    lines(x=horizons,y=y,col=col1)
  }
  if(length(q.true)==1){
    lines(horizon.ages.arb,y=rep(q.true,length(df$h)),col=col.tb,lwd=3,lty=3) # true extinction rate
  }
  else if(length(q.true)==2){
    lines(x=c(0,5),y=rep(q.true[1],2),col=col.tb,lwd=3,lty=3)
    lines(x=c(5.05,basin.age.max),y=rep(q.true[1],2),col=col.tb,lwd=3,lty=3) # I think x should really be based on horizon.ages.arb -- e.g. basin.age.max, basin.age.max-1, basin.age.max-2
  }
  else if(length(q.true)==3){
    lines(x=c(0,5,5.1,5.05,basin.age.max),y=c(q.true[1],q.true[1],q.true[2],q.true[1],q.true[1]),col=col.tb,lwd=3,lty=3) # x final should really be based on
  }
  # -------
  #eof
}

# function to generate subsets of polygons
#' @export
polygon.subsets<-function(input){
  df<-input

  k=1 # polygon counter
  polygons=list()
  g<-data.frame(h=numeric(),mean=numeric(),L=numeric(),U=numeric())

  for(i in 1:length(df$mean)) {

    m=df$mean[i]

    # if mean != NaN collect the polygon info
    if(is.nan(m)==0){
      # if mean < Inf: add values to polygon k
      if ((m > 0) & (m < Inf)){
        g=rbind(g,data.frame(h=df$h[i],mean=df$mean[i],L=df$L[i],U=df$U[i]))
      }
      # else create a new polygon
      else {
        polygons[[k]]=g
        k=k+1
        g<-data.frame(h=numeric(),mean=numeric(),L=numeric(),U=numeric())
      }
    }
    else {
      polygons[[k]]=g
      k=k+1
      g<-data.frame(h=numeric(),mean=numeric(),L=numeric(),U=numeric())
    }

    if(i==length(df$mean)) {
      polygons[[k]]=g
    }
  }
  return(polygons)
  #eof
}

# function to plot polygons
#' @export
draw.polygons<-function(polygons,col="grey75"){
  polygons<-polygons
  col<-col

  if (col=="grey75"){
    col1="grey75"
    col2="red"
  }
  else if(col=="red"){
    col1=rgb(255/255,51/255,51/255, alpha=0.5)
    col2="red"
  }
  else if(col=="blue"){
    col1=rgb(0/255,128/255,255/255, alpha=0.5)
    col2="blue"
  }
  else if(col=="green"){
    col1=rgb(76/255,153/255,0, alpha=0.5)
    col2="green"
  }
  else if(col=="orange"){
    col1=rgb(255/255,200/255,0/255, alpha=0.5)
    col2=rgb(255/255,205/255,0/255)
  }
  else if(col=="purple"){
    col1=rgb(102/255,0/255,51/255, alpha=0.5)
    col2=rgb(102/255,0/255,51/255)
  }
  else {
    col1=col
    col2=col
  }

  for (i in 1:length(polygons)) {

    sub.polygon=polygons[[i]]

    # skip polygon if n < 1 e.g. polygon is empty

    if(length(sub.polygon$h) == 1) {
      age1=sub.polygon$h-0.1
      age2=sub.polygon$h+0.1
      lower=sub.polygon$L
      upper=sub.polygon$U
      mean=sub.polygon$mean
      rect(age1,lower,age2,upper,col=col1,border=F)
      lines(c(age1,age2),c(mean,mean),lwd = 2)
      lines(c(age1,age2),c(lower,lower), col=col2)
      lines(c(age1,age2),c(upper,upper), col=col2)
    }

    if(length(sub.polygon$h) > 1) {
      polygon(c(sub.polygon$h,rev(sub.polygon$h)),c(sub.polygon$L,rev(sub.polygon$U)),col=col1, border = FALSE)
      lines(sub.polygon$h, sub.polygon$mean, lwd = 2)
      lines(sub.polygon$h, sub.polygon$U, col=col2,lty=2)
      lines(sub.polygon$h, sub.polygon$L, col=col2,lty=2)
    }

  }
  #eof
}

# function to plot sampling probability
#' @export
tree.reps.ps<-function(in1){
  ps.values<-in1

  # replace Inf estimates values with NaN
  ps.values$ps[ps.values$ps==Inf]<-NaN
  # calculate the mean
  m = round(mean(ps.values$ps,na.rm=T),2)
  # calculate the standard deviation
  s = round(sd(ps.values$ps,na.rm=T),2)
  tree.ps=paste("Pr(s) =",m,"+/-",s)

  return(tree.ps)
  #eof
}

is.wholenumber<-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


