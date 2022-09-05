#' Extended Concordance Index for Cause-specific Evaluation
#'
#' Estimate extended concordance index that evaluate the ability of a prognostic model to distinguish between the cause of interest and the healthy control and to predict event status under competing risks setting.
#'
#' @param time The observed time or followed up time for right censored data.
#' @param delta The observed status indicator. 0 for censoring.
#' @param cif1 Cause-1 CIF function estimated at predict time.
#' @param cif2 Cause-2 CIF function estimated at predict time.
#' @param surv Survival function estimated at predict time.
#' @param predicttime A single time point of prediction or evaluation.
#' @param cause Specify cause of interest.
#'
#' @return Estimated concordance index for cause-specific evaluation at predicttime.
#' @export
#'
#' @examples
ExCl=function(time, delta, cif1, cif2, surv, predicttime, cause){
  mydata=data.frame(time, delta, cif1, cif2, surv)
  colnames(mydata)=c("obstime","causesta", "cif1", "cif2", "s")
  mydata$trueclass=ifelse((mydata$obstime<=predicttime & mydata$causesta==1), 1, ifelse((mydata$obstime<=predicttime & mydata$causesta==2),2,ifelse((mydata$obstime>predicttime),3,NA)))
  mydata$predsta=apply(mydata[,3:5],1,which.max)
  mydata$correct=as.numeric(mydata$trueclass==mydata$predsta)
  cpmdata=mydata[order(mydata$obstime),]
  cpmdataobs<-na.omit(cpmdata)
  n.nmiss=length(cpmdataobs$trueclass)

  obsclass1=which(cpmdataobs$trueclass==1)                      ## T<t0 and event=1 and T<C
  obsclass2=which(cpmdataobs$trueclass==2)                      ## T<t0 and event=2 and T<C
  obsclass3=which(cpmdataobs$trueclass==3)

  nclass1=length(obsclass1)
  nclass2=length(obsclass2)
  nclass3=length(obsclass3)

  Tevent1=cpmdataobs$obstime[obsclass1] # event time for event=1
  Tevent2=cpmdataobs$obstime[obsclass2] # event time for event=2
  Tevent3=rep(predicttime, nclass3) # censored at predicttime
  Teventobs=c(Tevent1, Tevent2, Tevent3)
  obsevent.npts=length(Teventobs)

  dim.time = length(cpmdata$obstime)
  censoring.index=which(cpmdata$causesta==0)
  censoring.survival=NULL # KM estimator of censoring
  censoring.time=unique(sort(cpmdata$obstime[censoring.index])) # sorted unique time of censoring
  censoring.npts=length(censoring.time)

  timepool=matrix(rep(cpmdata$obstime, each=censoring.npts),censoring.npts, dim.time) # each row is a copy of obstime, censoring.npts copies/rows
  statuspool=matrix(rep(cpmdata$causesta, each=censoring.npts),censoring.npts, dim.time) # each row is a copy of causesta, censoring.npts copies/rows
  censor.Y=apply(ifelse(timepool>=censoring.time, 1,0), 1, sum) # at risk set
  censor.d=apply(ifelse(timepool==censoring.time & statuspool==0, 1,0),1,sum) # event
  censoring.survival=cumprod(1-censor.d/censor.Y)
  censoring.time=c(0, censoring.time)
  censoring.survival=c(1,censoring.survival)
  if(censoring.survival[length(censoring.survival)]==0){
    censoring.survival[length(censoring.survival)]=censoring.survival[length(censoring.survival)-1]
  }

  censortimepool=matrix(rep(censoring.time, each=obsevent.npts), obsevent.npts, length(censoring.time))
  max_index = function(uniquetime.event_time)
  {
    vec.len=length(uniquetime.event_time)
    uniquetime = uniquetime.event_time[1:vec.len-1]
    event_time = uniquetime.event_time[vec.len]
    if(event_time<=min(uniquetime))
      return(1)
    else
      return(max(which(uniquetime<event_time)))
  }
  # for n event times, from the first n-1 points, find the last one that is <= the nth time
  max_index_noevent = function(uniquetime.event_time)
  {
    vec.len=length(uniquetime.event_time)
    uniquetime = uniquetime.event_time[1:vec.len-1]
    event_time = uniquetime.event_time[vec.len]
    if(sum(uniquetime>=event_time)==0)
      return(vec.len-1)
    else
      return(min(which(uniquetime>=event_time)))

  }

  eventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index)
  probcen.event=censoring.survival[eventatcensordist.index]

  probcen.event1=probcen.event[1:nclass1]
  probcen.event2=probcen.event[(nclass1+1):(nclass1+nclass2)]
  noeventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index_noevent)
  probcen.noevent=censoring.survival[noeventatcensordist.index]
  probcen.event3=probcen.noevent[(nclass1+nclass2+1):(nclass1+nclass2+nclass3)]

  class1comp=cpmdataobs[obsclass1,]
  class2comp=cpmdataobs[obsclass2,]
  class3comp=cpmdataobs[obsclass3,]

  event1n=rep(1, nclass1)
  event2n=rep(1, nclass2)
  event3n=rep(1, nclass3)

  if(cause==1){
    Num.IPCW=0
    for(i in 1:nclass1){
      for(k in 1:nclass3){
        F1i=class1comp$cif1[i]
        Si=class1comp$s[i]
        F1k=class3comp$cif1[k]
        Sk=class3comp$s[k]
        # if(F1i>=F1k & Si<=Sk){
        #   Num.IPCW=Num.IPCW+1/(probcen.event1[i]*probcen.event3[k])
        # }
        # else{
        #   next
        # }
        ci=(F1i>=F1k)/(1+(F1i==F1k))
        ck=(Sk>=Si)/(1+(Sk==Si))
        Num.IPCW=Num.IPCW+ci*ck/(probcen.event1[i]*probcen.event3[k])
      }
    }
    probcensorweig=probcen.event1%*%t(probcen.event3)
    Denom.IPCW=sum((event1n%*%t(event3n))/probcensorweig)
    return(Num.IPCW/Denom.IPCW)
  }
  else if(cause==2){
    Num.IPCW=0
    for(j in 1:nclass2){
      for(k in 1:nclass3){
        F2j=class2comp$cif2[j]
        Sj=class2comp$s[j]
        F2k=class3comp$cif2[k]
        Sk=class3comp$s[k]
        # if(F2j>=F2k & Sj<=Sk){
        #   Num.IPCW=Num.IPCW+1/(probcen.event2[j]*probcen.event3[k])
        # }
        # else{
        #   next
        # }
        cj=(F2j>=F2k)/(1+(F2j==F2k))
        ck=(Sk>=Sj)/(1+(Sk==Sj))
        Num.IPCW=Num.IPCW+cj*ck/(probcen.event2[j]*probcen.event3[k])
      }
    }
    probcensorweig=probcen.event2%*%t(probcen.event3)
    Denom.IPCW=sum((event2n%*%t(event3n))/probcensorweig)
    return(Num.IPCW/Denom.IPCW)
  }
  else{
    return(NA)
  }
}
