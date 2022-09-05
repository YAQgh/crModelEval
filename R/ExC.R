#' Extended Concordance Index for Overall Evaluation
#'
#' Estimate extended concordance index that evaluate the ability of prognostic models to discriminate subjects from all groups simultaneously and to predict event status under competing risks setting.
#'
#' @param time The observed time or followed up time for right censored data.
#' @param delta The observed status indicator. 0 for censoring.
#' @param cif1 Cause-1 CIF function estimated at predict time.
#' @param cif2 Cause-2 CIF function estimated at predict time.
#' @param surv Survival function estimated at predict time.
#' @param predicttime A single time point of prediction or evaluation.
#'
#' @return Estimated concordance index for overall evaluation at predicttime.
#' @export
#'
#' @examples
ExC=function(time, delta, cif1, cif2, surv, predicttime){
  # mydata supposed to be a dataframe
  # mydata[,1] supposed to be observed time
  # mydata[,2] supposed to be cause indicator
  # mydata[, 3:5] supposed to be CIF1, CIF2, S at predicttime for each subject
  mydata=data.frame(time, delta, cif1, cif2, surv)

  colnames(mydata)=c("obstime","causesta","cif1","cif2","s")
  mydata$trueclass=ifelse((mydata$obstime<=predicttime & mydata$causesta==1), 1, ifelse((mydata$obstime<=predicttime & mydata$causesta==2),2,ifelse((mydata$obstime>predicttime),3,NA)))
  # mydata$predsta=apply(mydata[,3:5],1,which.max)
  # mydata$predsta[mydata$predsta==3]=0
  # mydata$correct=as.numeric(mydata$trueclass==mydata$predsta)
  # for those censored before predicttime, trueclass is NA
  cpmdata=mydata[order(mydata$obstime),] # order by obseved time
  cpmdataobs<-na.omit(cpmdata) # remove those censored before predicttime
  n.nmiss=length(cpmdataobs$trueclass) # number of subjects not missing

  obsclass1=which(cpmdataobs$trueclass==1)                      ## T<t0 and event=1 and T<C
  obsclass2=which(cpmdataobs$trueclass==2)                      ## T<t0 and event=2 and T<C
  obsclass3=which(cpmdataobs$trueclass==3)                      ## T>t0

  nclass1=length(obsclass1)
  nclass2=length(obsclass2)
  nclass3=length(obsclass3)

  Tevent1=cpmdataobs$obstime[obsclass1] # event time for event=1
  Tevent2=cpmdataobs$obstime[obsclass2] # event time for event=2
  Tevent3=rep(predicttime, nclass3) # censored at predicttime
  Teventobs=c(Tevent1, Tevent2, Tevent3)
  obsevent.npts=length(Teventobs)

  # estimate survival function of censoring
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

  censortimepool=matrix(rep(censoring.time, each=obsevent.npts), obsevent.npts, length(censoring.time))
  # each row is a copy of censoring.time, obsevent.npts number of rows

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
  # for n event times, from the first n-1 times, find the first one that is >= the nth time.

  eventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index)
  probcen.event=censoring.survival[eventatcensordist.index]

  probcen.event1=probcen.event[1:nclass1]
  probcen.event2=probcen.event[(nclass1+1):(nclass1+nclass2)]
  noeventatcensordist.index = apply(cbind(censortimepool, Teventobs), 1, max_index_noevent)
  probcen.noevent=censoring.survival[noeventatcensordist.index]
  probcen.event3=probcen.noevent[(nclass1+nclass2+1):(nclass1+nclass2+nclass3)]
  probcensorweig=kronecker(probcen.event1%*%t(probcen.event2),probcen.event3)

  # class1comp=cpmdataobs$correct[obsclass1]
  # class2comp=cpmdataobs$correct[obsclass2]
  # class3comp=cpmdataobs$correct[obsclass3]
  # Num.IPCW=sum(kronecker(class1comp%*%t(class2comp),class3comp)/probcensorweig)
  class1comp=cpmdataobs[obsclass1,]
  class2comp=cpmdataobs[obsclass2,]
  class3comp=cpmdataobs[obsclass3,]
  Num.IPCW=0
  for(i in 1:nclass1){
    for(j in 1:nclass2){
      for(k in 1:nclass3){
        F1i=class1comp$cif1[i]
        F2i=class1comp$cif2[i]
        Si=class1comp$s[i]
        F1j=class2comp$cif1[j]
        F2j=class2comp$cif2[j]
        Sj=class2comp$s[j]
        F1k=class3comp$cif1[k]
        F2k=class3comp$cif2[k]
        Sk=class3comp$s[k]
        # if(F1i>=F1j & F1i>=F1k & F2j>=F2i & F2j>=F2k & Sk>=Si & Sk>=Sj){
        #   # pdiout[(i-1)*nclass2+j,k]=1
        #   Num.IPCW=Num.IPCW+1/(probcen.event1[i]*probcen.event2[j]*probcen.event3[k])
        # }
        # else{
        #   next
        # }
        c1=(F1i>=F1j & F1i>=F1k)/(1+(F1i==F1j)+(F1i==F1k))
        c2=(F2j>=F2i & F2j>=F2k)/(1+(F2j==F2i)+(F2j==F2k))
        c3=(Sk>=Si & Sk>=Sj)/(1+(Sk==Si)+(Sk==Sj))
        Num.IPCW=Num.IPCW+c1*c2*c3/(probcen.event1[i]*probcen.event2[j]*probcen.event3[k])
      }
    }
  }

  event1n=rep(1, nclass1)
  event2n=rep(1, nclass2)
  event3n=rep(1, nclass3)
  Demon.IPCW=sum(kronecker(event1n%*%t(event2n),event3n)/probcensorweig)

  weightedHUM=Num.IPCW/Demon.IPCW
  return(weightedHUM)
}
