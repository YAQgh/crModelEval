riskEst=function(data, predicttime, formula, model){
  if(is.list(formula) & length(formula)==2){
    fm1=as.formula(formula[[1]])
    fm2=as.formula(formula[[2]])
  }
  else if(is.list(formula) & length(formula)==1){
    fm1=fm2=as.formula(formula[[1]])
  }
  else{
    fm1=fm2=as.formula(formula)
  }
  if(model=="cox"){
    mod1=mod2=riskRegression::CSC(list(fm1, fm2), data=data)
    names1=names(mod1$models$`Cause 1`$coefficients)
    names2=names(mod1$models$`Cause 2`$coefficients)
  }
  else if(model=="fg"){
    mod1=riskRegression::FGR(fm1, data=data, cause=1)
    mod2=riskRegression::FGR(fm2, data=data, cause=2)
    names1=names(mod1$crrFit$coef)
    names2=names(mod1$crrFit$coef)
  }
  else{
    stop("model should be cox or fg.")
  }
  cif1=riskRegression::predictRisk(mod1, newdata=data[, names1], times=predicttime, cause=1)
  cif2=riskRegression::predictRisk(mod2, newdata=data[, names2], times=predicttime, cause=2)
  surv=1-cif1-cif2
  return(cbind(cif1, cif2, surv))
}

bootfun=function(loc, data, predicttime, formula1, formula2, model, cause, method){
  obsdata=data[loc,]
  if(length(model)==1){
    risk1=riskEst(obsdata, predicttime, formula1, model)
    risk2=riskEst(obsdata, predicttime, formula2, model)
  }
  else if(length(model)==2){
    risk1=riskEst(obsdata, predicttime, formula1, model[1])
    risk2=riskEst(obsdata, predicttime, formula2, model[2])
  }
  else{
    stop("model can only be a vector of of one or two components.")
  }

  # extract names for time and delta
  if(is.list(formula1) & length(formula1)==2){
    fm1=as.formula(formula1[[1]])
    fm2=as.formula(formula1[[2]])
  }
  else if(is.list(formula1) & length(formula1)==1){
    fm1=fm2=as.formula(formula1[[1]])
  }
  else{
    fm1=fm2=as.formula(formula1)
  }
  tempmod=riskRegression::CSC(list(fm1, fm2), data=data)
  timename=attributes(tempmod$response)$dimnames[[2]][1]
  deltaname=attributes(tempmod$response)$dimnames[[2]][2]

  G_hat=comp_G_hat(dat=data.frame(x=obsdata[,timename], delt_eps=obsdata[,deltaname]), tau=predicttime)
  if(cause==0){
    if(method=="pdi"){
      pdi1=pdi(predicttime,risk1,obsdata[,timename],delt_eps=obsdata[,deltaname],G_hat,3)
      pdi2=pdi(predicttime,risk2,obsdata[,timename],delt_eps=obsdata[,deltaname],G_hat,3)
      diff=mean(pdi1)-mean(pdi2)
      return(diff)
    }
    else if(method=="ExC"){
      exc1=ExC(obsdata[,timename],obsdata[,deltaname],risk1[,1],risk1[,2],risk1[,3],predicttime)
      exc2=ExC(obsdata[,timename],obsdata[,deltaname],risk2[,1],risk2[,2],risk2[,3],predicttime)
      diff=exc1-exc2
      return(diff)
    }
    else{
      stop("method should be auc, bs, pdi or exc.")
    }
  }
  else{
    if(method=="pdi"){
      pdi1=pdi(predicttime,risk1,obsdata[,timename],delt_eps=obsdata[,deltaname],G_hat,3)
      pdi2=pdi(predicttime,risk2,obsdata[,timename],delt_eps=obsdata[,deltaname],G_hat,3)
      diff=pdi1[cause]-pdi2[cause]
      return(diff)
    }
    else if(method=="ExC"){
      exc1=ExCl(obsdata[,timename],obsdata[,deltaname],risk1[,1],risk1[,2],risk1[,3],predicttime,cause)
      exc2=ExCl(obsdata[,timename],obsdata[,deltaname],risk2[,1],risk2[,2],risk2[,3],predicttime,cause)
      diff=exc1-exc2
      return(diff)
    }
    else{
      stop("method should be pdi or exc.")
    }
  }
}

#' Bootstrap Confidence Interval for Model Comparison
#'
#' Confidence interval of difference between two models based on bias-corrected and accelerated bootstrap method
#'
#' @param data observed data. should include time, cause indicator with 0 for censoring, and all covariates to be included in the prognostic models.
#' @param predicttime time of model evaluation
#' @param formula1 list of Hist formulas for cause-1 and cause-2 events in the first model. cause-1 and cause-2 can have different formulas. If only one formula is given, will be used for both cause-1 and cause-2.
#' @param formula2 list of Hist formulas for cause-1 and cause-2 events in the second model. cause-1 and cause-2 can have different formulas. If only one formula is given, will be used for both cause-1 and cause-2.
#' @param model vector of two prognostic models to be fitted for comparison. "cox" and "fg" can be fitted. If only one model is specified, will be used for both models for comparison.
#' @param cause 1 or 2 for cause-specific evaluation. 0 for overall assessment.
#' @param method either "pdi" or "ExC"
#' @param nboot number of bootstrap samples to be taken.
#' @param alpha specify confidence level.
#'
#' @return bias-corrected and accelerated confidence interval of difference between two models at specified confidence level.
#' @export
#'
#' @examples
bcaCI=function(data, predicttime, formula1, formula2, model, cause, method, nboot, alpha){
  # formula1: list of two formulas used for cause-1 and cause-2 event for the first model
  # formula2: list of two formulas used for cause-1 and cause-2 event for the second model
  # if formula 1 and formula 2 has length one, then cause-1 and cause-2 rely on same group of cov
  # formula should be given by character string adapted for CSC and FGR in riskRegression
  # model should be a vector of character string, specifying cox or fg for model1 and model2
  # cause=0 for overall evaluation.
  return(bootstrap::bcanon(1:nrow(data),nboot=nboot,theta=bootfun,data=data,predicttime=predicttime, formula1=formula1, formula2=formula2, model=model, cause=cause, method=method, alpha=alpha)$confpoints)
}
