#' @title calib.stat
#' 
#' @description
#' Calibration statistic for censored time-to-event data
#' 
#' @details
#' Compute the calibration statistic for censored time-to-event data
#' 
#' @param p Predicted event probablities
#' @param T.test Observation times from the test set
#' @param C.test Event indicators from the test set
#' @param cutpts Vector of cutpoints defining risk categories
#' @param t Time horizon of predictions
#' 
#' @return The calibration statistic
calib.stat <- function(p,T.test,C.test,cutpts,t) {
  risk.class <- cut(p,cutpts,labels=FALSE)
  lev.stats <- sapply(1:(length(cutpts)-1),function(f) {
    ind <- which(risk.class==f)
    sf <- survfit(Surv(T.test[ind],C.test[ind])~1)
    ind.surv <- max(which(sf$time<=t))
    p.KM <- sf$surv[ind.surv]
    ##print(c(cutpts[f],p.KM,mean(p[ind],na.rm=TRUE),sqrt(S.KM$SKMV[ind.surv])))
    (mean(p[ind],na.rm=TRUE) - p.KM)^2/sf$std.err[ind.surv]
  })
  
  sum(lev.stats,na.rm=TRUE)
}

#' @title calib.table
#' 
#' @description
#' Calibration table for censored time-to-event data
#' 
#' @details
#' Compute the calibration table for censored time-to-event data
#' 
#' @param LP List of vectors of predicted event probablities.
#' @param T.test Observation times from the test set.
#' @param C.test Event indicators from the test set.
#' @param cutpts Vector of cutpoints defining risk categories.
#' @param t Time horizon of predictions.
#' @param modelnames (Optional) Names of the models in \code{LP}, for pretty printing. 
#' @return A calibration table, with one row for each category defined by \code{cutpts}.
calib.table <- function(LP,T.test,C.test,cutpts,t,modelnames=rep("",length(LP))) {
  LStats <- lapply(LP,function(p) {
    risk.class <- cut(p,cutpts,labels=FALSE)
    lev.stats <- sapply(1:(length(cutpts)-1),function(f) {
      ind <- which(risk.class==f)
      sf <- survfit(Surv(T.test[ind],C.test[ind])~1)
      ind.surv <- max(which(sf$time<=t))
      p.KM <- sf$surv[ind.surv]
      ##print(c(cutpts[f],p.KM,mean(p[ind],na.rm=TRUE),sqrt(S.KM$SKMV[ind.surv])))
      c(n=length(ind),pbar=mean(p[ind],na.rm=TRUE),pKM=1-p.KM)
    })
    return(t(lev.stats))    
  })

  final.tab <- do.call("cbind",LStats)
  return(final.tab)
}
