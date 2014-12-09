#' @title CNB
#' 
#' @description
#' Censored Naive Bayes
#' 
#' @details
#' This function implements two versions of the Censored Naive Bayes technique as described in the technical report available at http://z.umn.edu/julianw-cnb
#' 
#' @param X Matrix of predictors for training the model.
#' @param T Vector of observation times.
#' @param C Vector of event indicators (1 = event, 0 = censored).
#' @param SURVTIME Time horizon (on the same scale as the observation times in \code{T}) over which predictions of the survival probability are desired. 
#' For instance, if \code{SURVTIME} equals \eqn{t} and \code{X.test} equals \eqn{X}, then the function estimates \eqn{P(T \geq t | X)}.
#' @param X.test Matrix of predictors for which to obtain predictions.
#' @param scaleX Logical indicating whether or not to standardize \code{X} and \code{X.test} so that each column has zero mean and unit variance. Default is \code{TRUE}.
#' @param type Version of CNB to use. The default, \code{type="PC"}, implements Censored Naive Bayes after applying a principal components decomposition to the predictor matrix. \code{type="orig"} uses the original predictor matrix.
#' 
#' @return A vector of predicted survival probabilities, with one entry corresponding to each row of \code{X.test}.
#' 
#' @examples
#' library(survival)
#' X <- flchain[,c("age","kappa","lambda")]
#' frac.train <- 0.7
#' ind.train <- sample.int(nrow(flchain),nrow(flchain)*0.7)
#' ind.test <- setdiff(1:nrow(flchain),ind.train)
#' X.train <- X[ind.train,]
#' X.test <- X[ind.test,]
#' T.train <- as.vector(flchain[ind.train,"futime"])
#' C.train <- as.vector(flchain[ind.train,"death"])
#' T.test <- as.vector(flchain[ind.test,"futime"])
#' C.test <- as.vector(flchain[ind.test,"death"])
#' SURVTIME <- 5000
#' 
#' survpreds.cnb <- CNB(X.train,T.train,C.train,5000,X.test=X.test)
#' eventpreds.cnb <- 1 - survpreds.cnb
#' 
#' cutpts <- c(0,0.1,0.2,0.3,1)
#' calib.table(list(eventpreds.cnb),T.test,C.test,cutpts,SURVTIME)
#' 
#' surv <- Surv(T.test,C.test)
#' survConcordance(surv ~ eventpreds.cnb)
#' ## Compare concordance to a Cox Model
#' survConcordance(surv ~ predict(coxph(Surv(T.train,C.train)~age+kappa+lambda,data=X.train),X.test))
CNB <- function(X,T,C,SURVTIME,X.test=NULL,scaleX=TRUE,type="PC") {

  if(type=="orig") {
    if(is.null(X.test)) { xnew <- X } else { xnew <- X.test }
    if(scaleX) { X <- scale(X); xnew <- scale(xnew)}
  } else {
    if(!scaleX) { warning("Scaling the design matrix is highly recommended when using the principal components version of CNB!") }
    X.pc <- prcomp(X,center=TRUE,scale.=scaleX)
    X <- xnew <- X.pc$x
    if(!is.null(X.test)) {
      xnew <- predict(X.pc,X.test)
    }
  }
  
  sf <- survfit(Surv(T,C)~1)  
  sf.C <- survfit(Surv(T,1-C)~1)
  
  KM <- stepfun(sf$time,c(1,sf$surv))
  KM.C <- stepfun(sf.C$time,c(1,sf.C$surv))
  
  KM.surv <- KM(SURVTIME)
    
  ## Compute inverse probability of censoring weights
  wts.lt <- (T < SURVTIME & C == 1) * 1/KM.C(pmin(T,rep(SURVTIME,length(T))))
  
  EX.lt <- colSums(X*wts.lt)/sum(wts.lt)
  EX2.lt <- colSums(X^2*wts.lt)/sum(wts.lt)
  sd.lt <- sqrt(EX2.lt - EX.lt^2)
  
  ind.gt <- (T>=SURVTIME)
  EX.gt <- colMeans(X[ind.gt,])
  EX2.gt <- colMeans(X[ind.gt,]^2)
  sd.gt <- sqrt(EX2.gt - EX.gt^2)
  
  num <- exp(colSums(log(dnorm(t(xnew),mean=EX.gt,sd=sd.gt))))*KM.surv
  denom <- num + exp(colSums(log(dnorm(t(xnew),mean=EX.lt,sd=sd.lt))))*(1-KM.surv)
  
  return(num/denom)
}
