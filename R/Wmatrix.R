#'@title Build the Matrix of Environmental Covariables for Reaction Norm
#'
#'
#' @description Returns a covariable realzed matrix with dimensions of q x k, for q environments and k covariables.
#' @author Germano Costa Neto
#' @param df.cov data.frame of environmental variables
#' @param id.names vector (character). Indicates the name of the columns to be used as id for the environmental variables to be analysed.
#' @param env.id   vector (character). Indicates the name of the columns to be used as id for environments.
#' @param var.id  character. Indicates which variables will be used in the analysis.
#' @param statistic vector (character). Indicates what statistic must be runned, statistic = c('all','sum','mean','quantile'). Default: 'mean'.
#' @param probs vector(numeric). Indicates the probability quantiles, as probs = {0,1}, if cardinals is NULL. If is quantiles=NULL, quantiles = c(0.01,.25,.50,.99).
#' @param by.interval boolean. Indicates if temporal intervals must be computed insied of each environment. Default = FALSE.
#' @param time.window vector (numeric). If by.interval = TRUE, this argument indicates the temporal breaks for delimited intervals.
#' @param names.window vector(character). If by.interval = TRUE, this argument indicates the names of the desirable intervals.
#' @param scale boolean. If scale=TRUE, the variables (x) assumes a mean-centered scaled distribution, with x~N(0,1).
#' @param QT boolean. Indicates with Quality Control is applied. QC is based on the standard deviation tolerance (sd.tol), removing variables (x) with sd(x) > sd.tol
#' @param sd.tol numeric. Default value equal to 10.
#' @importFrom stats sd



W.matrix = function(df.cov, is.processed=FALSE,id.names=NULL,env.id=NULL,var.id=NULL,
                    probs=NULL,by.interval=NULL,time.window=NULL,names.window=NULL,
                    center=T,scale=T, sd.tol = 10,statistic=NULL,
                    tol=1E-3,QC=FALSE){

  require(reshape2)

  if(is.null(statistic)) statistic <-'mean'
  if(is.null(by.interval)) by.interval <- FALSE

  # if the env data are alredy processed in means, medians etc, ignore
  if(isFALSE(is.processed)){
    W<-summaryWTH(x=df.cov,id.names=id.names,
                  env.id = env.id,
                  statistic=statistic,
                  probs = probs,var.id=var.id,
                  by.interval=by.interval,
                  time.window=time.window,
                  names.window=names.window)


    colid <- c('env','variable','interval')
    W <- melt(W,measure.vars = names(W)[!names(W)%in%colid],variable.name = 'stat' )
    W$variable <-paste(W$variable,W$stat,sep = '_')
    if(isTRUE(by.interval)) W  <-acast(W,env~variable+interval,value.var = "value")
    if(isFALSE(by.interval)) W <-acast(W,env~variable,value.var = "value")
  }

  # parametrization for w~N(0,1)
  W <- W.scale(df.cov = W,center = center,scale = scale,sd.tol = sd.tol,QC = QC)

  return(W)

}





W.scale <-function(df.cov, center=T,scale=T, sd.tol = 4,tol=1E-3,QC=F){

  sdA   <- apply(df.cov,2,sd)
  t <- ncol(df.cov)
  removed <- names(sdA[sdA > sd.tol])
  if(!is.matrix(df.cov)){stop('df.cov must be a matrix')}
  df.cov<- scale(df.cov+tol,center = center,scale = scale)
  if(isTRUE(QC)) df.cov <- df.cov[,!colnames(df.cov) %in% removed]

  r <- length(removed)

  if(isTRUE(QC)){
    cat(paste0('------------------------------------------------','\n'))
    cat(paste0('Quality Control based on sd.tol=',sd.tol,'\n'))
    cat(paste0('Removed variables:','\n'))
    cat(paste0(r,' from ',t,'\n'))
    cat(paste0(removed,'\n'))
    cat(paste0('------------------------------------------------','\n'))
  }

  return(df.cov)
}
