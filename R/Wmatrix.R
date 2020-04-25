#'@title Build the Matrix of Environmental Covariables for Reaction Norm
#'
#'
#' @description Estimates a weather covariable realized matrix with dimensions of q x k, for q environments and k covariables.
#'
#' @author Germano Costa Neto
#'
#' @param weather.data data.frame of environmental variables gotten from \code{get_weather()}.
#' @param is.processed boolean. Indicates whether the dataframe was previously processed with \code{summaryWTH()} and contains means, medians, etc.
#' @param id.names character. Indicates the name of the columns to be used as id for the environmental variables to be analysed.
#' @param env.id   character. Indicates the name of the columns to be used as id for environments.
#' @param var.id  vector (character). Indicates which variables will be used in the analysis.
#' @param probs vector(numeric). Indicates the probability quantiles, as probs = {0,1}, if cardinals is NULL. If is quantiles=NULL, quantiles = c(0.01,.25,.50,.99).
#' @param by.interval boolean. Indicates if temporal intervals must be computed insied of each environment. Default = FALSE.
#' @param time.window vector (numeric). If \code{by.interval = TRUE}, this argument indicates the temporal breaks for delimited intervals.
#' @param names.window vector(character). If by.interval = TRUE, this argument indicates the names of the desirable intervals.
#' @param center boolean. Indicates whether the matrix should be centered.
#' @param scale boolean. If scale=TRUE, the variables assume a mean-centered scaled distribution, with \eqn{x~N(0,1)}.
#' @param sd.tol numeric. Default value equal to 10.
#' @param statistic vector (character). Indicates what statistic must be runned, statistic = c('all','sum','mean','quantile'). Default: 'mean'.
#' @param tol TODO.
#' @param QC boolean. Indicates with Quality Control is applied. QC is based on the standard deviation tolerance (sd.tol), removing variables (x) with sd(x) > sd.tol
#'
#' @return
#' Returns a weather covariable realized matrix with dimensions of q x k, for q environments and k covariables.
#'
#' @details
#' TODO
#'
#' @examples
#' ### Fetching weather information from NASA-POWER
#' weather.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA')
#'
#' ### Mean-centered and scaled matrix
#' W <- W.matrix(weather.data = weather.data, by.interval = FALSE)
#'
#' ## same as SummaryWTH, we can add time.windows
#' W <- W.matrix(weather.data = weather.data, by.interval = TRUE,
#'               time.window = c(0, 14, 35, 60, 90, 120))
#'
#' ## and select the statistic to be used
#' W <- W.matrix(weather.data = weather.data, by.interval = TRUE, statistic = 'mean',
#'               time.window = c(0, 14, 35, 60, 90, 120))
#'
#' W <- W.matrix(weather.data = weather.data, by.interval = TRUE, statistic = 'quantile',
#'               time.window = c(0, 14, 35, 60, 90, 120))
#'
#' ## We can perform a Quality Control (QC) based on the maximum sd tolerated
#' W <- W.matrix(weather.data = weather.data, by.interval = FALSE, QC = TRUE)
#'
#' ## We can perform a Quality Control (QC) based on the maximum sd tolerated
#' W <- W.matrix(weather.data = weather.data, by.interval = FALSE, QC = TRUE, sd.tol = 3)
#'
#' ## or create W for specific variables
#' id.var = c('T2M_MAX','T2M_MIN','T2M')
#' W <- W.matrix(weather.data = weather.data, var.id = id.var)
#'
#' ## or combine with summaryWTH by using is.processed = TRUE
#' data <- summaryWTH(weather.data, env.id = 'env', statistic = 'quantile')
#' W <- W.matrix(weather.data = data, is.processed = TRUE)
#'
#' @importFrom stats sd
#' @importFrom reshape2 acast
#'
#' @export

W.matrix = function(weather.data, is.processed=FALSE,id.names=NULL,env.id=NULL,var.id=NULL,
                    probs=NULL,by.interval=NULL,time.window=NULL,names.window=NULL,
                    center=TRUE,scale=TRUE, sd.tol = 10,statistic=NULL,
                    tol=1E-3, QC=FALSE){

  if(is.null(statistic)) statistic <-'mean'
  if(is.null(by.interval)) by.interval <- FALSE

  # if the env data are alredy processed in means, medians etc, ignore
  if(isFALSE(is.processed)){
    W<-summaryWTH(weather.data=weather.data,id.names=id.names,
                  env.id = env.id,
                  statistic=statistic,
                  probs = probs,var.id=var.id,
                  by.interval=by.interval,
                  time.window=time.window,
                  names.window=names.window)


    colid <- c('env','variable','interval')
    W <- melt(W,measure.vars = names(W)[!names(W)%in%colid],variable.name = 'stat' )
    W$variable <-paste(W$variable,W$stat,sep = '_')
    if(isTRUE(by.interval)) W  <- reshape2::acast(W,env~variable+interval,value.var = "value")
    if(isFALSE(by.interval)) W <- reshape2::acast(W,env~variable,value.var = "value")
  }

  # parametrization for w~N(0,1)
  W <- W.scale(weather.data = W,center = center,scale = scale,sd.tol = sd.tol,QC = QC)

  return(W)

}





W.scale <-function(weather.data, center=TRUE,scale=TRUE, sd.tol = 4,tol=1E-3,QC=FALSE){

  sdA   <- apply(weather.data,2,sd)
  t <- ncol(weather.data)
  removed <- names(sdA[sdA > sd.tol])
  if(!is.matrix(weather.data)){stop('weather.data must be a matrix')}
  weather.data<- scale(weather.data+tol,center = center,scale = scale)
  if(isTRUE(QC)) weather.data <- weather.data[,!colnames(weather.data) %in% removed]

  r <- length(removed)

  if(isTRUE(QC)){
    cat(paste0('------------------------------------------------','\n'))
    cat(paste0('Quality Control based on sd.tol=',sd.tol,'\n'))
    cat(paste0('Removed variables:','\n'))
    cat(paste0(r,' from ',t,'\n'))
    cat(paste0(removed,'\n'))
    cat(paste0('------------------------------------------------','\n'))
  }

  return(weather.data)
}
