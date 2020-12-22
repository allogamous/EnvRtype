#'@title  Basic Summary Statistics for Environmental Data
#'
#'
#' @description Summarize \code{get_weather()} outputs based on environments and defined time intervals (e.g.,phenology)
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather()} output.
#' @param id.names vector (character). Indicates the name of the columns to be used as id for the environmental variables to be analysed.
#' @param env.id   character. Name of the columns to be used as identification for environments.
#' @param days.id character. Name of the columns indicating the days from start.
#' @param var.id  character. Indicates which variables will be used in the analysis.
#' @param statistic vector (character). Indicates what statistic must be analysed, \code{statistic = c('all','sum','mean','quantile')}. Default: 'all'.
#' @param probs vector(numeric). Indicates the probability quantiles, as probs = {0,1}. If is NULL, \code{probs = c(0.25,.50,.75)}.
#' @param by.interval boolean. Indicates if temporal intervals must be computed within each environment. Default = FALSE.
#' @param time.window vector (numeric). If \code{by.interval = TRUE}, this argument indicates the temporal breaks for delimited intervals.
#' @param names.window vector(character). If \code{by.interval = TRUE}, this argument indicates the names of the desirable intervals.
#'
#' @details
#' TODO
#'
#' @examples
#' ### Fetching weather information from NASA-POWER
#' env.data = get_weather(lat = -13.05, lon = -56.05)
#'
#' ### Basic summary
#' summaryWTH(env.data)
#'
#' ### Returning only mean values
#' summaryWTH(env.data, env.id = 'env', statistic = 'mean')
#'
#' ### Summary by time intervals given by time.window and names.window
#' summaryWTH(env.data, env.id = 'env', by.interval = TRUE,
#'            time.window = c(0, 14, 35, 60, 90, 120),
#'            names.window = c('P-E', 'E-V1', 'V1-V4', 'V4-VT', 'VT-GF', 'GF-PM'))
#'
#'
#' @importFrom utils install.packages
#' @importFrom reshape2 melt dcast
#' @importFrom plyr ddply . summarise
#' @importFrom foreach %dopar% %:% foreach
#'
#' @export


summaryWTH <- function(env.data,id.names=NULL,
                       env.id = NULL,
                       days.id = NULL,
                       var.id=NULL,
                            statistic=NULL,
                            probs = NULL,
                            by.interval=FALSE,
                            time.window=NULL,
                            names.window=NULL){


  SumEweather <- function(.GetW,by.interval=FALSE){

    env <- variable  <- value <- interval <- NULL #supressor

    if(isFALSE(by.interval)){
      return(plyr::ddply(.GetW,plyr::.(env,variable),plyr::summarise,sum=sum(value,na.rm=TRUE)))
    }
    if(isTRUE(by.interval)){

      return(plyr::ddply(.GetW,plyr::.(env,interval,variable),plyr::summarise,sum=sum(value,na.rm=TRUE)))
    }

  }

  MeanEweather <- function(.GetW,by.interval=FALSE){

    env <- variable <- summarise <- value <- interval <- NULL

    if(isFALSE(by.interval)){
      return(plyr::ddply(.GetW,plyr::.(env,variable),plyr::summarise,mean=mean(value,na.rm=TRUE)))
    }
    if(isTRUE(by.interval)){

      return(plyr::ddply(.GetW,plyr::.(env,interval,variable),plyr::summarise,mean=mean(value,na.rm=TRUE)))
    }

  }

  QuantEweather <- function(.GetW,prob=c(.25,.5,.75),by.interval=FALSE){

    e <- v <- s <- i <- NULL #supressor

    #creating local functions based on '%:%' and '%dopar%'
    '%:%' <- foreach::'%:%'
    '%dopar%' <- foreach::'%dopar%'

    (envs <-unique(.GetW$env))
    (vars <- unique(.GetW$variable))

    if(isFALSE(by.interval)){

      Q.t <- foreach::foreach(e=1:length(envs), .combine = "rbind") %:%
        foreach::foreach(v=1:length(vars), .combine = "rbind") %:%
        foreach::foreach(s=1:length(prob), .combine = "rbind") %dopar% {


          quantil <- data.frame(qt=quantile(x=.GetW$value[.GetW$env %in% envs[e] & .GetW$variable %in% vars[v]],prob[s]),
                                prob=paste0('prob_',prob[s]),
                                env=envs[e],var=vars[v])
          return(quantil)
        }
      Q.t <- reshape2::dcast(Q.t ,env+var~prob,value.var='qt')
      names(Q.t)[1:2] <- c('env','variable')
      return(Q.t)
    }
    if(isTRUE(by.interval)){

      (inters <- unique(.GetW$interval))

      Q.t <- foreach::foreach(e=1:length(envs),    .combine = "rbind") %:%
        foreach::foreach(i=1:length(inters),  .combine = "rbind") %:%
        foreach::foreach(v=1:length(vars) ,   .combine = "rbind") %:%
        foreach::foreach(s=1:length(prob),    .combine = "rbind") %dopar% {


          quantil <- data.frame(qt=quantile(x=.GetW$value[.GetW$env %in% envs[e] & .GetW$interval %in% inters[i] & .GetW$variable %in% vars[v]],prob[s]),
                                prob=paste0('prob_',prob[s]),
                                env=envs[e],var=vars[v],
                                interval=inters[i])
          return(quantil)
        }
      Q.t <- reshape2::dcast(Q.t ,env+var+interval~prob,value.var='qt')
      names(Q.t)[1:3] <- c('env','variable','interval')
      return(Q.t)
    }

  }

  meltWTH <- function(.GeTw,var.id=NULL, by.interval=FALSE,days=NULL,time.window=NULL,names.window=NULL,id.names=NULL,env.id=NULL){

    if(is.null(env.id)) env.id <- '.id'
    if(is.null(days)) days <- "daysFromStart"

    names(.GeTw)[names(.GeTw) %in% env.id] <- 'env'
    names(.GeTw)[names(.GeTw) %in% days  ] <- "daysFromStart"
    if(is.null(id.names)) id.names <- c('env','LON','LAT','YEAR','MM','DD','DOY','YYYYMMDD',"daysFromStart")
    if(is.null(var.id)) var.id <- names(.GeTw)[!names(.GeTw)%in%id.names]
    if(isFALSE(by.interval)){

      ts <- suppressWarnings(reshape2::melt(.GeTw,measure.vars=var.id,id.vars = id.names))
      ts$value <- suppressWarnings(as.numeric(ts$value))
      return(ts)
    }
    if(isTRUE(by.interval)){

      .GeTw$interval <- stage.by.dae(.dae = .GeTw[,'daysFromStart'],
                                     .breaks = time.window,
                                     .names = names.window)

      id.names <-c(id.names,'interval')
      ts <- suppressWarnings(reshape2::melt(.GeTw,measure.vars=var.id,id.vars = id.names))
      ts$value <- suppressWarnings(as.numeric(ts$value))
      return(ts)
    }
  }

  stage.by.dae = function(.dae=NULL, .breaks=NULL, .names=NULL){
    if(is.null(.dae)) stop(".dae is missing")
    if(is.null(.breaks)) .breaks <-seq(from=1-min(.dae),to=max(.dae)+10,by=10)
    if(is.null(.names))  .names  <-paste0("Interval_",.breaks)
    .breaks <- c(.breaks,Inf)
    pstage = cut(x = .dae,breaks=.breaks,right = FALSE)
    levels(pstage) = .names
    return(pstage)
  }

  if (!requireNamespace('reshape2')) utils::install.packages("reshape2")
  if (!requireNamespace('plyr')) utils::install.packages("plyr")
  class(env.data) <- 'data.frame'
  if(!any(statistic %in% c('all','sum','mean','quantile'))) statistic <- 'all'

  .GET <- meltWTH(.GeTw = env.data,days = days.id,
                     by.interval = by.interval,
                     time.window = time.window,
                     names.window = names.window,
                     id.names =  id.names,var.id = var.id,
                     env.id=env.id)

 # .GET <- droplevels(.GET)
  if(is.null(probs)) probs<-c(0.25,0.50,0.75)
  if(is.null(statistic)) statistic <- 'all'
  if(statistic == 'quantile') return(QuantEweather(.GetW = .GET,prob = probs,by.interval = by.interval))
  if(statistic == 'mean'    ) return(MeanEweather (.GetW = .GET,by.interval = by.interval))
  if(statistic == 'sum'    ) return(SumEweather  (.GetW = .GET,by.interval = by.interval))
  if(statistic == 'all'){

    id<-c('env','variable')

    if(isTRUE(by.interval)) id<-c(id,'interval')
    out<-merge(
      merge(
        MeanEweather (.GetW = .GET,by.interval = by.interval),
        SumEweather  (.GetW = .GET,by.interval = by.interval),
        by=id),
      QuantEweather(.GetW = .GET,prob = probs,by.interval = by.interval),
      by=id)
    suppressWarnings(return(out))
  }
}


