#'@title  Basic Summary Statistics for Environmental Data
#'
#'
#' @description Summarize get_weather() outputs based on environments and defined time intervals (e.g.,phenology)
#' @author Germano Costa Neto
#' @param id.names vector (character). Indicates the name of the columns to be used as id for the environmental variables to be analysed.
#' @param env.id   vector (character). Indicates the name of the columns to be used as id for environments.
#' @param var.id  character. Indicates which variables will be used in the analysis.
#' @param statistic vector (character). Indicates what statistic must be runned, statistic = c('all','sum','mean','quantile'). Default: 'all'.
#' @param probs vector(numeric). Indicates the probability quantiles, as probs = {0,1}. If is NULL, probs = c(0.25,.50,.75).
#' @param by.interval boolean. Indicates if temporal intervals must be computed insied of each environment. Default = FALSE.
#' @param time.window vector (numeric). If by.interval = TRUE, this argument indicates the temporal breaks for delimited intervals.
#' @param names.window vector(character). If by.interval = TRUE, this argument indicates the names of the desirable intervals.
#' @importFrom utils write.csv
#' @importFrom nasapower get_power
#' @importFrom utils install.packages
#' @importFrom reshape2 melt dcast
#' @importFrom plyr ddply . summarise


summaryWTH <- function(x,id.names=NULL,
                       env.id = NULL,
                       days.id = NULL,
                       var.id=NULL,
                            statistic=NULL,
                            probs = NULL,
                            by.interval=FALSE,
                            time.window=NULL,
                            names.window=NULL){
  if (!require(reshape2)) install.packages("reshape2"); require(reshape2)
  if (!require(plyr)) install.packages("plyr") ; require(plyr)
  class(x) <- 'data.frame'
  if(!any(statistic %in% c('all','sum','mean','quantile'))) statistic <- 'all'

  .GET <- meltWTH(.GeTw = x,days = days.id,
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



SumEweather <- function(.GetW,by.interval=FALSE){

  env <- variable  <- value <- interval <- NULL #supressor

  if(isFALSE(by.interval)){
    return(ddply(.GetW,.(env,variable),summarise,sum=sum(value,na.rm=T)))
  }
  if(isTRUE(by.interval)){

    return(ddply(.GetW,.(env,interval,variable),summarise,sum=sum(value,na.rm=T)))
  }

}

MeanEweather <- function(.GetW,by.interval=FALSE){

  env <- variable <- summarise <- value <- interval <- NULL

  if(isFALSE(by.interval)){
    return(ddply(.GetW,.(env,variable),summarise,mean=mean(value,na.rm=T)))
  }
  if(isTRUE(by.interval)){

    return(ddply(.GetW,.(env,interval,variable),summarise,mean=mean(value,na.rm=T)))
  }

}





QuantEweather <- function(.GetW,prob=c(.25,.5,.75),by.interval=FALSE){

  e <- v <- s <- i <- NULL #supressor

  require(doParallel)
  require(reshape2)
  (envs <-unique(.GetW$env))
  (vars <- unique(.GetW$variable))

  if(isFALSE(by.interval)){

    Q.t <- foreach(e=1:length(envs), .combine = "rbind") %:%
      foreach(v=1:length(vars), .combine = "rbind") %:%
      foreach(s=1:length(prob), .combine = "rbind") %dopar% {


        quantil <- data.frame(qt=quantile(x=.GetW$value[.GetW$env %in% envs[e] & .GetW$variable %in% vars[v]],prob[s]),
                              prob=paste0('prob_',prob[s]),
                              env=envs[e],var=vars[v])
        return(quantil)
      }
    Q.t <- dcast(Q.t ,env+var~prob,value.var='qt')
    names(Q.t)[1:2] <- c('env','variable')
    return(Q.t)
  }
  if(isTRUE(by.interval)){

    (inters <- unique(.GetW$interval))

    Q.t <- foreach(e=1:length(envs),    .combine = "rbind") %:%
      foreach(i=1:length(inters),  .combine = "rbind") %:%
      foreach(v=1:length(vars) ,   .combine = "rbind") %:%
      foreach(s=1:length(prob),    .combine = "rbind") %dopar% {


        quantil <- data.frame(qt=quantile(x=.GetW$value[.GetW$env %in% envs[e] & .GetW$interval %in% inters[i] & .GetW$variable %in% vars[v]],prob[s]),
                              prob=paste0('prob_',prob[s]),
                              env=envs[e],var=vars[v],
                              interval=inters[i])
        return(quantil)
      }
    Q.t <- dcast(Q.t ,env+var+interval~prob,value.var='qt')
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

    ts <- suppressWarnings(melt(.GeTw,measure.vars=var.id,id.vars = id.names))
    ts$value <- suppressWarnings(as.numeric(ts$value))
    return(ts)
  }
  if(isTRUE(by.interval)){

    .GeTw$interval <- stage.by.dae(.dae = .GeTw[,'daysFromStart'],
                                   .breaks = time.window,
                                   .names = names.window)

    id.names <-c(id.names,'interval')
    ts <- suppressWarnings(melt(.GeTw,measure.vars=var.id,id.vars = id.names))
    ts$value <- suppressWarnings(as.numeric(ts$value))
    return(ts)
  }
}

stage.by.dae = function(.dae=NULL, .breaks=NULL, .names=NULL){
  if(is.null(.dae)) stop(".dae is missing")
  if(is.null(.breaks)) .breaks <-seq(from=1-min(.dae),to=max(.dae)+10,by=10)
  if(is.null(.names))  .names  <-paste0("Interval_",.breaks)
  .breaks <- c(.breaks,Inf)
  pstage = cut(x = .dae,breaks=.breaks,right = F)
  levels(pstage) = .names
  return(pstage)
}
