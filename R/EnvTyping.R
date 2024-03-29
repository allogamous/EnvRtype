#'@title Environmental Typologies based on Cardinal or Quantilic Limits
#'
#'
#' @description Returns environmental typologies that can be used as envirotype markers. This typologies are given based on cardinals (discrete intervals for each variable). The user must inform the name of the environmental variables of interest and inform the cardinal option as list of cardinal vectors for each variable. Its also possible to fix intervals or scale data.
#' @author Germano Costa Neto
#'
#' @param env.data data.frame of environmental variables gotten from \code{get_weather()}.
#' @param id.names vector (character). Indicates the name of the columns to be used as id for the environmental variables to be analysed.
#' @param env.id   vector (character). Indicates the name of the columns to be used as id for environments.
#' @param var.id  character. Indicates which variables will be used in the analysis.
#' @param cardinals list (numeric). A list of cardinals (vector of numeric thresholds) for each variable. Indicates the cardinal limtis for each environmental type. If is NULL, see quantles argument.
#' @param days.id character. Name of the columns indicating the days from start.
#' @param quantiles vector (numeric). Indicates the probability quantiles, as probs = {0,1}, if cardinals is \code{NULL}. If is \code{quantiles=NULL}, \code{quantiles = c(0.01,.25,.50,.99)}.
#' @param by.interval boolean. Indicates if temporal intervals must be computed insied of each environment. Default = FALSE.
#' @param time.window vector (numeric). If \code{by.interval = TRUE}, this argument indicates the temporal breaks for delimited intervals.
#' @param names.window vector(character). If \code{by.interval = TRUE}, this argument indicates the names of the desirable intervals.
#' @param scale boolean. If \code{scale=TRUE}, the variables (x) assumes a mean-centered scaled distribution, with x~N(0,1).
#' @param format character. The shape of the output, assuming \code{format = c('long','wide')}. Default is 'long'.
#'
#' @return
#' A dataframe with environmental typologies. Event frequencies are provided for each variable at each environment within each interval.
#'
#' @details
#' TODO
#'
#' @examples
#'\dontrun{
#' ### Fetching weather information from NASA-POWER
#' env.data = get_weather(lat = -13.05, lon = -56.05)
#'
#' ### By.intervals (generic time intervals)
#' env_typing(env.data = env.data, env.id = 'env', var.id = 'T2M', by.interval = TRUE)
#'
#' ### By.intervals (specific time intervals)
#' env_typing(env.data = env.data,
#'           env.id = 'env',
#'           var.id = 'T2M',
#'           by.interval = TRUE,
#'           time.window = c(0, 15, 35, 65, 90, 120))
#'
#' ### By.intervals (specific time intervals and with specific names
#' env_typing(env.data = env.data,
#'                  env.id = 'env',
#'                  var.id = 'T2M',
#'                  by.interval = TRUE,
#'                  time.window = c(0, 15, 35, 65, 90, 120),
#'                  names.window = c('1-intial growing',
#'                                   '2-leaf expansion I',
#'                                   '3-leaf expansion II',
#'                                   '4-flowering',
#'                                   '5-grain filling',
#'                                   '6-maturation'))
#'
#' ### With set cardinals.
#' env_typing(env.data = env.data,
#'           var.id = c('T2M','PRECTOT','WS2M'),
#'           cardinals = list(T2M = c(0, 9, 22, 32, 45),
#'                            PRECTOT = c(0, 5, 10),
#'                            WS2M = NULL),
#'           env.id = 'env')
#'}
#' @importFrom stats quantile
#' @importFrom foreach %dopar% %:% foreach
#' @importFrom reshape2 acast
#'
#' @export

env_typing <- function(env.data,var.id,env.id,cardinals=NULL,days.id=NULL,
                      time.window=NULL,names.window=NULL,quantiles=NULL,
                      id.names=NULL,by.interval=FALSE,scale=FALSE,
                      format=NULL){

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

      Q.t <-
        foreach::foreach(e=1:length(envs),    .combine = "rbind") %:%
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

  x <- j <- s <- median <-  NULL #supressor

  #creating local functions based on '%:%' and '%dopar%'
  '%:%' <- foreach::'%:%'
  '%dopar%' <- foreach::'%dopar%'

  if(isTRUE(scale)) env.data[,var.id] <- scale(env.data[,var.id],center = TRUE,scale = TRUE)

  .GET <- meltWTH(.GeTw = env.data,days = days.id,
                  by.interval = by.interval,
                  time.window = time.window,
                  names.window = names.window,
                  id.names =  id.names,var.id = var.id,
                  env.id=env.id)

  if(isFALSE(by.interval)) .GET <-data.frame(.GET,interval='by environment')
  env <- as.character(unique(.GET$env))
  int <- as.character(unique(.GET$interval))
  vars <- as.character(unique(.GET$variable))

  if(is.null(quantiles)) quantiles <- c(.01,.25,.50,.75,.99)

  if(is.null(cardinals)){
    cardinals <- vector('list',length = length(vars))
    names(cardinals) = vars
  }
  # testing for missing cardinals (if null, c(min, 1st qu, 2 qt, 3qt, max))
  for(i in 1:length(vars)){
    if(is.null(cardinals[[i]])){
      SM <- round(as.numeric(quantile(.GET$value[.GET$variable %in% vars[i]],quantiles)),3)
      SM <- SM[!duplicated(SM)]
      #.NULL <- which(SM == 0) # there is not allowed 2 equal breaks
      #if(length(.NULL) > 1) SM <- SM[-sample(.NULL,length(.NULL) -1,replace = T)]
      cardinals[[i]] <- SM
    }
  }


  results <- foreach::foreach(s=1:length(vars), .combine = "rbind") %:% foreach::foreach(j=1:length(env), .combine = "rbind") %:% foreach::foreach(t=1:length(int), .combine = "rbind") %dopar% {

    .out=data.frame(table(cut(x =  .GET$value[which(.GET$env == env[j] & .GET$variable == vars[s] & .GET$interval == int[t])],
                              breaks = cardinals[[s]],right = TRUE)),env=env[j],interval=int[t],var=vars[s])

    return(.out)
  }
  names(results)[1] <- 'env.variable'
  results$env.variable <- paste0(results$var,'_',results$env.variable)
  if(isTRUE(by.interval))  results$env.variable <- paste0(results$env.variable,'_',results$interval)
  if(is.null(format)) format <- 'long'
  #if(!any(format %in% c('wide','long'))) cat('format as long'); format <- 'long'
  if(format == 'wide') results <- reshape2::acast(results,env~env.variable,median,value.var = 'Freq'); results[is.na(results)] <-0
  return(results)
}
