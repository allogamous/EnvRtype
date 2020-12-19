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
#'
#' @importFrom stats quantile
#' @importFrom foreach %dopar% %:% foreach
#' @importFrom reshape2 acast
#'
#' @export

env_typing <- function(env.data,var.id,env.id,cardinals=NULL,days.id=NULL,
                      time.window=NULL,names.window=NULL,quantiles=NULL,
                      id.names=NULL,by.interval=FALSE,scale=FALSE,
                      format=NULL){

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

  if(isFALSE(by.interval)) .GET <-data.frame(.GET,interval='full-cycle')
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
