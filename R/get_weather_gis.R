#'@title  Easily Collection of Worldwide Daily Weather Data
#'
#'
#' @description Imports daily-scale weather data from the NASA POWER GIS database
#' @author Germano Costa Neto
#' @param env.id vector (character or level). Indicates site/environment identification
#' @param lat vector (numeric). Indicates de latitude valuesort
#' @param lon vector (numeric). containing longitude values
#' @param variables.names vector of variables names. Should be "T2M","T2M_MAX","T2M_MIN","PRECTOT", "WS2M","RH2M","T2MDEW", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN", and/or "ALLSKY_TOA_SW_DWN"
#' @param start.day start point
#' @param end.day end point
#' @param dir.path output directorie
#' @param save bollean. If TRUE, save each environmental data.frame as .csv in dir.path
#' @param temporal.scale character. Default = 'DAILY'. See get_power() function in nasapower package for more details
#' @importFrom utils write.csv
#' @importFrom nasapower get_power
#' @importFrom plyr ldply
#' @importFrom utils install.packages

#----------------------------------------------------------------------------------------
# getting weather data from NASAPOWER GIS database
# adaptation from nansapower package's get_power function
#----------------------------------------------------------------------------------------

get_weather = function(env.id = NULL,lat   = NULL,lon   = NULL,
                           start.day = NULL,end.day = NULL,
                           variables.names = NULL, dir.path = NULL,
                           save=FALSE,asdataframe=TRUE,temporal.scale = 'DAILY'){
  cat('------------------------------------------------ \n')
  cat('ATTENTION: This function requires internet access \n')
  cat('------------------------------------------------  \n')

  # checking in inputs
#  if(is.null(env.id)){env.id <- paste0("env", seq_along(lat))} #creates a name for the enviroments (if null)
#  if(!(is.character(env.id) || is.factor(env.id)))
#     {stop("env.id should be a vector of characters (e.g. 'env1') or factors")} #checks if env.id is character of factors

  # if nasapower is not installed, the CRAN version will be installed
  if (!requireNamespace('nasapower', quietly = TRUE)) {install.packages("nasapower")} #talvez seja interessante instalar junto com o pacote
  if (!requireNamespace('plyr', quietly = TRUE)) {install.packages("plyr")}

  # if output dir.path is null, the current directorie folder will be used
  if(is.null(dir.path)){dir.path = getwd()}

  # if start.day is null, current day - 1000 days
  if(is.null(start.day))
    {
    start.day<- Sys.Date()-1000
    cat(paste0('start.day is NULL','\n'))
    cat(paste0('matched as ',start.day,'\n'))
    }

  # if end.day is null, start.day + 30
  if(is.null(end.day))
  {
    end.day<- start.day+30
    cat(paste0('end.day is NULL','\n'))
    cat(paste0('matched as ',end.day,'\n'))
  }

  # if variables is null, the default list will be used
  if(is.null(variables.names)){
    variables.names = c("T2M","T2M_MAX","T2M_MIN","PRECTOT",
                   "WS2M","RH2M","T2MDEW","ALLSKY_SFC_LW_DWN",
                   "ALLSKY_SFC_SW_DWN","ALLSKY_TOA_SW_DWN")
  }

  # preapring outputs
  env.id = as.factor(env.id)
  .Ne    = length(env.id)
  .C     = vector(length = .Ne,"list")

  for(.E in 1:.Ne){
    CL = data.frame(nasapower::get_power(community = "AG",lonlat = c(lon[.E], lat[.E]),
                                         pars = variables.names,
                                         dates = c(start.day[.E],end.day[.E]),
                                         temporal_average = temporal.scale))
    CL$daysFromStart = 1:nrow(CL)
    .C[[.E]] = CL
    names(.C)[[.E]] = env.id[.E]

    # if save is true, write the weather into csv files
    if(isTRUE(save)){write.csv(file=paste(env.id[.E],".csv",sep=""), x = CL)}

  }

  names(.C) = env.id
  # if asdataframe is true, df_convert function will be used
  if(isTRUE(asdataframe)) .C = plyr::ldply(.C)
  names(.C)[names(.C) %in% '.id'] = 'env'

  return(.C)
}


