#'@title  Support functions to estimate atmospheric parameters
#'
#' @description Core of functions to estimate atmospheric parameters related to evapotranspiration. Fore more details about the equations see de FAO-evapotranspiration publication http://www.fao.org/3/X0490E/x0490e07.htm#atmospheric%20parameters.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather()}-like output with additional solar radiation information. If a \code{get_weather()} is provided, no further argument is required but \code{G} and \code{alpha}.
#' @param PREC character. Indicates the column of precipitation.
#' @param Tdew character. Indicates the column of dew/frost.
#' @param Tmin character. Indicates the column of minimum temperature.
#' @param Tmax character. Indicates the column of maximum temperature.
#' @param RH character. Indicates the column of reative humidity.
#' @param Rad character. Indicates the column of solar radiation. This parameter can be calculated from \code{Param_Radiation()}.
#' @param G numeric. Flux of heat conducted into the ground. Default is 0.
#' @param alpha Alpha of Priestley & Taylor's (1972) equation. Default is 1.26, which fit data from most sources.
#' @param Alt numeric. Elevation above sea level (meters)
#' @param merge boolean. If \code{TRUE}, calculated variables are merged to the original \code{get_weather()} dataframe.
#'
#' @return
#' A dataframe with parameters related to evapotranspiration. See details for further information.
#'
#' @details
#' This function requires a dataframe with all parameters listed above. If any is missing, an error will be returned.
#' The calculated variables are:
#' \itemize{
#'  \item SPV: Slope of saturation vapour pressure curve (kPa.Celsius)
#'  \item VPD: Vapour pressure deficit (kPa)
#'  \item ETP: Potential Evapotranspiration (mm.day)
#'  \item PEPT: Deficit by Precipitation (mm.day)
#'  }
#'
#' @examples
#' ### Fetching weather information from NASA-POWER
#' env.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA')
#'
#' ### Calculating solar radiation
#' env.data = Param_Radiation(env.data, merge = TRUE)
#'
#' ### Calculating solar radiation
#' param_atmospheric(env.data)

#' @export

# http://www.fao.org/3/X0490E/x0490e07.htm#atmospheric%20parameters

param_atmospheric <- function(env.data, PREC=NULL, Tdew=NULL,
                              Tmin=NULL, Tmax=NULL, RH=NULL,
                              Rad=NULL, G=NULL, Alt=600, alpha=1.26,
                              merge=FALSE){


  teten <- function(Temp) return(.61078*exp((17.27*Temp)/(Temp+237.3)))
  psyco <- function(atm) return((((1.013*1E-3)*atm)/.622*2.45)*.665*atm*1E-3)
  AtmP <-  function(elevation=600) return(101.3*((293-.0065*elevation)/293)^5.26)
  slope.vapor <- function(Tmed) return(4098*(.6108*exp((17.27*Tmed)/(Tmed+237.3)))/(Tmed+237.2)^2)

  # Pristley-Taylor Equation
  EToPT<-function(alfa=1.26,Srad,G=NULL,slope,psyc){
    if(is.null(G)){G=0}
    W = slope/(slope+psyc)
    return(alfa*W*(Srad-G)*.408)}



  if(is.null(PREC)) PREC <-'PRECTOT'; PREC <- env.data[,PREC]
  if(is.null(Tdew)) Tdew <-'T2MDEW';Tdew <- env.data[,Tdew]
  if(is.null(Tmin)) Tmin <-'T2M_MIN';Tmin<- env.data[,Tmin]
  if(is.null(Tmax)) Tmax <-'T2M_MAX';Tmax <- env.data[,Tmax]
  if(is.null(Alt)){ Alt <-'ALT';Alt <- env.data[,Alt]}
  # if(is.null(Alt)) Alt <- 600
  #if(isFALSE(Alt  %in% names(env.data))){Alt <- 600; cat('Missing ALT value. We adopted 600m. Please use the Extract_GIS funciton to collect ALT from SRTM files \n')}

  if(is.null(RH))  RH <-'RH2M'; RH<- env.data[,RH]
  if(is.null(Rad)) Rad <- 'SRAD';Rad <- env.data[,Rad]


  Tmed <- (Tmin+Tmax)/2
  ATP <- AtmP(Alt)
  Es <- teten(Tdew)
  Eamin <- teten(Tmin)
  Eamax <- teten(Tmax)
  VPD <- ((Eamin +Eamax)/2)-Es
  #VPD <- Es*(1-RH/100)
  Psy <- psyco(atm = ATP)
  Slope <- 4098*(.6108*exp((17.27*Tmed)/(Tmed+237.3)))/(Tmed+237.2)^2
  ETo <- EToPT(alfa=alpha,Srad=Rad,G=G,slope=Slope,psyc=Psy)
  PETo <- PREC-ETo

  cat('---------------------------------------------------------------------- \n')
  cat('Slope of saturation vapour pressure curve (SPV, in kPa.Celsius) \n')
  cat('Vapour pressure deficit (VPD, in kPa) \n')
  cat('Potential Evapotranspiration (ETP, in mm.day) \n')
  cat('Deficit by Precipitation - ETP (PETP, in mm.day) \n')
  cat('---------------------------------------------------------------------- \n')
  cat('\n')

  if(isFALSE(merge)) return(data.frame(VPD=VPD,SPV=Slope,ETP=ETo,PETP=PETo))
  if(isTRUE(merge)) return(data.frame(env.data,data.frame(VPD=VPD,SPV=Slope,ETP=ETo,PETP=PETo)))


}







