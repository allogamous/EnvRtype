#'@title  Support functions to estimate basic radiation parameters
#'
#' @description Core of functions to estimate atmospheric parameters related to solar radiation. Fore more details about the equations see de FAO-evapotranspiration publication http://www.fao.org/3/X0490E/x0490e07.htm#radiation.
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather()} output or A \code{get_weather()}-like dataframe.
#' @param merge boolean. If \code{TRUE}, calculated variables are merged to the original dataframe. Default is TRUE
#'
#' @return Returns a dataframe with parameters related to solar radiation. See details for further information.
#'
#' @details
#' This functions requires day of the year, latitude, thermal radiative flux, and sky insolation which are provided by \code{get_weather()}. If one of these parameters are missing, an error will be returned.
#' \itemize{
#'  \item n: Actual duration of sunshine (hour)
#'  \item N: Daylight hours (hour)
#'  \item RTA: Extraterrestrial radiation (MJ/m^2/day)
#'  \item SRAD: Solar radiation (MJ/m^2/day)
#'  }
#'
#' @examples
#'\dontrun{
#' ### Fetching weather information from NASA-POWER
#' env.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA')
#'
#' ### Calculating solar radiation
#' param_radiation(env.data)
#'
#' env.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA') %>%
#'            param_radiation()
#'
#' #' env.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA') %>%
#'            param_radiation(day.id='DOY',latitude='LAT')
#'}
#' @importFrom stats median
#'
#' @export

# http://www.fao.org/3/X0490E/x0490e07.htm#radiation
param_radiation <-function(env.data, day.id = 'DOY', latitude = 'LAT',solarRadiation=NULL, merge=TRUE){


  DOY <- env.data[,which(colnames(env.data) %in% day.id)]
  LAT <- env.data[,which(colnames(env.data) %in% latitude)]

  Ra <- function(J,lat){
    rlat = deg2rad(lat)
    fi = .409*sin((2*pi/365)*J-1.39)
    dr = 1+ .033*cos((2*pi/365)*J)
    ws = acos(-tan(rlat)*tan(fi))
    Ra = (1440/pi)*.0820*dr*(ws*sin(rlat)*sin(fi)+cos(rlat)*cos(fi)*sin(ws))
    N = (24/pi)*ws
    cat('---------------------------------------------------------------------- \n')
    cat('Using NASA POWER API outputs to compute:  \n')
    cat('Extraterrestrial radiation (RTA, MJ/m^2/day)  \n')
    cat('Daylight hours (N, hours) \n')

    return(data.frame(Ra=Ra,N=N))
  }


  RadN <-Ra(J = DOY,lat = LAT)

  if(is.null(solarRadiation)) solarRadiation = 'ALLSKY_SFC_SW_DWN'

  DWN <- env.data[,which(colnames(env.data) %in% solarRadiation)]
#  LWD <- env.data$ALLSKY_SFC_LW_DWN

 # LWD[LWD == -999] <- NA
 # DWN[DWN == -99] <- NA
 # LWD[is.na(LWD)] <- median(LWD,na.rm=TRUE)
  DWN[is.na(DWN)] <- median(DWN,na.rm=TRUE)

 # Srad <- LWD-DWN
  n <- RadN$N*(DWN/RadN$Ra)

  cat('Actual duration of sunshine (n, hours) \n')
 # cat('Solar Radiation (SRAD, MJ/m^2/day) \n') # wrong concept!!
  cat('---------------------------------------------------------------------- \n')
  cat('\n')

  #if(!isTRUE(merge)) return(data.frame(n=n,N=RadN$N,RTA=RadN$Ra, SRAD=Srad))
  #if(isTRUE(merge)) return(data.frame(env.data,data.frame(n=n,N=RadN$N,RTA=RadN$Ra, SRAD=Srad)))

  if(!isTRUE(merge)) return(data.frame(n=n,N=RadN$N,RTA=RadN$Ra))
  if(isTRUE(merge)) return(data.frame(env.data,data.frame(n=n,N=RadN$N,RTA=RadN$Ra)))

}







