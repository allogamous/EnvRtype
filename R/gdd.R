#'@title Support functions to estimate thermal parameters
#'
#'
#' @description Calculation Of Heat According To The Growing Degree Day Model. From maximum (tmax) and minimum (tmin) temperature values, and according to the cardinal values of base temperature (Tbase) computes the daily thermal sum in ° C / day
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather()} output or A \code{get_weather()}-like dataframe.
#' @param Tmax character. Indicates the column of maximum air temperature (Celsius).
#' @param Tmin character. Indicates the column of minimum air temperature (Celsius).
#' @param Tbase1 numeric. Minimum cardinal value for temperature base for phenological development (Celsius).
#' @param Tbase2 numeric. Maximum cardinal value for temperature base for phenological development (Celsius).
#' @param Topt1 numeric. Lower temperature bound for phenological development (Celsius).
#' @param Topt2 numeric. Upper temperature bound for phenological development (Celsius).
#' @param merge boolean. If \code{TRUE}, calculated variables are merged to the original \code{get_weather()} dataframe.
#'
#' @return
#' A dataframe with parameters related to atmospheric temperature. See details for further information.
#'
#' @details
#' This function requires a dataframe with all parameters listed above. If any is missing, an error will be returned.
#' The calculated variables are:
#' \itemize{
#'  \item GDD: Growing Degree Day (oC/day)
#'  \item FRUE: Effect of temperature on radiation use efficiency (from 0 to 1)
#'  \item T2M_RANGE: Daily Temperature Range (oC day)
#'  }
#'
#' @examples
#' ### Fetching weather information from NASA-POWER
#' env.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA')
#'
#' ### Calculating solar radiation
#' Param_Temperature(env.data)
#'
#' @export

Param_Temperature <- function(env.data,Tmax=NULL, Tmin=NULL,
                              Tbase1=9,Tbase2=45, Topt1=26,Topt2=32,merge=FALSE) {

  if(is.null(Tmin)) Tmin <-'T2M_MIN';Tmin<- env.data[,Tmin]
  if(is.null(Tmax)) Tmax <-'T2M_MAX';Tmax <- env.data[,Tmax]

  Tmed <- (Tmin+Tmax)/2


  adjust_for_Tbase <- function(x, Tbase1) {
    ifelse(test = x < Tbase1, yes = Tbase1, no = x)
  }
  adjust_for_Tbase2 <- function(x, Tbase2) {
    ifelse(test = x > Tbase2, yes = Tbase2, no = x)
  }

  F.RUE.Temperature = function(Tbase1,Tbase2,Tmed,Topt1,Topt2){

    n = length(Tmed)
    FTAR = rep(1,n)
    for(i in 1:n){
      if(Tmed[i] < Topt1){FTAR[i] = (Tmed[i]-Tbase1)/(Topt1-Tbase1)}
      if(Tmed[i] > Topt2){FTAR[i] = (Tmed[i]-Topt2)/(Tbase2-Topt2)}
    }
    return(FTAR)
  }

  FRUE <- F.RUE.Temperature(Tbase1=Tbase1,Tmed=Tmed,Tbase2=Tbase2,Topt1=Topt1,Topt2=Topt2)

  tmax_adjusted <- adjust_for_Tbase(Tmax, Tbase1)
  tmin_adjusted <- adjust_for_Tbase(Tmin, Tbase1)

  tmax_adjusted <- adjust_for_Tbase2(tmax_adjusted, Tbase2)
  tmin_adjusted <- adjust_for_Tbase2(tmin_adjusted, Tbase2)

  gdd_temp <- (tmax_adjusted + tmin_adjusted) / 2 - Tbase1

  cat('---------------------------------------------------------------------- \n')
  cat('Growing Degree Day (GDD,oC/day) \n')
  cat('Effect of temperature on radiation use efficiency (FRUE, from 0 to 1)  \n')
  cat('Daily Temperature Range, (T2M_RANGE, oC day)\n')
  cat('---------------------------------------------------------------------- \n')

  if(isFALSE(merge))   return(data.frame(GDD=gdd_temp,FRUE = FRUE, T2M_RANGE=Tmax-Tmin ))
  if(isTRUE(merge))   return(data.frame(env.data,GDD=gdd_temp,FRUE = FRUE, T2M_RANGE=Tmax-Tmin ))

}
