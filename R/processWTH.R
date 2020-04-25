#'@title  Basic Processing Tools for Enriching get_weather outputs
#'
#'
#' @description A wraper of \code{Param_Radiation()}, \code{Param_Atmospheric()}, and \code{Param_Temperature()} for evaluating get_weather() datasets. Calculates a series of parameters based on the \code{get_weather()} object.
#' @author Germano Costa Neto
#'
#' @param weather.data data.frame. A \code{get_weather()} output.
#'
#' @return
#' Returns a \code{get_wheather()} object with an additional set of parameters calculated from the nasapower data.
#'
#' @details
#' This function requires a dataframe with all parameters listed above. If any is missing, an error will be returned.
#' The estimated variables are:
#' \itemize{
#'  \item n: Actual duration of sunshine (hour)
#'  \item N: Daylight hours (hour)
#'  \item RTA: Extraterrestrial radiation (MJ/m^2/day)
#'  \item SRAD: Solar radiation (MJ/m^2/day)
#'  \item SPV: Slope of saturation vapour pressure curve (kPa.Celsius)
#'  \item VPD: Vapour pressure deficit (kPa)
#'  \item ETP: Potential Evapotranspiration (mm.day)
#'  \item PEPT: Deficit by Precipitation (mm.day)
#'  \item GDD: Growing Degree Day (oC/day)
#'  \item FRUE: Effect of temperature on radiation use efficiency (from 0 to 1)
#'  \item T2M_RANGE: Daily Temperature Range (oC day)
#' }
#'
#' @examples
#' ### Fetching weather information from NASA-POWER
#' weather.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA')
#'
#' ### Returning calculated parameters merged to the \code{get_weather()} dataframe
#' processWTH(weather.data)
#'
#' @export

processWTH <- function(weather.data){

  # computing radiation paramters
  weather.data <-Param_Radiation(weather.data=weather.data, merge = TRUE)
  # computing atmospheric paramters
  weather.data <-Param_Atmospheric(weather.data=weather.data, merge = TRUE)
  weather.data <- Param_Temperature(weather.data=weather.data,merge = TRUE)

  return(weather.data)
}


